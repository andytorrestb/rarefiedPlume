#!/bin/bash
#
# Submit the whole Kn sweep unattended, within VEGA's 4-queued-jobs-per-user
# limit. Runs on the LOGIN node in the background; it only sleeps and calls
# showq/msub every POLL seconds (default 1 h; negligible CPU). Start it from here:
#
#     cd cases/cai2012-unigas
#     nohup setsid bash runSweep.sh > sweep/sweep.out 2>&1 < /dev/null &
#
#     tail -f sweep/sweep.log        # progress
#     cat sweep/*.state              # per-case state
#     kill "$(cat sweep/sweep.pid)"  # stop submitting (queued jobs stay queued)
#
# Order:
#   1. runMeshAll.sh (generate, validate, mesh all seven); abort if it fails
#   2. runCase.sh per case, in CASES order, as queue slots free up
#      (Kn0p1 on normalq 24 h, the rest on the script's longq 120 h)
#   3. when a run job leaves the queue: finished -> runPost.sh (Cai comparison,
#      normalq, 32 cores); not finished (walltime) -> resubmit runCase.sh,
#      at most MAX_TRIES runs in total; then the case is marked failed
#   4. exit when every case is done or failed
#
# Safe to restart: state lives in sweep/; a re-run picks up where it stopped
# (delete sweep/ to start over). Never touches a job it did not submit.
set -o pipefail   # not -u: hpc_env.sh (Lmod, OpenFOAM) refuses it
cd "$(dirname "$0")" || exit 1
FAMILY=$PWD
USER_NAME=$(id -un)
POLL=${POLL:-3600}   # seconds between queue checks
MAX_QUEUED=${MAX_QUEUED:-4}
MAX_TRIES=${MAX_TRIES:-15}   # Kn0p0001 needs ~10 x 120 h runs (2026-11-24 estimate)
GRACE=180   # seconds after msub before a job missing from showq counts as ended
CASES=(Kn0p1_np1x Kn1_np1x Kn10_np1x Kn100_np1x Kn0p01_np1x Kn0p001_np1x Kn0p0001_np1x)

STATE=$FAMILY/sweep
mkdir -p "$STATE"
LOG=$STATE/sweep.log
log() { printf '%s  %s\n' "$(date '+%F %T')" "$*" >> "$LOG"; }

if [ -f "$STATE/sweep.pid" ] && kill -0 "$(cat "$STATE/sweep.pid")" 2>/dev/null; then
    echo "runSweep: already running as PID $(cat "$STATE/sweep.pid")" >&2; exit 1
fi
echo $$ > "$STATE/sweep.pid"
trap 'log "submitter stopped (signal)"; rm -f "$STATE/sweep.pid"; exit 1' INT TERM

# The environment, once (python + venv for the 'finished' check; msub/showq).
source ./hpc_env.sh >> "$STATE/env.log" 2>&1 || { log "FAILED: hpc_env.sh"; exit 1; }

# ---- helpers ---------------------------------------------------------------
# state file per job key: "phase jobid tries submitted_epoch"
get()  { [ -f "$STATE/$1.state" ] && cat "$STATE/$1.state" || echo "pending - 0 0"; }
put()  { echo "$2 $3 $4 $5" > "$STATE/$1.state"; }

queued_ids() { showq -u "$USER_NAME" 2>/dev/null | awk -v u="$USER_NAME" '$2==u {print $1}'; }
n_queued()   { queued_ids | wc -l; }
in_queue()   { queued_ids | grep -qx "$1"; }
ended() {    # jobid submitted_epoch
    [ "$1" = "-" ] && return 1
    in_queue "$1" && return 1
    [ $(( $(date +%s) - $2 )) -ge $GRACE ]
}

submit() {   # dir, msub args... ; prints the job id
    local dir=$1; shift
    (cd "$dir" && msub "$@" 2>>"$STATE/msub.err") | tr -d '[:space:]'
}

# ---- 1. mesh ----------------------------------------------------------------
read -r phase id tries t0 < <(get mesh)
if [ "$phase" = pending ]; then
    id=$(submit "$FAMILY" runMeshAll.sh)
    [ -n "$id" ] || { log "FAILED: msub runMeshAll.sh (see msub.err)"; exit 1; }
    put mesh running "$id" 1 "$(date +%s)"; t0=$(date +%s); phase=running
    log "submitted runMeshAll.sh as $id"
fi
while [ "$phase" = running ]; do
    if ended "$id" "$t0"; then
        if grep -q '^FAILED' "$FAMILY"/ugf-meshall.e"$id" 2>/dev/null || [ ! -f manifest-unigas.yaml ]; then
            put mesh failed "$id" 1 "$t0"; log "FAILED: runMeshAll.sh $id; see ugf-meshall.e$id. Stopping."
            rm -f "$STATE/sweep.pid"; exit 1
        fi
        put mesh done "$id" 1 "$t0"; phase=done; log "mesh $id done"
    else
        sleep "$POLL"
    fi
done
[ "$phase" = done ] || { log "mesh is '$phase'; stopping"; exit 1; }

# ---- 2-3. runs and posts ----------------------------------------------------
run_args() {  # case -> msub args for runCase.sh
    local short=${1%_np1x}
    if [ "$1" = Kn0p1_np1x ]; then
        echo "-q normalq -l walltime=24:00:00 -N ugf-$short runCase.sh"
    else
        echo "-N ugf-$short runCase.sh"
    fi
}
post_args() { echo "-q normalq -l nodes=1:ppn=32,walltime=24:00:00 -N ugf-post-${1%_np1x} runPost.sh"; }
finished() { python -m plumetools.cai2012ugf.hpc finished "$FAMILY/Cases/$1" >/dev/null 2>&1; }

while :; do
    remaining=0
    nq=$(n_queued)   # counted once per poll; showq can lag a fresh msub
    for c in "${CASES[@]}"; do
        read -r phase id tries t0 < <(get "$c")
        dir=$FAMILY/Cases/$c
        case "$phase" in
            pending|resume|needpost)
                remaining=1
                [ "$nq" -lt "$MAX_QUEUED" ] || continue
                nq=$((nq + 1))
                if [ "$phase" = needpost ]; then
                    # shellcheck disable=SC2046
                    new=$(submit "$dir" $(post_args "$c"))
                    [ -n "$new" ] && { put "$c" post "$new" "$tries" "$(date +%s)"; log "$c: submitted runPost.sh as $new"; } \
                                  || log "$c: msub runPost.sh failed (see msub.err); will retry"
                else
                    # shellcheck disable=SC2046
                    new=$(submit "$dir" $(run_args "$c"))
                    [ -n "$new" ] && { put "$c" running "$new" $((tries + 1)) "$(date +%s)"; log "$c: submitted runCase.sh as $new (try $((tries + 1)))"; } \
                                  || log "$c: msub runCase.sh failed (see msub.err); will retry"
                fi ;;
            running)
                remaining=1
                ended "$id" "$t0" || continue
                if finished "$c"; then
                    put "$c" needpost - "$tries" 0; log "$c: run $id finished"
                elif [ "$tries" -lt "$MAX_TRIES" ]; then
                    put "$c" resume - "$tries" 0; log "$c: run $id ended unfinished (walltime?); will resubmit"
                else
                    put "$c" failed "$id" "$tries" "$t0"; log "$c: FAILED after $tries run(s); see $dir/error.log and ugf-*.e$id"
                fi ;;
            post)
                remaining=1
                ended "$id" "$t0" || continue
                if grep -qs '^FAILED' "$dir"/*.e"$id"; then
                    put "$c" postfailed "$id" "$tries" "$t0"; log "$c: runPost $id FAILED; see $dir/error_post.log"
                else
                    put "$c" done "$id" "$tries" "$t0"; log "$c: done; results in $dir/results"
                fi ;;
        esac
    done
    [ "$remaining" -eq 1 ] || break
    sleep "$POLL"
done
log "sweep complete: $(cd "$STATE" && for c in "${CASES[@]}"; do printf '%s=%s ' "${c%_np1x}" "$(cut -d' ' -f1 "$c.state")"; done)"
rm -f "$STATE/sweep.pid"
