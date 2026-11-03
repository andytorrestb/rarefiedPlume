# cases/cai2012-health — statistical health of the Cai Kn = 100 case

A 3 × 3 sweep over **parcels per exit cell** (5 / 20 / 40) and **averaging
duration** (0.5 / 1.5 / 4.5 domain transits) at fixed Kn = 100, judged through
imagery and the centreline error. `Cases/ppc020/s1p5` is `cases/cai2012`'s Kn100
case re-run and is the control. Everything else is fixed in
[`baseCase/case.yaml`](baseCase/case.yaml); `study.yaml` holds the two axes and
nothing else.

```bash
./generate_cases.py --dry-run     # the matrix and its predicted cost; no CaseFoam needed
./generate_cases.py               # nine cases via CaseFoam (pip install -e ".[cases]")
./AllmeshCases
./AllrunCases --only Cases/ppc005/s0p5 --only Cases/ppc005/s1p5   # cheapest pair first
./AllrunCases                     # all nine, cheapest first, timed into results/cost-model.yaml
./AllpostCases                    # validate, audit, render, dissect (ffmpeg, pvpython)
```

`AllrunCases` prunes reconstructed parcel clouds from all but the newest time
(`--keep-clouds` to keep them). `results/cost-model.yaml` is tracked: it is the
measured wall clock of every case run so far.

- Tutorial: [`docs/tutorials/cai2012-health.md`](../../docs/tutorials/cai2012-health.md)
- Design, costs, imagery and recorded results (moved from this README):
  [`docs/cai2012-health.md`](../../docs/cai2012-health.md)
