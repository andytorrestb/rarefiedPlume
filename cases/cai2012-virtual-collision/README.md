# cases/cai2012-virtual-collision — Cai, DSMC and a reduced collision model

A read-only assessment that compares three predictions of the Cai free plume:
Cai's collisionless evaluator, DSMC results of the `cai2012` and Kn-extension
families, and a **reduced virtual-collision model**
(`plumetools.virtualcollision`). It never generates cases, starts DSMC, submits
a job, or changes another family's policy gates. Comparisons are exploratory:
spatial convergence of the comparison has not been established.

```bash
python cases/cai2012-virtual-collision/assess.py inventory --output <scratch-dir>
python cases/cai2012-virtual-collision/assess.py smoke --output <scratch-dir> --no-plots
python cases/cai2012-virtual-collision/compare_exports.py --kn 10 1 --samples 2048
```

`results/smoke/` is the **committed** bounded smoke result; writing `smoke`
output there rewrites tracked files. `compare_exports.py`, `audit_hazard.py` and
`review_results.py` write into `cases/cai2012-kn-extension/…/results/virtual_collision/`.
`qsub runAssess.sh` runs `full` on a batch node (set `VC_HPC_ENV` for your site).

- Tutorial: [`docs/tutorials/cai2012-virtual-collision.md`](../../docs/tutorials/cai2012-virtual-collision.md)
- Mathematical contract and assessment protocol:
  [`docs/cai2012-virtual-collision.md`](../../docs/cai2012-virtual-collision.md)
- Working records: [`run-report.md`](run-report.md), [`refinement-plan.md`](refinement-plan.md),
  [`sampling-refinement-design.md`](sampling-refinement-design.md),
  [`sampling-implementation-prompt.md`](sampling-implementation-prompt.md)
