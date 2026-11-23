"""The Cai plume Knudsen sweep on uniGasFoam's hybrid USP-DSMC solver.

The ``cases/cai2012-unigas`` family runs the same seven cases as
``cases/cai2012-kn-extension`` -- the same physics, meshes, particle weight and
run window -- on a patched uniGasFoam instead of standard ``dsmcFoam``.

This package adds only what uniGasFoam needs: its dictionaries, a lossless
conversion of its output fields to the names the existing analysis reads, a
post-run suitability verdict, restart and leg control, and scaling and report
adapters. It imports :mod:`plumetools.cai2012` and :mod:`plumetools.cai2012ext`
and never edits them: those packages are the frozen record of what was run.

Design: docs/superpowers/specs/2026-11-19-cai-kn-sweep-unigasfoam-design.md.
Nothing produced with it is physically validated.
"""
