"""Plate-impingement Knudsen sweep on uniGasFoam (cases/plate-unigas).

Glue between two packages this one never edits:

* :mod:`plumetools.nozzleplume` -- the ``case.yaml`` schema, the plate and its
  symmetry planes, the snappyHexMesh octree, the nozzle patch, the run window;
* :mod:`plumetools.cai2012ugf` -- the uniGasFoam dictionaries, ``unigas.yaml``
  and the restart helpers.

Design: docs/superpowers/specs/2026-11-24-plate-impingement-kn-sweep-unigasfoam-design.md
"""
