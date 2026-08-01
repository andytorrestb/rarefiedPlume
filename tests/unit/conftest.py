"""Synthetic OpenFOAM meshes for parser tests.

Deliberately small and deliberately awkward: each fixture reproduces one shape of
input that broke the pre-refactor parsers.
"""

from __future__ import annotations

import pytest

from _helpers import (
    MIXED_FACES,
    QUAD_FACES,
    SQUARE_POINTS,
    TRI_FACES,
    make_case,
    patch,
)


@pytest.fixture
def tri_case(tmp_path):
    """Triangular faces only -- the shape every committed active mesh has."""
    return make_case(tmp_path / "tri", points=SQUARE_POINTS, faces=TRI_FACES,
                     patches=[patch("inflow", 2, 0)], sets={"inflow": [0, 1]})


@pytest.fixture
def quad_case(tmp_path):
    """Quadrilateral faces -- what blockMesh produces, and what AR-01 broke on."""
    return make_case(tmp_path / "quad", points=SQUARE_POINTS, faces=QUAD_FACES,
                     patches=[patch("inflow", 2, 0)], sets={"inflow": [0, 1]})


@pytest.fixture
def mixed_case(tmp_path):
    """Mixed 3- and 4-vertex faces, as in the archived 2d-wedge mesh."""
    return make_case(tmp_path / "mixed", points=SQUARE_POINTS, faces=MIXED_FACES,
                     patches=[patch("inflow", 3, 0)])


@pytest.fixture
def ingroups_case(tmp_path):
    """`inGroups` before `nFaces` -- blockMesh writes this; the legacy 3-line skip broke."""
    return make_case(tmp_path / "ingroups", points=SQUARE_POINTS, faces=TRI_FACES,
                     patches=[patch("inflow", 2, 0, extra="inGroups        1(inlet);\n")])


@pytest.fixture
def unspecified_case(tmp_path):
    """`type Unspecified;` -- the Pointwise export artifact of HA-07."""
    return make_case(tmp_path / "unspec", points=SQUARE_POINTS, faces=TRI_FACES,
                     patches=[patch("inflow", 2, 0, ptype="Unspecified")])


@pytest.fixture
def crlf_case(tmp_path):
    """CRLF line endings -- the legacy fixed-width slice mis-parsed nFaces (RB-02)."""
    return make_case(tmp_path / "crlf", points=SQUARE_POINTS, faces=TRI_FACES,
                     patches=[patch("inflow", 2, 0)], newline="\r\n",
                     sets={"inflow": [0, 1]})


@pytest.fixture
def prefix_case(tmp_path):
    """A patch whose name prefixes another -- the legacy substring search matched
    whichever came first (RB-03)."""
    return make_case(tmp_path / "prefix", points=SQUARE_POINTS, faces=MIXED_FACES,
                     patches=[patch("inflowOuter", 1, 0), patch("inflow", 2, 1)])


@pytest.fixture
def panel_case(tmp_path):
    """An extra patch, as iss_solar_panels_3 adds -- must not shift the others."""
    return make_case(tmp_path / "panel", points=SQUARE_POINTS, faces=MIXED_FACES,
                     patches=[patch("inflow", 1, 0),
                              patch("panel", 1, 1, ptype="Unspecified"),
                              patch("vacuum", 1, 2)])
