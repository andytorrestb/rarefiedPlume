"""The shared OpenFOAM readers: plumetools.foamio.times, and their error.

Twelve loops over a case's subdirectories decided, each in its own way, which
were time directories; foamio.times is the one rule now. FoamDataError is the
error the shared readers raise, and the names each family's handlers have
always caught still catch it.
"""

from __future__ import annotations

import pytest

from plumetools.foamio import times
from plumetools.foamio.errors import FoamDataError


# --------------------------------------------------------------------------- #
# time directories
# --------------------------------------------------------------------------- #

def test_time_directories_are_numeric_and_after_zero(tmp_path):
    for name in ("0", "0.0001", "2e-05", "0.05", "constant", "system", "0.orig",
                 "processor0"):
        (tmp_path / name).mkdir()
    (tmp_path / "0.07").write_text("a file, not a time")
    assert times.time_names(tmp_path) == ["2e-05", "0.0001", "0.05"]
    assert times.time_names(tmp_path, after=0.0001) == ["0.05"]
    assert times.latest_time(tmp_path) == "0.05"
    assert times.time_value("1e-3") == 1e-3 and times.time_value("constant") is None


def test_no_solver_output_is_an_error_that_says_so(tmp_path):
    (tmp_path / "0").mkdir()
    with pytest.raises(FoamDataError, match="solver wrote nothing"):
        times.latest_time(tmp_path)


def test_a_missing_directory_raises_unless_allowed(tmp_path):
    with pytest.raises(FileNotFoundError):
        times.time_directories(tmp_path / "absent")
    assert times.time_directories(tmp_path / "absent", missing_ok=True) == []


# --------------------------------------------------------------------------- #
# the old names still work, and are the shared objects
# --------------------------------------------------------------------------- #

def test_error_classes_keep_what_each_handler_caught():
    from plumetools.cai2012 import post as cai_post
    from plumetools.nozzleplume import post as nozzle_post
    from plumetools.rarefaction import foam

    assert nozzle_post.PostError is FoamDataError is foam.PostError
    assert issubclass(cai_post.PostError, FoamDataError)


def test_cai2012_readers_still_raise_their_own_error(tmp_path):
    from plumetools.cai2012 import post as cai_post

    with pytest.raises(cai_post.PostError, match="no field file"):
        cai_post.read_internal_field(tmp_path / "absent")
    (tmp_path / "0").mkdir()
    with pytest.raises(cai_post.PostError, match="wrote nothing"):
        cai_post.latest_time(tmp_path)
