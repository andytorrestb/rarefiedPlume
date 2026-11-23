"""Case figures: a check that could not judge anything draws no figure.

A run that keeps only its last couple of sampled saves (the laptop runs, which
have no disk for more) leaves the stationarity check ``insufficient-frames``
with no regions judged. Post-processing must carry on past that, not crash.
"""

from __future__ import annotations

import pytest

from plumetools.cai2012ext import convergence, figures

pytest.importorskip("matplotlib")


def test_stationarity_figure_skipped_when_too_few_frames(tmp_path):
    frames = [{"time": t, "mean": {"core": 1.0}, "instant": {"core": 1.0}}
              for t in (0.062689754, 0.06304822)]
    assessment = convergence.assess(frames)
    assert assessment["verdict"] == "insufficient-frames"

    written = figures.plot_stationarity(tmp_path, frames, assessment,
                                        transit_s=2.8e-3, start_s=0.0501,
                                        label="Kn1_np1x")

    assert written == []
    assert not (tmp_path / "stationarity.png").exists()
