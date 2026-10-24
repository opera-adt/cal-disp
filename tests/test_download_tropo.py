"""Tests for TROPO staging: an empty scene search must not "complete"."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

pytest.importorskip("asf_search")

from cal_disp.download import _stage_tropo  # noqa: E402
from cal_disp.download._errors import DownloadError  # noqa: E402


def test_no_scenes_found_raises(tmp_path: Path):
    times = [datetime(2022, 1, 11, 0, 26), datetime(2022, 7, 22, 0, 26)]
    empty = _stage_tropo.asf.ASFSearchResults([])
    with patch.object(_stage_tropo, "find_nearest_scenes", return_value=empty):
        with patch.object(_stage_tropo.asf, "ASFSearchResults") as results_cls:
            with pytest.raises(DownloadError, match="No TROPO scenes found"):
                _stage_tropo.download_tropo(times, tmp_path)
    # Nothing was "downloaded" (the old code downloaded 0 scenes and returned)
    results_cls.return_value.download.assert_not_called()


def test_missing_times_listed_in_error(tmp_path: Path):
    found = MagicMock()
    found.__len__.return_value = 1
    found.__iter__.return_value = iter(["scene"])
    empty = _stage_tropo.asf.ASFSearchResults([])
    times = [datetime(2022, 1, 11), datetime(2022, 7, 22)]
    with patch.object(_stage_tropo, "find_nearest_scenes", side_effect=[found, empty]):
        with pytest.raises(DownloadError, match="1 of 2 sensing time") as exc_info:
            _stage_tropo.download_tropo(times, tmp_path)
    assert "2022-07-22" in str(exc_info.value)
    assert "2022-01-11" not in str(exc_info.value)


def test_scenes_found_are_downloaded(tmp_path: Path):
    found = MagicMock()
    found.__len__.return_value = 2
    found.__iter__.return_value = iter(["before", "after"])
    with patch.object(_stage_tropo, "find_nearest_scenes", return_value=found):
        with patch.object(_stage_tropo.asf, "ASFSearchResults") as results_cls:
            _stage_tropo.download_tropo(
                [datetime(2022, 1, 11)], tmp_path, num_workers=3
            )
    results_cls.assert_called_once_with(["before", "after"])
    results_cls.return_value.download.assert_called_once_with(
        path=tmp_path, processes=3
    )
