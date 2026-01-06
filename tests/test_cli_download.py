"""Tests for cal-disp download CLI commands."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from cal_disp.cli.download import burst_bounds, disp_s1, download_group, tropo, unr


class TestDispS1Download:
    """Tests for disp-s1 download command."""

    def test_requires_frame_id(self, cli_runner: CliRunner):
        """Should require --frame-id option."""
        result = cli_runner.invoke(disp_s1, [])

        assert result.exit_code != 0
        assert "frame" in result.output.lower() or "required" in result.output.lower()

    def test_requires_output_dir(self, cli_runner: CliRunner):
        """Should require --output-dir option."""
        result = cli_runner.invoke(disp_s1, ["--frame-id", "8882"])

        assert result.exit_code != 0

    def test_basic_invocation(self, cli_runner: CliRunner, tmp_path: Path):
        """Should run with minimal required arguments."""
        with patch("cal_disp.download.download_disp") as mock_download:
            result = cli_runner.invoke(
                disp_s1,
                [
                    "--frame-id",
                    "8882",
                    "--output-dir",
                    str(tmp_path),
                    "-s",
                    "2022-07-22",
                    "-e",
                    "2022-07-22",
                ],
            )

            assert result.exit_code == 0
            assert "Download complete" in result.output
            mock_download.assert_called_once()

    def test_with_date_range(self, cli_runner: CliRunner, tmp_path: Path):
        """Should accept start and end dates."""
        with patch("cal_disp.download.download_disp") as mock_download:
            result = cli_runner.invoke(
                disp_s1,
                [
                    "--frame-id",
                    "8882",
                    "--output-dir",
                    str(tmp_path),
                    "--start",
                    "2022-07-01",
                    "--end",
                    "2022-07-31",
                ],
            )

            assert result.exit_code == 0
            call_kwargs = mock_download.call_args.kwargs
            assert isinstance(call_kwargs["start"], datetime)
            assert isinstance(call_kwargs["end"], datetime)

    def test_with_multiple_workers(self, cli_runner: CliRunner, tmp_path: Path):
        """Should accept --num-workers option."""
        with patch("cal_disp.download.download_disp") as mock_download:
            result = cli_runner.invoke(
                disp_s1,
                [
                    "--frame-id",
                    "8882",
                    "--output-dir",
                    str(tmp_path),
                    "--num-workers",
                    "8",
                    "--start",
                    "2022-07-01",
                    "--end",
                    "2022-08-31",
                ],
            )

            assert result.exit_code == 0
            call_kwargs = mock_download.call_args.kwargs
            assert call_kwargs["num_workers"] == 8

    def test_invalid_date_format(self, cli_runner: CliRunner, tmp_path: Path):
        """Should reject invalid date format."""
        result = cli_runner.invoke(
            disp_s1,
            [
                "--frame-id",
                "8882",
                "--output-dir",
                str(tmp_path),
                "--start",
                "invalid-date",
            ],
        )

        assert result.exit_code != 0


class TestUNRDownload:
    """Tests for unr download command."""

    def test_requires_frame_id(self, cli_runner: CliRunner):
        """Should require --frame-id option."""
        result = cli_runner.invoke(unr, [])

        assert result.exit_code != 0

    def test_requires_output_dir(self, cli_runner: CliRunner):
        """Should require --output-dir option."""
        result = cli_runner.invoke(unr, ["--frame-id", "8882"])

        assert result.exit_code != 0

    def test_basic_invocation(self, cli_runner: CliRunner, tmp_path: Path):
        """Should run with minimal arguments."""
        with patch("cal_disp.download.download_unr_grid") as mock_download:
            result = cli_runner.invoke(
                unr,
                [
                    "--frame-id",
                    "8882",
                    "--output-dir",
                    str(tmp_path),
                ],
            )

            assert result.exit_code == 0
            assert "Download complete" in result.output
            mock_download.assert_called_once()

    def test_with_custom_margin(self, cli_runner: CliRunner, tmp_path: Path):
        """Should accept --margin option."""
        with patch("cal_disp.download.download_unr_grid") as mock_download:
            result = cli_runner.invoke(
                unr,
                [
                    "--frame-id",
                    "8882",
                    "--output-dir",
                    str(tmp_path),
                    "--margin",
                    "1.0",
                ],
            )

            assert result.exit_code == 0
            call_kwargs = mock_download.call_args.kwargs
            assert call_kwargs["margin_deg"] == 1.0

    def test_creates_output_directory(self, cli_runner: CliRunner, tmp_path: Path):
        """Should create output directory if it doesn't exist."""
        output_dir = tmp_path / "nonexistent" / "dir"
        assert not output_dir.exists()

        with patch("cal_disp.download.download_unr_grid"):
            result = cli_runner.invoke(
                unr,
                [
                    "--frame-id",
                    "8882",
                    "--output-dir",
                    str(output_dir),
                ],
            )

            assert result.exit_code == 0
            assert output_dir.exists()


class TestTropoDownload:
    """Tests for tropo download command."""

    def test_requires_input_file(self, cli_runner: CliRunner):
        """Should require --input-file option."""
        result = cli_runner.invoke(tropo, [])

        assert result.exit_code != 0

    def test_requires_output_dir(
        self, cli_runner: CliRunner, sample_disp_product: Path
    ):
        """Should require --output-dir option."""
        result = cli_runner.invoke(
            tropo,
            ["--input-file", str(sample_disp_product)],
        )

        assert result.exit_code != 0

    def test_basic_invocation(
        self,
        cli_runner: CliRunner,
        tmp_path: Path,
        sample_disp_product: Path,
    ):
        """Should run with minimal arguments."""
        with patch("cal_disp.download.download_tropo") as mock_download:
            with patch(
                "cal_disp.download.utils.extract_sensing_times_from_file"
            ) as mock_extract:
                mock_extract.return_value = [
                    datetime(2022, 1, 11),
                    datetime(2022, 7, 22),
                ]

                result = cli_runner.invoke(
                    tropo,
                    [
                        "--input-file",
                        str(sample_disp_product),
                        "--output-dir",
                        str(tmp_path),
                    ],
                )

                assert result.exit_code == 0
                assert "Download complete" in result.output
                mock_download.assert_called_once()

    def test_with_interp_flag(
        self,
        cli_runner: CliRunner,
        tmp_path: Path,
        sample_disp_product: Path,
    ):
        """Should accept --interp flag."""
        with patch("cal_disp.download.download_tropo") as mock_download:
            with patch(
                "cal_disp.download.utils.extract_sensing_times_from_file"
            ) as mock_extract:
                mock_extract.return_value = [datetime(2022, 1, 11)]

                result = cli_runner.invoke(
                    tropo,
                    [
                        "--input-file",
                        str(sample_disp_product),
                        "--output-dir",
                        str(tmp_path),
                        "--interp",
                    ],
                )

                assert result.exit_code == 0
                call_kwargs = mock_download.call_args.kwargs
                assert call_kwargs["interp"] is True

    def test_with_custom_workers(
        self,
        cli_runner: CliRunner,
        tmp_path: Path,
        sample_disp_product: Path,
    ):
        """Should accept --num-workers option."""
        with patch("cal_disp.download.download_tropo") as mock_download:
            with patch(
                "cal_disp.download.utils.extract_sensing_times_from_file"
            ) as mock_extract:
                mock_extract.return_value = [datetime(2022, 1, 11)]

                result = cli_runner.invoke(
                    tropo,
                    [
                        "--input-file",
                        str(sample_disp_product),
                        "--output-dir",
                        str(tmp_path),
                        "--num-workers",
                        "8",
                    ],
                )

                assert result.exit_code == 0
                call_kwargs = mock_download.call_args.kwargs
                assert call_kwargs["num_workers"] == 8


class TestBurstBoundsDownload:
    """Tests for burst-bounds download command."""

    def test_requires_input_file(self, cli_runner: CliRunner):
        """Should require --input-file option."""
        result = cli_runner.invoke(burst_bounds, [])

        assert result.exit_code != 0

    def test_requires_output_dir(
        self, cli_runner: CliRunner, sample_disp_product: Path
    ):
        """Should require --output-dir option."""
        result = cli_runner.invoke(
            burst_bounds,
            ["--input-file", str(sample_disp_product)],
        )

        assert result.exit_code != 0

    def test_basic_invocation(
        self,
        cli_runner: CliRunner,
        tmp_path: Path,
        sample_disp_product: Path,
    ):
        """Should run with minimal arguments."""
        with patch("cal_disp.download.generate_s1_burst_tiles") as mock_download:
            with patch(
                "cal_disp.download.utils.extract_sensing_times_from_file"
            ) as mock_extract:
                mock_extract.return_value = [
                    datetime(2022, 1, 11),
                    datetime(2022, 7, 22),
                ]

                result = cli_runner.invoke(
                    burst_bounds,
                    [
                        "--input-file",
                        str(sample_disp_product),
                        "--output-dir",
                        str(tmp_path),
                    ],
                )

                assert result.exit_code == 0
                assert "Download complete" in result.output
                # Should be called twice (once per sensing time)
                assert mock_download.call_count == 2

    def test_rejects_non_s1_product(
        self,
        cli_runner: CliRunner,
        tmp_path: Path,
    ):
        """Should reject non-S1 products."""
        # Create a non-S1 filename
        non_s1_file = (
            tmp_path
            / "OPERA_L3_DISP-NI_IW_F08882_VV_20220111T002651Z_20220722T002657Z_v1.0.nc"
        )
        non_s1_file.touch()

        result = cli_runner.invoke(
            burst_bounds,
            [
                "--input-file",
                str(non_s1_file),
                "--output-dir",
                str(tmp_path),
            ],
        )

        assert result.exit_code != 0
        assert "Only DISP-S1 products supported" in result.output

    def test_extracts_frame_id_from_filename(
        self,
        cli_runner: CliRunner,
        tmp_path: Path,
        sample_disp_product: Path,
    ):
        """Should correctly extract frame ID from filename."""
        with patch("cal_disp.download.generate_s1_burst_tiles") as mock_download:
            with patch(
                "cal_disp.download.utils.extract_sensing_times_from_file"
            ) as mock_extract:
                mock_extract.return_value = [datetime(2022, 1, 11)]

                cli_runner.invoke(
                    burst_bounds,
                    [
                        "--input-file",
                        str(sample_disp_product),
                        "--output-dir",
                        str(tmp_path),
                    ],
                )

                call_kwargs = mock_download.call_args.kwargs
                assert call_kwargs["frame_id"] == 8882


class TestDownloadGroup:
    """Tests for download command group."""

    def test_shows_help(self, cli_runner: CliRunner):
        """Should display help for download group."""
        result = cli_runner.invoke(download_group, ["--help"])

        assert result.exit_code == 0
        assert "disp-s1" in result.output
        assert "unr" in result.output
        assert "tropo" in result.output
        assert "burst-bounds" in result.output

    def test_all_subcommands_have_help(self, cli_runner: CliRunner):
        """All subcommands should have help text."""
        for cmd in [disp_s1, unr, tropo, burst_bounds]:
            result = cli_runner.invoke(cmd, ["--help"])
            assert result.exit_code == 0
            assert len(result.output) > 0
