"""Tests for base configuration models (InputFileGroup, ancillary groups)."""

from __future__ import annotations

from pathlib import Path

import pytest

from cal_disp.config._base import (
    DynamicAncillaryFileGroup,
    InputFileGroup,
    StaticAncillaryFileGroup,
)


class TestInputFileGroup:
    """Tests for InputFileGroup configuration."""

    def test_creates_with_required_fields(
        self,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
    ):
        """Should create with all required fields."""
        config = InputFileGroup(
            disp_file=sample_disp_product,
            unr_grid_latlon_file=sample_unr_grid_latlon,
            unr_timeseries_dir=sample_unr_timeseries_dir,
            frame_id=8882,
            unr_grid_version="0.2",
            unr_grid_type="constant",
        )

        assert config.disp_file == sample_disp_product
        assert config.unr_grid_latlon_file == sample_unr_grid_latlon
        assert config.frame_id == 8882

    def test_validates_disp_extension(self, tmp_path: Path):
        """Should reject invalid DISP file extensions."""
        bad_file = tmp_path / "disp.txt"
        bad_file.touch()

        with pytest.raises(ValueError, match="must be .nc or .h5"):
            InputFileGroup(
                disp_file=bad_file,
                unr_grid_latlon_file=tmp_path / "grid_latlon_lookup_v0.2.txt",
                unr_timeseries_dir=tmp_path,
                frame_id=8882,
                unr_grid_version="0.2",
                unr_grid_type="constant",
            )

    def test_accepts_nc_extension(self, sample_unr_grid_latlon: Path, tmp_path: Path):
        """Should accept .nc extension."""
        nc_file = tmp_path / "disp.nc"
        nc_file.touch()
        (tmp_path / "test.tenv8").touch()

        config = InputFileGroup(
            disp_file=nc_file,
            unr_grid_latlon_file=sample_unr_grid_latlon,
            unr_timeseries_dir=tmp_path,
            frame_id=8882,
            unr_grid_version="0.2",
            unr_grid_type="constant",
        )

        assert config.disp_file == nc_file

    def test_accepts_h5_extension(self, sample_unr_grid_latlon: Path, tmp_path: Path):
        """Should accept .h5 extension."""
        h5_file = tmp_path / "disp.h5"
        h5_file.touch()
        (tmp_path / "test.tenv8").touch()

        config = InputFileGroup(
            disp_file=h5_file,
            unr_grid_latlon_file=sample_unr_grid_latlon,
            unr_timeseries_dir=tmp_path,
            frame_id=8882,
            unr_grid_version="0.2",
            unr_grid_type="constant",
        )

        assert config.disp_file == h5_file

    def test_validates_latlon_filename(
        self,
        sample_disp_product: Path,
        sample_unr_timeseries_dir: Path,
        tmp_path: Path,
    ):
        """Should validate latlon file starts with grid_latlon_lookup."""
        bad_file = tmp_path / "wrong_name.txt"
        bad_file.touch()

        with pytest.raises(ValueError, match="Expected grid_latlon_lookup"):
            InputFileGroup(
                disp_file=sample_disp_product,
                unr_grid_latlon_file=bad_file,
                unr_timeseries_dir=sample_unr_timeseries_dir,
                frame_id=8882,
                unr_grid_version="0.2",
                unr_grid_type="constant",
            )

    def test_validates_latlon_extension(
        self,
        sample_disp_product: Path,
        sample_unr_timeseries_dir: Path,
        tmp_path: Path,
    ):
        """Should validate latlon file has .txt extension."""
        bad_file = tmp_path / "grid_latlon_lookup_v0.2.csv"
        bad_file.touch()

        with pytest.raises(ValueError, match="must be .txt"):
            InputFileGroup(
                disp_file=sample_disp_product,
                unr_grid_latlon_file=bad_file,
                unr_timeseries_dir=sample_unr_timeseries_dir,
                frame_id=8882,
                unr_grid_version="0.2",
                unr_grid_type="constant",
            )

    def test_validates_grid_dir_exists(
        self,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        tmp_path: Path,
    ):
        """Should validate grid directory exists."""
        nonexistent = tmp_path / "nonexistent"

        with pytest.raises(ValueError, match="does not exist"):
            InputFileGroup(
                disp_file=sample_disp_product,
                unr_grid_latlon_file=sample_unr_grid_latlon,
                unr_timeseries_dir=nonexistent,
                frame_id=8882,
                unr_grid_version="0.2",
                unr_grid_type="constant",
            )

    def test_validates_grid_dir_has_tenv8_files(
        self,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        tmp_path: Path,
    ):
        """Should validate grid directory contains .tenv8 files."""
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()

        with pytest.raises(ValueError, match="No .tenv8 files found"):
            InputFileGroup(
                disp_file=sample_disp_product,
                unr_grid_latlon_file=sample_unr_grid_latlon,
                unr_timeseries_dir=empty_dir,
                frame_id=8882,
                unr_grid_version="0.2",
                unr_grid_type="constant",
            )

    def test_validates_frame_id_range(
        self,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
    ):
        """Should validate frame ID is within reasonable range."""
        with pytest.raises(ValueError, match="between 1 and 99999"):
            InputFileGroup(
                disp_file=sample_disp_product,
                unr_grid_latlon_file=sample_unr_grid_latlon,
                unr_timeseries_dir=sample_unr_timeseries_dir,
                frame_id=0,
                unr_grid_version="0.2",
                unr_grid_type="constant",
            )

        with pytest.raises(ValueError, match="between 1 and 99999"):
            InputFileGroup(
                disp_file=sample_disp_product,
                unr_grid_latlon_file=sample_unr_grid_latlon,
                unr_timeseries_dir=sample_unr_timeseries_dir,
                frame_id=100000,
                unr_grid_version="0.2",
                unr_grid_type="constant",
            )

    def test_accepts_valid_frame_ids(
        self,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
    ):
        """Should accept valid frame IDs."""
        for frame_id in [1, 100, 8882, 99999]:
            config = InputFileGroup(
                disp_file=sample_disp_product,
                unr_grid_latlon_file=sample_unr_grid_latlon,
                unr_timeseries_dir=sample_unr_timeseries_dir,
                frame_id=frame_id,
                unr_grid_version="0.2",
                unr_grid_type="constant",
            )
            assert config.frame_id == frame_id

    def test_forbids_extra_fields(
        self,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
    ):
        """Should forbid extra fields not in schema."""
        with pytest.raises(ValueError):
            InputFileGroup(
                disp_file=sample_disp_product,
                unr_grid_latlon_file=sample_unr_grid_latlon,
                unr_timeseries_dir=sample_unr_timeseries_dir,
                frame_id=8882,
                unr_grid_version="0.2",
                unr_grid_type="constant",
                extra_field="not allowed",
            )


class TestDynamicAncillaryFileGroup:
    """Tests for DynamicAncillaryFileGroup configuration."""

    def test_creates_with_required_fields(
        self,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
    ):
        """Should create with only required fields."""
        los_file, dem_file = sample_static_layers

        config = DynamicAncillaryFileGroup(
            algorithm_parameters_file=sample_algorithm_params,
            static_los_file=los_file,
            static_dem_file=dem_file,
        )

        assert config.algorithm_parameters_file == sample_algorithm_params
        assert config.los_file == los_file
        assert config.dem_file == dem_file

    def test_los_file_alias(
        self,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
    ):
        """Should accept static_los_file as alias."""
        los_file, dem_file = sample_static_layers

        config = DynamicAncillaryFileGroup(
            algorithm_parameters_file=sample_algorithm_params,
            static_los_file=los_file,  # Using alias
            static_dem_file=dem_file,
        )

        assert config.los_file == los_file

    def test_dem_file_alias(
        self,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
    ):
        """Should accept static_dem_file as alias."""
        los_file, dem_file = sample_static_layers

        config = DynamicAncillaryFileGroup(
            algorithm_parameters_file=sample_algorithm_params,
            static_los_file=los_file,
            static_dem_file=dem_file,  # Using alias
        )

        assert config.dem_file == dem_file

    def test_optional_mask_file(
        self,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
        sample_mask_file: Path,
    ):
        """Should accept optional mask file."""
        los_file, dem_file = sample_static_layers

        config = DynamicAncillaryFileGroup(
            algorithm_parameters_file=sample_algorithm_params,
            static_los_file=los_file,
            static_dem_file=dem_file,
            mask_file=sample_mask_file,
        )

        assert config.mask_file == sample_mask_file

    def test_optional_tropo_files(
        self,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
        tmp_path: Path,
    ):
        """Should accept optional tropo file lists."""
        los_file, dem_file = sample_static_layers

        # Create mock tropo files
        ref_tropo1 = tmp_path / "ref_tropo_1.h5"
        ref_tropo2 = tmp_path / "ref_tropo_2.h5"
        sec_tropo = tmp_path / "sec_tropo.h5"
        ref_tropo1.touch()
        ref_tropo2.touch()
        sec_tropo.touch()

        config = DynamicAncillaryFileGroup(
            algorithm_parameters_file=sample_algorithm_params,
            static_los_file=los_file,
            static_dem_file=dem_file,
            reference_tropo_files=[ref_tropo1, ref_tropo2],
            secondary_tropo_files=[sec_tropo],
        )

        assert len(config.reference_tropo_files) == 2
        assert len(config.secondary_tropo_files) == 1

    def test_tropo_alias(
        self,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
        tmp_path: Path,
    ):
        """Should accept ref_tropo_files and sec_tropo_files aliases."""
        los_file, dem_file = sample_static_layers

        tropo_file = tmp_path / "tropo.h5"
        tropo_file.touch()

        config = DynamicAncillaryFileGroup(
            algorithm_parameters_file=sample_algorithm_params,
            static_los_file=los_file,
            static_dem_file=dem_file,
            ref_tropo_files=[tropo_file],  # Using alias
            sec_tropo_files=[tropo_file],  # Using alias
        )

        assert config.reference_tropo_files == [tropo_file]
        assert config.secondary_tropo_files == [tropo_file]

    def test_optional_iono_files(
        self,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
        tmp_path: Path,
    ):
        """Should accept optional ionospheric correction files."""
        los_file, dem_file = sample_static_layers

        iono_file = tmp_path / "iono.h5"
        iono_file.touch()

        config = DynamicAncillaryFileGroup(
            algorithm_parameters_file=sample_algorithm_params,
            static_los_file=los_file,
            static_dem_file=dem_file,
            iono_files=[iono_file],
        )

        assert config.iono_files == [iono_file]

    def test_optional_tiles_files(
        self,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
        tmp_path: Path,
    ):
        """Should accept optional tile boundary files."""
        los_file, dem_file = sample_static_layers

        tile_file = tmp_path / "burst_bounds.geojson"
        tile_file.touch()

        config = DynamicAncillaryFileGroup(
            algorithm_parameters_file=sample_algorithm_params,
            static_los_file=los_file,
            static_dem_file=dem_file,
            tiles_files=[tile_file],
        )

        assert config.tiles_files == [tile_file]

    def test_get_all_files(
        self,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
        tmp_path: Path,
    ):
        """Should return all file paths including lists."""
        los_file, dem_file = sample_static_layers

        tropo_file = tmp_path / "tropo.h5"
        tropo_file.touch()

        config = DynamicAncillaryFileGroup(
            algorithm_parameters_file=sample_algorithm_params,
            static_los_file=los_file,
            static_dem_file=dem_file,
            reference_tropo_files=[tropo_file],
        )

        files = config.get_all_files()

        assert "algorithm_parameters_file" in files
        assert "los_file" in files
        assert "dem_file" in files
        assert "reference_tropo_files[0]" in files  # Flattened

    def test_none_values_excluded_from_get_all_files(
        self,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
    ):
        """Should not include None optional fields in get_all_files."""
        los_file, dem_file = sample_static_layers

        config = DynamicAncillaryFileGroup(
            algorithm_parameters_file=sample_algorithm_params,
            static_los_file=los_file,
            static_dem_file=dem_file,
            # All optional fields are None
        )

        files = config.get_all_files()

        assert "mask_file" not in files
        assert "reference_tropo_files" not in files


class TestStaticAncillaryFileGroup:
    """Tests for StaticAncillaryFileGroup configuration."""

    def test_all_fields_optional(self):
        """Should create with all fields as None."""
        config = StaticAncillaryFileGroup()

        assert config.algorithm_parameters_overrides_json is None
        assert config.deformation_area_database_json is None
        assert config.event_database_json is None

    def test_algorithm_overrides(self, tmp_path: Path):
        """Should accept algorithm overrides JSON file."""
        overrides_file = tmp_path / "overrides.json"
        overrides_file.write_text('{"frame_8882": {"downsample_factor": 5}}')

        config = StaticAncillaryFileGroup(
            algorithm_parameters_overrides_json=overrides_file
        )

        assert config.algorithm_parameters_overrides_json == overrides_file

    def test_deformation_database(self, tmp_path: Path):
        """Should accept deformation area database."""
        defo_file = tmp_path / "defo_areas.geojson"
        defo_file.write_text('{"type": "FeatureCollection", "features": []}')

        config = StaticAncillaryFileGroup(deformation_area_database_json=defo_file)

        assert config.deformation_area_database_json == defo_file

    def test_defo_area_alias(self, tmp_path: Path):
        """Should accept defo_area_db_json as alias."""
        defo_file = tmp_path / "defo_areas.geojson"
        defo_file.write_text('{"type": "FeatureCollection", "features": []}')

        config = StaticAncillaryFileGroup(defo_area_db_json=defo_file)

        assert config.deformation_area_database_json == defo_file

    def test_event_database(self, tmp_path: Path):
        """Should accept event database."""
        event_file = tmp_path / "events.geojson"
        event_file.write_text('{"type": "FeatureCollection", "features": []}')

        config = StaticAncillaryFileGroup(event_database_json=event_file)

        assert config.event_database_json == event_file

    def test_event_db_alias(self, tmp_path: Path):
        """Should accept event_db_json as alias."""
        event_file = tmp_path / "events.geojson"
        event_file.write_text('{"type": "FeatureCollection", "features": []}')

        config = StaticAncillaryFileGroup(event_db_json=event_file)

        assert config.event_database_json == event_file

    def test_has_algorithm_overrides(self, tmp_path: Path):
        """Should correctly report if algorithm overrides are present."""
        config_without = StaticAncillaryFileGroup()
        assert not config_without.has_algorithm_overrides()

        overrides_file = tmp_path / "overrides.json"
        overrides_file.write_text("{}")

        config_with = StaticAncillaryFileGroup(
            algorithm_parameters_overrides_json=overrides_file
        )
        assert config_with.has_algorithm_overrides()

    def test_has_deformation_database(self, tmp_path: Path):
        """Should correctly report if deformation database is present."""
        config_without = StaticAncillaryFileGroup()
        assert not config_without.has_deformation_database()

        defo_file = tmp_path / "defo.geojson"
        defo_file.write_text("{}")

        config_with = StaticAncillaryFileGroup(deformation_area_database_json=defo_file)
        assert config_with.has_deformation_database()

    def test_has_event_database(self, tmp_path: Path):
        """Should correctly report if event database is present."""
        config_without = StaticAncillaryFileGroup()
        assert not config_without.has_event_database()

        event_file = tmp_path / "events.geojson"
        event_file.write_text("{}")

        config_with = StaticAncillaryFileGroup(event_database_json=event_file)
        assert config_with.has_event_database()

    def test_get_all_files_empty(self):
        """Should return empty dict when no files configured."""
        config = StaticAncillaryFileGroup()

        files = config.get_all_files()

        assert files == {}

    def test_get_all_files_with_values(self, tmp_path: Path):
        """Should return all configured files."""
        overrides = tmp_path / "overrides.json"
        defo = tmp_path / "defo.geojson"
        events = tmp_path / "events.geojson"

        for f in [overrides, defo, events]:
            f.write_text("{}")

        config = StaticAncillaryFileGroup(
            algorithm_parameters_overrides_json=overrides,
            deformation_area_database_json=defo,
            event_database_json=events,
        )

        files = config.get_all_files()

        assert len(files) == 3
        assert "algorithm_parameters_overrides_json" in files
        assert "deformation_area_database_json" in files
        assert "event_database_json" in files
