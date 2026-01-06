"""Tests for WorkerSettings configuration."""

from __future__ import annotations

import multiprocessing

import pytest

from cal_disp.config._workers import WorkerSettings


class TestWorkerSettingsBasics:
    """Tests for basic WorkerSettings functionality."""

    def test_defaults(self):
        """Should have expected default values."""
        settings = WorkerSettings()

        assert settings.n_workers == 4
        assert settings.threads_per_worker == 2
        assert settings.block_shape == (128, 128)

    def test_custom_values(self):
        """Should accept custom values."""
        settings = WorkerSettings(
            n_workers=8,
            threads_per_worker=4,
            block_shape=(256, 256),
        )

        assert settings.n_workers == 8
        assert settings.threads_per_worker == 4
        assert settings.block_shape == (256, 256)

    def test_total_threads_property(self):
        """Should correctly calculate total threads."""
        settings = WorkerSettings(n_workers=4, threads_per_worker=3)

        assert settings.total_threads == 12

    def test_block_size_property(self):
        """Should correctly calculate block size."""
        settings = WorkerSettings(block_shape=(100, 200))

        assert settings.block_size == 20000


class TestWorkerSettingsValidation:
    """Tests for WorkerSettings validation."""

    def test_rejects_zero_workers(self):
        """Should reject zero workers."""
        with pytest.raises(ValueError, match="greater than or equal to 1"):
            WorkerSettings(n_workers=0)

    def test_rejects_negative_workers(self):
        """Should reject negative workers."""
        with pytest.raises(ValueError, match="greater than or equal to 1"):
            WorkerSettings(n_workers=-1)

    def test_rejects_too_many_workers(self):
        """Should reject unreasonably high worker count."""
        with pytest.raises(ValueError, match="less than or equal to 256"):
            WorkerSettings(n_workers=300)

    def test_rejects_zero_threads(self):
        """Should reject zero threads per worker."""
        with pytest.raises(ValueError, match="greater than or equal to 1"):
            WorkerSettings(threads_per_worker=0)

    def test_rejects_too_many_threads(self):
        """Should reject unreasonably high threads per worker."""
        with pytest.raises(ValueError, match="less than or equal to 32"):
            WorkerSettings(threads_per_worker=50)

    def test_validates_block_shape_length(self):
        """Should require exactly 2 dimensions."""
        with pytest.raises(ValueError, match="exactly 2 dimensions"):
            WorkerSettings(block_shape=(128,))  # type: ignore

        with pytest.raises(ValueError, match="exactly 2 dimensions"):
            WorkerSettings(block_shape=(128, 128, 128))  # type: ignore

    def test_validates_block_shape_positive(self):
        """Should require positive block dimensions."""
        with pytest.raises(ValueError, match="must be positive"):
            WorkerSettings(block_shape=(0, 128))

        with pytest.raises(ValueError, match="must be positive"):
            WorkerSettings(block_shape=(128, -10))

    def test_validates_block_shape_integers(self):
        """Should require integer block dimensions."""
        with pytest.raises(ValueError, match="must be integers"):
            WorkerSettings(block_shape=(128.5, 128))  # type: ignore

    def test_accepts_list_for_block_shape(self):
        """Should convert list to tuple for block_shape."""
        settings = WorkerSettings(block_shape=[100, 200])  # type: ignore

        assert settings.block_shape == (100, 200)
        assert isinstance(settings.block_shape, tuple)


class TestWorkerSettingsMemoryEstimation:
    """Tests for memory estimation methods."""

    def test_estimate_memory_per_block_default(self):
        """Should estimate memory for default dtype."""
        settings = WorkerSettings(block_shape=(256, 256))

        # 256*256 * 8 bytes = 524288 bytes = 0.5 MB
        mem_mb = settings.estimate_memory_per_block()

        assert mem_mb == pytest.approx(0.5, rel=0.01)

    def test_estimate_memory_per_block_float32(self):
        """Should estimate memory for float32."""
        settings = WorkerSettings(block_shape=(256, 256))

        # 256*256 * 4 bytes = 0.25 MB
        mem_mb = settings.estimate_memory_per_block(dtype_size=4)

        assert mem_mb == pytest.approx(0.25, rel=0.01)

    def test_estimate_memory_per_block_multiple_bands(self):
        """Should estimate memory for multiple bands."""
        settings = WorkerSettings(block_shape=(256, 256))

        # 256*256 * 8 bytes * 3 bands = 1.5 MB
        mem_mb = settings.estimate_memory_per_block(dtype_size=8, n_bands=3)

        assert mem_mb == pytest.approx(1.5, rel=0.01)

    def test_estimate_total_memory(self):
        """Should estimate total memory across workers."""
        settings = WorkerSettings(n_workers=4, block_shape=(256, 256))

        # 0.5 MB/block * 4 workers * 1.5 overhead = 3 MB = 0.003 GB
        total_gb = settings.estimate_total_memory()

        assert total_gb == pytest.approx(0.003, rel=0.05)

    def test_estimate_total_memory_custom_overhead(self):
        """Should accept custom overhead factor."""
        settings = WorkerSettings(n_workers=4, block_shape=(256, 256))

        total_gb = settings.estimate_total_memory(overhead_factor=2.0)

        # Should be higher with 2.0 vs 1.5 overhead
        assert total_gb > settings.estimate_total_memory(overhead_factor=1.5)


class TestWorkerSettingsClassMethods:
    """Tests for WorkerSettings factory methods."""

    def test_create_lightweight(self):
        """Should create lightweight configuration."""
        settings = WorkerSettings.create_lightweight()

        assert settings.n_workers == 2
        assert settings.threads_per_worker == 1
        assert settings.block_shape == (64, 64)

    def test_create_standard(self):
        """Should create standard configuration."""
        settings = WorkerSettings.create_standard()

        assert settings.n_workers == 4
        assert settings.threads_per_worker == 2
        assert settings.block_shape == (128, 128)

    def test_create_heavy(self):
        """Should create heavy-duty configuration."""
        settings = WorkerSettings.create_heavy()

        assert settings.n_workers == 8
        assert settings.threads_per_worker == 4
        assert settings.block_shape == (256, 256)

    def test_create_from_cpu_count_default(self):
        """Should create based on CPU count."""
        settings = WorkerSettings.create_from_cpu_count()

        cpu_count = multiprocessing.cpu_count()
        expected_workers = max(1, int(cpu_count * 0.75 / 2))

        assert settings.n_workers == expected_workers
        assert settings.threads_per_worker == 2

    def test_create_from_cpu_count_custom_fraction(self):
        """Should accept custom CPU usage fraction."""
        settings = WorkerSettings.create_from_cpu_count(use_fraction=0.5)

        cpu_count = multiprocessing.cpu_count()
        expected_workers = max(1, int(cpu_count * 0.5 / 2))

        assert settings.n_workers == expected_workers

    def test_create_from_cpu_count_custom_threads(self):
        """Should accept custom threads per worker."""
        settings = WorkerSettings.create_from_cpu_count(threads_per_worker=4)

        assert settings.threads_per_worker == 4

    def test_create_from_cpu_count_custom_blocks(self):
        """Should accept custom block shape."""
        settings = WorkerSettings.create_from_cpu_count(block_shape=(512, 512))

        assert settings.block_shape == (512, 512)


class TestWorkerSettingsSystemValidation:
    """Tests for validate_against_system method."""

    def test_validate_against_system_normal(self):
        """Should pass validation for reasonable settings."""
        settings = WorkerSettings(n_workers=2, threads_per_worker=2)

        validation = settings.validate_against_system()

        assert "valid" in validation
        assert "warnings" in validation
        assert "system_cpu_count" in validation

    def test_warns_on_thread_oversubscription(self):
        """Should warn when total threads exceed CPU count."""
        cpu_count = multiprocessing.cpu_count()

        settings = WorkerSettings(
            n_workers=cpu_count,
            threads_per_worker=4,  # Likely to exceed
        )

        validation = settings.validate_against_system()

        if settings.total_threads > cpu_count:
            assert len(validation["warnings"]) > 0
            assert any("oversubscription" in w.lower() for w in validation["warnings"])

    def test_warns_on_large_blocks(self):
        """Should warn about very large block sizes."""
        settings = WorkerSettings(block_shape=(2000, 2000))  # 4M elements

        validation = settings.validate_against_system()

        assert any("very large" in w for w in validation["warnings"])

    def test_warns_on_small_blocks(self):
        """Should warn about very small block sizes."""
        settings = WorkerSettings(block_shape=(10, 10))  # 100 elements

        validation = settings.validate_against_system()

        assert any("very small" in w for w in validation["warnings"])

    def test_reports_cpu_utilization(self):
        """Should report CPU utilization percentage."""
        settings = WorkerSettings(n_workers=2, threads_per_worker=1)

        validation = settings.validate_against_system()

        assert "cpu_utilization" in validation
        assert 0 <= validation["cpu_utilization"] <= 2.0  # Could exceed 1.0


class TestWorkerSettingsSummary:
    """Tests for summary method."""

    def test_summary_contains_key_info(self):
        """Should include all key settings in summary."""
        settings = WorkerSettings(
            n_workers=8,
            threads_per_worker=4,
            block_shape=(256, 256),
        )

        summary = settings.summary()

        assert "8" in summary  # n_workers
        assert "4" in summary  # threads_per_worker
        assert "32" in summary  # total threads
        assert "256 x 256" in summary  # block shape

    def test_summary_is_multiline(self):
        """Should return multi-line string."""
        settings = WorkerSettings()

        summary = settings.summary()

        assert "\n" in summary
        assert len(summary.split("\n")) > 5


class TestWorkerSettingsPydantic:
    """Tests for Pydantic-specific features."""

    def test_forbids_extra_fields(self):
        """Should forbid extra fields."""
        with pytest.raises(ValueError):
            WorkerSettings(
                n_workers=4,
                extra_field="not allowed",  # type: ignore
            )

    def test_validates_on_assignment(self):
        """Should validate when fields are modified."""
        settings = WorkerSettings()

        with pytest.raises(ValueError):
            settings.n_workers = 0

        with pytest.raises(ValueError):
            settings.block_shape = (0, 128)  # type: ignore

    def test_serialization(self):
        """Should serialize to dict correctly."""
        settings = WorkerSettings(
            n_workers=8,
            threads_per_worker=4,
            block_shape=(256, 256),
        )

        data = settings.model_dump()

        assert data["n_workers"] == 8
        assert data["threads_per_worker"] == 4
        assert data["block_shape"] == (256, 256)

    def test_json_schema_examples(self):
        """Should have examples in JSON schema."""
        schema = WorkerSettings.model_json_schema()

        assert "examples" in schema
        assert len(schema["examples"]) > 0


@pytest.mark.parametrize(
    "n_workers,threads,expected_total",
    [
        (1, 1, 1),
        (2, 2, 4),
        (4, 4, 16),
        (8, 2, 16),
    ],
)
def test_total_threads_calculation(n_workers: int, threads: int, expected_total: int):
    """Should correctly calculate total threads for various configurations."""
    settings = WorkerSettings(n_workers=n_workers, threads_per_worker=threads)

    assert settings.total_threads == expected_total


@pytest.mark.parametrize(
    "block_shape,expected_size",
    [
        ((64, 64), 4096),
        ((128, 128), 16384),
        ((256, 256), 65536),
        ((100, 200), 20000),
    ],
)
def test_block_size_calculation(block_shape: tuple[int, int], expected_size: int):
    """Should correctly calculate block size for various shapes."""
    settings = WorkerSettings(block_shape=block_shape)

    assert settings.block_size == expected_size
