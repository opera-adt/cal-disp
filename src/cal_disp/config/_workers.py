import multiprocessing
from typing import Any, Dict, Tuple

from pydantic import ConfigDict, Field, field_validator

from ._yaml import YamlModel


class WorkerSettings(YamlModel):
    """Settings for controlling CPU/threading and parallelism.

    This class configures Dask distributed computing settings including
    worker count, threads per worker, and data block sizes for chunked
    processing.

    Attributes
    ----------
    n_workers : int
        Number of Dask workers to spawn. Default is 4.
    threads_per_worker : int
        Number of threads per worker. Default is 2.
    block_shape : Tuple[int, int]
        Block size (rows, columns) for chunked data loading.

    Examples
    --------
    >>> # Default settings
    >>> settings = WorkerSettings()
    >>>
    >>> # Custom configuration
    >>> settings = WorkerSettings(
    ...     n_workers=8,
    ...     threads_per_worker=4,
    ...     block_shape=(256, 256)
    ... )
    >>> print(f"Total threads: {settings.total_threads}")

    """

    n_workers: int = Field(
        default=4,
        ge=1,
        le=256,  # Reasonable upper limit
        description=(
            "Number of workers to use in dask.Client. Typically set to the number "
            "of physical CPU cores or less. Default: 4"
        ),
    )

    threads_per_worker: int = Field(
        default=2,
        ge=1,
        le=32,  # Reasonable upper limit
        description=(
            "Number of threads per worker in dask.Client. Controls parallelism "
            "within each worker. Default: 2"
        ),
    )

    block_shape: Tuple[int, int] = Field(
        default=(128, 128),
        description=(
            "Size (rows, columns) of blocks of data to load at a time. "
            "Larger blocks use more memory but may be more efficient. "
            "Must be positive integers. Default: (128, 128)"
        ),
    )

    @field_validator("block_shape", mode="before")
    @classmethod
    def _validate_block_shape(cls, v) -> Tuple[int, int]:
        """Validate block shape has positive dimensions.

        Parameters
        ----------
        v : tuple | list
            Block shape to validate.

        Returns
        -------
        Tuple[int, int]
            Validated block shape.

        Raises
        ------
        ValueError
            If block dimensions are not positive or not exactly 2 values.

        """
        # Convert to tuple if list
        if isinstance(v, list):
            v = tuple(v)

        # Check it's a tuple/list with 2 elements
        if not isinstance(v, tuple) or len(v) != 2:
            raise ValueError(
                f"block_shape must have exactly 2 dimensions (rows, cols), got {v}"
            )

        # Check both are integers
        if not all(isinstance(x, int) for x in v):
            raise ValueError(f"block_shape dimensions must be integers, got {v}")

        # Check both are positive
        if not all(x > 0 for x in v):
            raise ValueError(f"block_shape dimensions must be positive, got {v}")

        return v

    # CHANGED: Use @property instead of @computed_field to avoid mypy issues
    @property
    def total_threads(self) -> int:
        """Total number of threads across all workers.

        Returns
        -------
        int
            n_workers * threads_per_worker

        """
        return self.n_workers * self.threads_per_worker

    @property
    def block_size(self) -> int:
        """Total number of elements per block.

        Returns
        -------
        int
            Product of block_shape dimensions (rows * columns).

        """
        return self.block_shape[0] * self.block_shape[1]

    def estimate_memory_per_block(self, dtype_size: int = 8, n_bands: int = 1) -> float:
        """Estimate memory usage per block in MB.

        Parameters
        ----------
        dtype_size : int, default=8
            Size of data type in bytes (e.g., 8 for float64, 4 for float32).
        n_bands : int, default=1
            Number of bands/layers in the data.

        Returns
        -------
        float
            Estimated memory in megabytes.

        Examples
        --------
        >>> settings = WorkerSettings(block_shape=(256, 256))
        >>> mem_mb = settings.estimate_memory_per_block(dtype_size=8, n_bands=2)
        >>> print(f"Estimated memory: {mem_mb:.2f} MB")

        """
        bytes_per_block = self.block_size * dtype_size * n_bands
        return bytes_per_block / (1024 * 1024)  # Convert to MB

    def estimate_total_memory(
        self, dtype_size: int = 8, n_bands: int = 1, overhead_factor: float = 1.5
    ) -> float:
        """Estimate total memory usage across all workers in GB.

        Parameters
        ----------
        dtype_size : int, default=8
            Size of data type in bytes.
        n_bands : int, default=1
            Number of bands/layers in the data.
        overhead_factor : float, default=1.5
            Multiplier for overhead (copies, intermediate results).

        Returns
        -------
        float
            Estimated total memory in gigabytes.

        """
        mb_per_block = self.estimate_memory_per_block(dtype_size, n_bands)
        # Assume each worker might hold multiple blocks
        total_mb = mb_per_block * self.n_workers * overhead_factor
        return total_mb / 1024  # Convert to GB

    def summary(self) -> str:
        """Generate a human-readable summary of settings.

        Returns
        -------
        str
            Multi-line summary string.

        """
        lines = [
            "WorkerSettings Configuration:",
            "=" * 50,
            f"Workers:              {self.n_workers}",
            f"Threads per worker:   {self.threads_per_worker}",
            f"Total threads:        {self.total_threads}",
            f"Block shape:          {self.block_shape[0]} x {self.block_shape[1]}",
            f"Block size:           {self.block_size:,} elements",
            f"Estimated mem/block:  {self.estimate_memory_per_block():.2f} MB",
            f"Estimated total mem:  {self.estimate_total_memory():.2f} GB",
            "=" * 50,
        ]
        return "\n".join(lines)

    @classmethod
    def create_lightweight(cls) -> "WorkerSettings":
        """Create lightweight settings for small datasets or testing.

        Returns
        -------
        WorkerSettings
            Configuration with 2 workers, 1 thread each, small blocks.

        """
        return cls(n_workers=2, threads_per_worker=1, block_shape=(64, 64))

    @classmethod
    def create_standard(cls) -> "WorkerSettings":
        """Create standard settings for typical workloads.

        Returns
        -------
        WorkerSettings
            Configuration with 4 workers, 2 threads each, medium blocks.

        """
        return cls(n_workers=4, threads_per_worker=2, block_shape=(128, 128))

    @classmethod
    def create_heavy(cls) -> "WorkerSettings":
        """Create heavy-duty settings for large datasets.

        Returns
        -------
        WorkerSettings
            Configuration with 8 workers, 4 threads each, large blocks.

        """
        return cls(n_workers=8, threads_per_worker=4, block_shape=(256, 256))

    @classmethod
    def create_from_cpu_count(
        cls,
        use_fraction: float = 0.75,
        threads_per_worker: int = 2,
        block_shape: Tuple[int, int] = (128, 128),
    ) -> "WorkerSettings":
        """Create settings based on available CPU count.

        Parameters
        ----------
        use_fraction : float, default=0.75
            Fraction of available CPUs to use (0.0 to 1.0).
        threads_per_worker : int, default=2
            Threads per worker.
        block_shape : Tuple[int, int], default=(128, 128)
            Block shape for data loading.

        Returns
        -------
        WorkerSettings
            Configuration tuned to system CPU count.

        Examples
        --------
        >>> # Use 75% of available CPUs
        >>> settings = WorkerSettings.create_from_cpu_count(use_fraction=0.75)

        """
        cpu_count = multiprocessing.cpu_count()
        n_workers = max(1, int(cpu_count * use_fraction / threads_per_worker))

        return cls(
            n_workers=n_workers,
            threads_per_worker=threads_per_worker,
            block_shape=block_shape,
        )

    def validate_against_system(self) -> Dict[str, Any]:
        """Validate settings against system resources.

        Returns
        -------
        dict
            Dictionary with validation results and warnings.

        Examples
        --------
        >>> settings = WorkerSettings(n_workers=16, threads_per_worker=8)
        >>> validation = settings.validate_against_system()
        >>> if validation['warnings']:
        ...     for warning in validation['warnings']:
        ...         print(f"Warning: {warning}")

        """
        cpu_count = multiprocessing.cpu_count()
        warnings = []

        # Check if total threads exceed CPU count
        if self.total_threads > cpu_count:
            warnings.append(
                f"Total threads ({self.total_threads}) exceeds available "
                f"CPU cores ({cpu_count}). May cause oversubscription."
            )

        # Check if block size is very large
        if self.block_size > 1_000_000:
            warnings.append(
                f"Block size ({self.block_size:,}) is very large. "
                "May cause high memory usage."
            )

        # Check if block size is very small
        if self.block_size < 1_000:
            warnings.append(
                f"Block size ({self.block_size:,}) is very small. "
                "May cause poor performance due to overhead."
            )

        return {
            "valid": len(warnings) == 0,
            "warnings": warnings,
            "system_cpu_count": cpu_count,
            "configured_total_threads": self.total_threads,
            "cpu_utilization": self.total_threads / cpu_count if cpu_count > 0 else 0,
        }

    model_config = ConfigDict(
        extra="forbid",
        validate_assignment=True,
        json_schema_extra={
            "examples": [
                {"n_workers": 4, "threads_per_worker": 2, "block_shape": [128, 128]},
                {"n_workers": 8, "threads_per_worker": 4, "block_shape": [256, 256]},
            ]
        },
    )
