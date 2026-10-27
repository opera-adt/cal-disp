"""Data staging helpers (the ``download`` extra).

Submodules that need ``asf_search``/``opera_utils`` are imported lazily, so
the core package (``cal_disp.product`` imports ``_stage_unr``) works without
the extra installed.
"""

from __future__ import annotations

from importlib import import_module

# NOTE add option to use s3 paths

_LAZY_ATTRS = {
    "download_disp": "._stage_disp",
    "download_tropo": "._stage_tropo",
    "download_unr_grid": "._stage_unr",
    "generate_s1_burst_tiles": "._stage_burst_bounds",
    "utils": ".utils",
}

__all__ = sorted(_LAZY_ATTRS)


def __getattr__(name: str):  # noqa: ANN202
    """Import the requested submodule/function on first access."""
    if name not in _LAZY_ATTRS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module = import_module(_LAZY_ATTRS[name], __name__)
    return module if name == "utils" else getattr(module, name)


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(_LAZY_ATTRS))
