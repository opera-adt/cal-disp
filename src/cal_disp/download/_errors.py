"""Exceptions raised by the staging (download) helpers."""

from __future__ import annotations

__all__ = ["DownloadError"]


class DownloadError(RuntimeError):
    """A download finished without a usable file.

    Raised when a response is not the requested data (an HTML error page,
    a body shorter than its ``Content-Length``) or when a search that must
    return something (e.g. tropospheric scenes for a sensing time) is empty.
    """
