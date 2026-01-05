from __future__ import annotations

import logging
import resource
import sys
import time
from collections.abc import Callable
from functools import wraps
from pathlib import Path
from typing import TypeVar

if sys.version_info >= (3, 10):
    from typing import ParamSpec
else:
    from typing_extensions import ParamSpec

from cal_disp._types import PathOrStr

# Used for callable types
T = TypeVar("T")
P = ParamSpec("P")


def setup_logging(
    logger_name: str = "cal_disp",
    *,
    level: str = "INFO",
    filename: PathOrStr | None = None,
) -> logging.Logger:
    """Configure logging with optional file output.

    Sets up a logger with console output and optionally writes to a file.
    Automatically suppresses verbose third-party library logging.

    Parameters
    ----------
    logger_name : str, optional
        Name of the logger to configure. Default is "cal_disp".
    level : str, optional
        Logging level (DEBUG, INFO, WARNING, ERROR). Default is "INFO".
    filename : PathOrStr or None, optional
        If provided, also log to this file. Parent directories are created
        if they don't exist. Default is None (console only).

    Returns
    -------
    logging.Logger
        Configured logger instance.

    Examples
    --------
    Console logging only:
    >>> logger = setup_logging(level="DEBUG")

    With file output:
    >>> logger = setup_logging(filename="output.log")

    """
    logger = logging.getLogger(logger_name)
    logger.setLevel(getattr(logging, level.upper()))
    logger.handlers.clear()

    # Console handler
    console = logging.StreamHandler()
    console.setLevel(getattr(logging, level.upper()))
    fmt = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    console.setFormatter(fmt)
    logger.addHandler(console)

    # File handler
    if filename is not None:
        Path(filename).parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(filename)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(fmt)
        logger.addHandler(file_handler)

    # Suppress noisy third-party loggers
    for lib in ["asf_search", "botocore", "boto3", "urllib3", "s3fs"]:
        logging.getLogger(lib).setLevel(logging.WARNING)

    return logger


def log_runtime(f: Callable[P, T]) -> Callable[P, T]:
    """Decorate a function to time how long it takes to run.

    Examples
    --------
    >>> @log_runtime
    ... def test_func():
    ...     return 2 + 4

    """
    logger = logging.getLogger(__name__)

    @wraps(f)
    def wrapper(*args: P.args, **kwargs: P.kwargs) -> T:
        t1 = time.time()
        result = f(*args, **kwargs)
        t2 = time.time()

        elapsed_seconds = t2 - t1
        elapsed_minutes = elapsed_seconds / 60.0
        logger.debug(
            f"Total elapsed time for {f.__module__}.{f.__name__}: "
            f"{elapsed_minutes:.2f} minutes ({elapsed_seconds:.2f} seconds)"
        )

        return result

    return wrapper


def get_max_memory_usage(units: str = "GB", children: bool = True) -> float:
    """Get the maximum memory usage of the current process.

    Parameters
    ----------
    units : str, optional
        The units to return ("GB", "MB", "KB", "byte"). Default is "GB".
    children : bool, optional
        Whether to include child process memory. Default is True.

    Returns
    -------
    float
        The maximum memory usage in the specified units.

    Raises
    ------
    ValueError
        If the units are not recognized.

    References
    ----------
    .. [1] https://stackoverflow.com/a/7669279/4174466
    .. [2] https://unix.stackexchange.com/a/30941/295194
    .. [3] https://manpages.debian.org/bullseye/manpages-dev/getrusage.2.en.html

    """
    max_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if children:
        max_mem += resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss

    if units.lower().startswith("g"):
        factor = 1e9
    elif units.lower().startswith("m"):
        factor = 1e6
    elif units.lower().startswith("k"):
        factor = 1e3
    elif units.lower().startswith("byte"):
        factor = 1.0
    else:
        msg = f"Unknown units: {units}"
        raise ValueError(msg)

    # On Linux, ru_maxrss is in kilobytes; on macOS it's in bytes
    if sys.platform.startswith("linux"):
        factor /= 1e3

    return max_mem / factor
