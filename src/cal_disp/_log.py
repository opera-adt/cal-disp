import logging
from pathlib import Path

from cal_disp._types import PathOrStr


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
    # Create or get logger
    logger = logging.getLogger(logger_name)
    logger.setLevel(getattr(logging, level.upper()))
    logger.handlers.clear()  # Remove existing handlers

    # Console handler
    console = logging.StreamHandler()
    console.setLevel(getattr(logging, level.upper()))
    fmt = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    console.setFormatter(fmt)
    logger.addHandler(console)

    # Optional file handler
    if filename is not None:
        Path(filename).parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(filename)
        file_handler.setLevel(logging.DEBUG)  # Always capture DEBUG to file
        file_handler.setFormatter(fmt)
        logger.addHandler(file_handler)

    # Suppress noisy third-party loggers
    for lib in ["asf_search", "botocore", "boto3", "urllib3", "s3fs"]:
        logging.getLogger(lib).setLevel(logging.WARNING)

    return logger
