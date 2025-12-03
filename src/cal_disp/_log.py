import logging
from logging import Formatter
from pathlib import Path

from cal_disp._types import PathOrStr


def setup_file_logging(filename: PathOrStr) -> None:
    """Redirect all logging to a file."""
    logger = logging.getLogger()
    # In addition to stderr, log to a file if requested
    Path(filename).parent.mkdir(parents=True, exist_ok=True)
    file_handler = logging.FileHandler(filename)
    file_handler.setLevel(logging.DEBUG)
    formatter = Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
