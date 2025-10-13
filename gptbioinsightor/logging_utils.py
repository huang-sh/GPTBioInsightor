import logging
import sys
from typing import Optional


def _configure_default_logger(name: str = "gptbioinsightor") -> logging.Logger:
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stdout)
        handler.set_name("gptbioinsightor-default")
        formatter = logging.Formatter("[GPTBioInsightor] %(levelname)s - %(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        logger.propagate = False
    return logger


logger = _configure_default_logger()


def configure_logging(
    level: Optional[int] = None,
    handler: Optional[logging.Handler] = None,
    fmt: Optional[str] = None,
) -> logging.Logger:
    """Allow users to customize the packaged logger."""
    if level is not None:
        logger.setLevel(level)
    if handler is not None:
        # Remove default handler if a custom handler is provided
        logger.handlers = []
        logger.addHandler(handler)
    if fmt is not None:
        for handler in logger.handlers:
            handler.setFormatter(logging.Formatter(fmt))
    return logger
