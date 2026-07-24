from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
import logging

LOG_NAMESPACE = "af3_screening"


@contextmanager
def rule_logger(
    logfile: str | Path,
    *,
    rule_name: str,
    level: int = logging.INFO,
    mode: str = "w",
) -> Iterator[logging.Logger]:
    path = Path(logfile)
    path.parent.mkdir(parents=True, exist_ok=True)

    parent = logging.getLogger(LOG_NAMESPACE)

    previous_level = parent.level
    previous_propagate = parent.propagate

    handler = logging.FileHandler(
        path,
        mode=mode,
        encoding="utf-8",
    )
    handler.setLevel(level)
    handler.setFormatter(
        logging.Formatter(
            "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
    )

    parent.setLevel(level)
    parent.propagate = False
    parent.addHandler(handler)

    logger = logging.getLogger(f"{LOG_NAMESPACE}.rules.{rule_name}")

    try:
        yield logger
    except Exception:
        logger.exception("Rule failed")
        raise
    finally:
        parent.removeHandler(handler)
        handler.close()

        parent.setLevel(previous_level)
        parent.propagate = previous_propagate
