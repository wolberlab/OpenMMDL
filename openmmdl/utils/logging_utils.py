import logging
import os
from pathlib import Path
from typing import Optional, Union


PathLike = Union[str, os.PathLike]


def setup_logging(
    verbose: bool = False,
    log_dir: Optional[PathLike] = None,
    app_logger_name: str = "openmmdl",
    log_prefix: str = "openmmdl",
    append: bool = False,
) -> None:
    """
    Configure shared OpenMMDL logging.

    Console
    -------
    - normal mode: INFO+ from OpenMMDL only
    - verbose mode: DEBUG+ from OpenMMDL with logger names

    Files
    -----
    - <log_prefix>.log:
        clean INFO+ OpenMMDL messages only
    - <log_prefix>_verbose.log:
        DEBUG+ messages from OpenMMDL and dependencies

    Notes
    -----
    - OpenMMDL loggers should use names under `app_logger_name`
      (for example, "openmmdl.openmmdl_analysis...").
    - Dependency INFO chatter is captured only in the verbose logfile,
      not in the console or the normal logfile.
    """
    if log_dir is None:
        log_dir = Path.cwd()
    else:
        log_dir = Path(log_dir)

    log_dir.mkdir(parents=True, exist_ok=True)

    file_mode = "a" if append else "w"

    normal_log_path = log_dir / f"{log_prefix}.log"
    verbose_log_path = log_dir / f"{log_prefix}_verbose.log"

    console_formatter = logging.Formatter(
        "%(levelname)s:%(name)s:%(message)s"
        if verbose
        else "%(levelname)s: %(message)s"
    )
    normal_file_formatter = logging.Formatter("%(levelname)s: %(message)s")
    verbose_file_formatter = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(name)s | %(processName)s | %(message)s"
    )

    root_logger = logging.getLogger()
    _clear_handlers(root_logger)
    root_logger.setLevel(logging.DEBUG)

    verbose_file_handler = logging.FileHandler(verbose_log_path, mode=file_mode)
    verbose_file_handler.setLevel(logging.DEBUG)
    verbose_file_handler.setFormatter(verbose_file_formatter)
    root_logger.addHandler(verbose_file_handler)

    app_logger = logging.getLogger(app_logger_name)
    _clear_handlers(app_logger)
    app_logger.setLevel(logging.DEBUG)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    console_handler.setFormatter(console_formatter)

    normal_file_handler = logging.FileHandler(normal_log_path, mode=file_mode)
    normal_file_handler.setLevel(logging.INFO)
    normal_file_handler.setFormatter(normal_file_formatter)

    app_logger.addHandler(console_handler)
    app_logger.addHandler(normal_file_handler)

    # Let OpenMMDL logs also go to root, so they appear in the verbose logfile.
    app_logger.propagate = True

    # Dependency logs should be captured in the verbose logfile.
    # They will not appear in console or normal log because those handlers
    # are attached only to the OpenMMDL logger.
    logging.getLogger("MDAnalysis").setLevel(logging.INFO)
    logging.getLogger("plip").setLevel(logging.INFO)
    logging.getLogger("prolif").setLevel(logging.INFO)
    logging.getLogger("matplotlib").setLevel(logging.WARNING)


def get_logger(name: str) -> logging.Logger:
    return logging.getLogger(name)


def _clear_handlers(logger: logging.Logger) -> None:
    for handler in list(logger.handlers):
        try:
            handler.flush()
            handler.close()
        except Exception:
            pass
        logger.removeHandler(handler)