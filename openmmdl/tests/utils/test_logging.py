import logging

import pytest

from openmmdl.utils.logging_utils import setup_logging


pytestmark = pytest.mark.serial


@pytest.fixture(autouse=True)
def reset_logging_state():
    """
    Keep tests isolated because setup_logging() reconfigures global logging.
    """
    root_logger = logging.getLogger()
    app_logger = logging.getLogger("openmmdl")

    old_root_handlers = list(root_logger.handlers)
    old_root_level = root_logger.level

    old_app_handlers = list(app_logger.handlers)
    old_app_level = app_logger.level
    old_app_propagate = app_logger.propagate

    yield

    for logger in (root_logger, app_logger):
        for handler in list(logger.handlers):
            try:
                handler.flush()
                handler.close()
            except Exception:
                pass
            logger.removeHandler(handler)

    for handler in old_root_handlers:
        root_logger.addHandler(handler)
    root_logger.setLevel(old_root_level)

    for handler in old_app_handlers:
        app_logger.addHandler(handler)
    app_logger.setLevel(old_app_level)
    app_logger.propagate = old_app_propagate


def _flush_all_handlers():
    for logger_name in ("", "openmmdl"):
        logger = logging.getLogger(logger_name)
        for handler in logger.handlers:
            handler.flush()


def test_setup_logging_creates_both_log_files_and_writes_openmmdl_messages(tmp_path, capsys):
    setup_logging(verbose=False, log_dir=tmp_path, log_prefix="openmmdl")

    logger = logging.getLogger("openmmdl.tests.analysis")
    logger.info("analysis started")

    _flush_all_handlers()

    normal_log = tmp_path / "openmmdl.log"
    verbose_log = tmp_path / "openmmdl_verbose.log"

    assert normal_log.exists()
    assert verbose_log.exists()

    normal_text = normal_log.read_text()
    verbose_text = verbose_log.read_text()

    captured = capsys.readouterr()
    console_text = captured.out + captured.err

    assert "INFO: analysis started" in normal_text
    assert "analysis started" in verbose_text
    assert "INFO: analysis started" in console_text


def test_dependency_info_is_hidden_from_console_and_normal_log_but_kept_in_verbose_log(tmp_path, capsys):
    setup_logging(verbose=False, log_dir=tmp_path, log_prefix="openmmdl")

    app_logger = logging.getLogger("openmmdl.tests.analysis")
    dependency_logger = logging.getLogger("plip.structure.preparation")

    app_logger.info("top level analysis message")
    dependency_logger.info("analyzing 9 ligands")
    dependency_logger.info("protonated structure written to /tmp/frame_1_protonated.pdb")

    _flush_all_handlers()

    normal_text = (tmp_path / "openmmdl.log").read_text()
    verbose_text = (tmp_path / "openmmdl_verbose.log").read_text()

    captured = capsys.readouterr()
    console_text = captured.out + captured.err

    assert "top level analysis message" in normal_text
    assert "top level analysis message" in verbose_text
    assert "top level analysis message" in console_text

    assert "analyzing 9 ligands" not in normal_text
    assert "protonated structure written" not in normal_text
    assert "analyzing 9 ligands" not in console_text
    assert "protonated structure written" not in console_text

    assert "analyzing 9 ligands" in verbose_text
    assert "protonated structure written to /tmp/frame_1_protonated.pdb" in verbose_text


def test_verbose_mode_shows_debug_and_logger_name_in_console(tmp_path, capsys):
    setup_logging(verbose=True, log_dir=tmp_path, log_prefix="openmmdl")

    logger = logging.getLogger("openmmdl.tests.analysis")
    logger.debug("debug details")

    _flush_all_handlers()

    normal_text = (tmp_path / "openmmdl.log").read_text()
    verbose_text = (tmp_path / "openmmdl_verbose.log").read_text()

    captured = capsys.readouterr()
    console_text = captured.out + captured.err

    assert "debug details" not in normal_text
    assert "debug details" in verbose_text
    assert "DEBUG:openmmdl.tests.analysis:debug details" in console_text


def test_normal_mode_hides_openmmdl_debug_from_console_and_normal_log_but_keeps_it_in_verbose_log(
    tmp_path, capsys
):
    setup_logging(verbose=False, log_dir=tmp_path, log_prefix="openmmdl")

    logger = logging.getLogger("openmmdl.tests.analysis")
    logger.debug("hidden debug")

    _flush_all_handlers()

    normal_text = (tmp_path / "openmmdl.log").read_text()
    verbose_text = (tmp_path / "openmmdl_verbose.log").read_text()

    captured = capsys.readouterr()
    console_text = captured.out + captured.err

    assert "hidden debug" not in normal_text
    assert "hidden debug" in verbose_text
    assert "hidden debug" not in console_text


def test_custom_log_prefix_changes_output_filenames(tmp_path):
    setup_logging(verbose=False, log_dir=tmp_path, log_prefix="openmmdl_simulation")

    logger = logging.getLogger("openmmdl.tests.simulation")
    logger.info("simulation started")

    _flush_all_handlers()

    assert (tmp_path / "openmmdl_simulation.log").exists()
    assert (tmp_path / "openmmdl_simulation_verbose.log").exists()
    assert "INFO: simulation started" in (tmp_path / "openmmdl_simulation.log").read_text()