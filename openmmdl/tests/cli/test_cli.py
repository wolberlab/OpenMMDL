import openmmdl.cli.cli as cli


def test_main_without_args_prints_help_and_returns_0(capsys):
    rc = cli.main([])

    captured = capsys.readouterr()
    output = captured.out + captured.err

    assert rc == 0
    assert "OpenMMDL command-line interface." in output
    assert "analysis" in output
    assert "simulation" in output
    assert "setup" in output
    assert "visualization" in output


def test_unknown_command_returns_2_and_prints_help(capsys):
    rc = cli.main(["not-a-real-command"])

    captured = capsys.readouterr()
    stderr = captured.err
    stdout = captured.out

    assert rc == 2
    assert "Unknown command: 'not-a-real-command'" in stderr
    assert "OpenMMDL command-line interface." in stdout


def test_version_flag_prints_version_and_returns_0(capsys):
    rc = cli.main(["--version"])

    captured = capsys.readouterr()
    output = captured.out.strip()

    assert rc == 0
    assert output == cli.__version__


def test_top_level_help_alias_returns_0(capsys):
    rc = cli.main(["help"])

    captured = capsys.readouterr()
    output = captured.out + captured.err

    assert rc == 0
    assert "OpenMMDL command-line interface." in output


def test_analysis_command_forwards_arguments_to_analysis_entrypoint(monkeypatch):
    called = {}

    def fake_run_target(target, forwarded, prog):
        called["target"] = target
        called["forwarded"] = forwarded
        called["prog"] = prog
        return 0

    monkeypatch.setattr(cli, "_run_target", fake_run_target)

    rc = cli.main(["analysis", "--verbose", "-t", "top.pdb", "-d", "traj.dcd", "-n", "UNK"])

    assert rc == 0
    assert called["target"] == cli.COMMANDS["analysis"][0]
    assert called["forwarded"] == ["--verbose", "-t", "top.pdb", "-d", "traj.dcd", "-n", "UNK"]
    assert called["prog"] == "openmmdl analysis"


def test_command_help_alias_is_forwarded_to_subcommand(monkeypatch):
    called = {}

    def fake_run_target(target, forwarded, prog):
        called["target"] = target
        called["forwarded"] = forwarded
        called["prog"] = prog
        return 0

    monkeypatch.setattr(cli, "_run_target", fake_run_target)

    rc = cli.main(["analysis", "help"])

    assert rc == 0
    assert called["target"] == cli.COMMANDS["analysis"][0]
    assert called["forwarded"] == ["--help"]
    assert called["prog"] == "openmmdl analysis"


def test_systemexit_from_subcommand_is_converted_to_return_code(monkeypatch):
    def fake_run_target(target, forwarded, prog):
        raise SystemExit(2)

    monkeypatch.setattr(cli, "_run_target", fake_run_target)

    rc = cli.main(["analysis", "--bad-arg"])

    assert rc == 2