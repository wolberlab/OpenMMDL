from __future__ import annotations

import argparse
import importlib
import sys
from typing import Dict, List

COMMANDS: Dict[str, tuple[str, str]] = {
    "setup": (
        "openmmdl.openmmdl_setup.openmmdlsetup:main",
        "Start the OpenMMDL setup UI (prepare inputs for simulation)",
    ),
    "simulation": (
        "openmmdl.openmmdl_simulation.openmmdlsimulation:main",
        "Run OpenMM protein-ligand MD simulation",
    ),
    "analysis": (
        "openmmdl.openmmdl_analysis.openmmdlanalysis:main",
        "Analyze an OpenMMDL MD trajectory",
    ),
    "visualization": (
        "openmmdl.openmmdl_analysis.visualization.visualization:run_visualization",
        "Launch the visualization notebook",
    ),
}


HELP_ALIASES = {"--help", "-h", "help", "-help"}


def _build_top_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="openmmdl",
        description="OpenMMDL command-line interface.",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # We add subparsers only so argparse prints a clean command list in help.
    sub = parser.add_subparsers(dest="command", metavar="<command>")
    for name, (_, short_help) in COMMANDS.items():
        sub.add_parser(name, help=short_help)

    return parser


def _normalize_help_tokens(argv: List[str]) -> List[str]:
    """
    Normalize 'help' and '-help' to '--help'
    """
    if not argv:
        return argv

    # Top-level help aliases
    if argv[0] in {"help", "-help"}:
        return ["--help", *argv[1:]]

    # Command-level help aliases (2nd token)
    if len(argv) >= 2 and argv[1] in {"help", "-help"}:
        return [argv[0], "--help", *argv[2:]]

    return argv


def _run_target(target: str, forwarded: List[str], prog: str) -> int:
    """
    Import target and call either:
      - module:function     (callable)
      - module              (expects main())
    Forward argv via sys.argv.
    """
    old_argv = sys.argv
    sys.argv = [prog, *forwarded]
    try:
        if ":" in target:
            mod_name, func_name = target.split(":", 1)
            mod = importlib.import_module(mod_name)
            fn = getattr(mod, func_name, None)
            if fn is None or not callable(fn):
                raise RuntimeError(f"Module '{mod_name}' does not expose callable '{func_name}'.")
            rc = fn()
            return 0 if rc is None else int(rc)

        mod = importlib.import_module(target)
        if not hasattr(mod, "main"):
            raise RuntimeError(
                f"Module '{target}' does not expose a main() function. "
                "Add def main(argv=None) and call it from __main__."
            )
        rc = mod.main()  # type: ignore[attr-defined]
        return 0 if rc is None else int(rc)
    finally:
        sys.argv = old_argv


def main(argv: List[str] | None = None) -> int:
    """
    Dispatch: openmmdl <command> [args...]
    """
    if argv is None:
        argv = sys.argv[1:]

    argv = _normalize_help_tokens(argv)

    # Top-level help
    if not argv or argv[0] in HELP_ALIASES:
        _build_top_parser().print_help()
        return 0

    command = argv[0]
    if command not in COMMANDS:
        sys.stderr.write(f"Unknown command: {command!r}\n\n")
        _build_top_parser().print_help()
        return 2

    target, _ = COMMANDS[command]
    forwarded = argv[1:]

    # If user typed: openmmdl <command> --help, it will be forwarded
    # and the underlying argparse will handle it.
    prog = f"openmmdl {command}"
    try:
        return _run_target(target, forwarded, prog=prog)
    except SystemExit as e:
        # Preserve argparse exit behavior from downstream CLIs
        return int(e.code) if isinstance(e.code, int) else 0


if __name__ == "__main__":
    raise SystemExit(main())
