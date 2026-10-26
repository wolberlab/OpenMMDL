from __future__ import annotations

import argparse


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="openmmdl analysis-gui",
        description=(
            "Start the OpenMMDL Analysis GUI: configure and run an analysis, "
            "or explore the results folder of a finished one."
        ),
    )
    parser.add_argument(
        "results",
        nargs="?",
        default=None,
        help="Optional results folder of a finished analysis to open right away.",
    )
    parser.add_argument(
        "--host",
        default="127.0.0.1",
        help="Host interface for the GUI web server. Default: 127.0.0.1",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=5050,
        help="Port for the GUI web server. Default: 5050",
    )
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="Do not open the GUI in a browser automatically.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    from openmmdl.openmmdl_analysis_gui.server import main as run_server

    run_server(
        host=args.host,
        port=args.port,
        open_browser=not args.no_browser,
        results_folder=args.results,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
