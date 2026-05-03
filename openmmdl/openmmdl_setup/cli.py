from __future__ import annotations

import argparse


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="openmmdl setup",
        description="Start the OpenMMDL setup UI.",
    )
    parser.add_argument(
        "--host",
        default="127.0.0.1",
        help="Host interface for the setup web server. Default: 127.0.0.1",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=5000,
        help="Port for the setup web server. Default: 5000",
    )
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="Do not automatically open the setup UI in a browser.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    from openmmdl.openmmdl_setup.openmmdlsetup import main as run_setup

    run_setup(
        host=args.host,
        port=args.port,
        open_browser=not args.no_browser,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())