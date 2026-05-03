from __future__ import annotations

import argparse


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="openmmdl visualization",
        description="Launch the OpenMMDL visualization notebook.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    parser.parse_args(argv)

    from openmmdl.openmmdl_analysis.visualization.visualization import run_visualization

    run_visualization()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())