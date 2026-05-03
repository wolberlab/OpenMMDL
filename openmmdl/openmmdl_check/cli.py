from __future__ import annotations

import argparse
import importlib.metadata
import importlib.util
import platform
import shutil
import sys
from dataclasses import dataclass

from openmmdl import __version__


@dataclass(frozen=True)
class CheckItem:
    name: str
    module: str
    package: str | None = None
    required: bool = True


CHECKS = [
    CheckItem("OpenMM", "openmm", "openmm", required=True),
    CheckItem("MDAnalysis", "MDAnalysis", "MDAnalysis", required=True),
    CheckItem("RDKit", "rdkit", "rdkit", required=True),
    CheckItem("MDTraj", "mdtraj", "mdtraj", required=True),
    CheckItem("PDBFixer", "pdbfixer", "pdbfixer", required=True),
    CheckItem("Flask", "flask", "flask", required=True),
    CheckItem("CairoSVG", "cairosvg", "CairoSVG", required=True),
    CheckItem("ProLIF", "prolif", "prolif", required=True),
    CheckItem("PLIP", "plip", "plip", required=False),
    CheckItem("NGLView", "nglview", "nglview", required=False),
    CheckItem("OpenFF Toolkit", "openff.toolkit", "openff-toolkit", required=False),
    CheckItem("OpenMM Force Fields", "openmmforcefields", "openmmforcefields", required=False),
]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="openmmdl check",
        description="Check the OpenMMDL installation and runtime environment.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Return a non-zero exit code if a required dependency is missing.",
    )
    return parser


def _version_for(package_name: str | None) -> str:
    if package_name is None:
        return "unknown"

    try:
        return importlib.metadata.version(package_name)
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _module_available(module_name: str) -> bool:
    return importlib.util.find_spec(module_name) is not None


def _format_status(available: bool) -> str:
    return "OK" if available else "MISSING"


def _print_header() -> None:
    print("OpenMMDL environment check")
    print("==========================")
    print(f"OpenMMDL: {__version__}")
    print(f"Python:   {platform.python_version()} ({sys.executable})")
    print(f"Platform: {platform.platform()}")
    print()


def _print_dependency_checks() -> bool:
    print("Dependencies")
    print("------------")

    missing_required = False

    for item in CHECKS:
        available = _module_available(item.module)
        version = _version_for(item.package) if available else "-"
        kind = "required" if item.required else "optional"

        if item.required and not available:
            missing_required = True

        print(f"{_format_status(available):<8} {item.name:<22} {version:<15} {kind}")

    print()
    return missing_required


def _print_command_checks() -> None:
    print("Commands")
    print("--------")

    for command in ["jupyter"]:
        path = shutil.which(command)
        status = "OK" if path else "MISSING"
        print(f"{status:<8} {command:<22} {path or '-'}")

    print()


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    _print_header()
    missing_required = _print_dependency_checks()
    _print_command_checks()

    if missing_required:
        print("Some required dependencies are missing.")
        if args.strict:
            return 1
    else:
        print("OpenMMDL check completed successfully.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())