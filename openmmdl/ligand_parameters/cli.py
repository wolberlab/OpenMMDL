"""`openmmdl ligand`: fit or import bespoke ligand parameters, and restore them after tleap."""

from __future__ import annotations

import argparse

from openmmdl.ligand_parameters.core import FitOptions, prepare_parameters


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="openmmdl ligand",
        description="Bespoke ligand parameters with OpenFF BespokeFit for the OpenMM and Amber paths.",
    )
    commands = parser.add_subparsers(dest="command", metavar="<command>", required=True)

    fit = commands.add_parser("fit", help="Fit bespoke torsions for a ligand with BespokeFit")
    fit.add_argument("--executable", default="openff-bespoke", help="openff-bespoke executable (default: from PATH)")
    fit.add_argument(
        "--qc-spec",
        nargs=3,
        default=["psi4", "b3lyp-d3bj", "dzvp"],
        metavar=("PROGRAM", "METHOD", "BASIS"),
        help="QC level for the torsion scans (default: psi4 b3lyp-d3bj dzvp; e.g. xtb gfn2xtb none for a quick test)",
    )
    fit.add_argument("--workers", type=int, default=1, help="QC compute workers (default: 1)")
    fit.add_argument("--cores", type=int, default=2, help="Cores per QC worker (default: 2)")
    fit.add_argument("--memory-gb", type=int, default=4, help="Memory per QC worker in GB (default: 4)")

    import_ = commands.add_parser("import", help="Import an existing SMIRNOFF force field for a ligand")
    import_.add_argument("--offxml", required=True, help="SMIRNOFF force field (.offxml) covering the ligand")

    for command in (fit, import_):
        command.add_argument("--ligand", required=True, help="Ligand SDF/MOL file with explicit hydrogens")
        command.add_argument("--output", required=True, help="Output directory for the parameters")
        command.add_argument(
            "--target",
            choices=["openmm", "amber"],
            default="openmm",
            help="Also write a tleap library for the Amber path (default: openmm)",
        )

    restore = commands.add_parser("restore", help="Put the exact ligand parameters back into a tleap topology")
    restore.add_argument("--prmtop", required=True, help="Topology written by tleap (replaced in place)")
    restore.add_argument("--inpcrd", required=True, help="Coordinates written by tleap (replaced in place)")
    restore.add_argument("--ligand-dir", required=True, help="Output directory of 'openmmdl ligand fit/import --target amber'")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "restore":
        from openmmdl.ligand_parameters.amber import restore_ligand_parameters

        restore_ligand_parameters(args.prmtop, args.inpcrd, args.ligand_dir)
        print(f"Restored the bespoke ligand parameters in {args.prmtop}")
        return 0
    if args.command == "fit":
        program, method, basis = args.qc_spec
        options = FitOptions(
            executable=args.executable,
            qc_program=program,
            qc_method=method,
            qc_basis=basis,
            workers=args.workers,
            cores=args.cores,
            memory_gb=args.memory_gb,
        )
        result = prepare_parameters(args.ligand, args.output, options=options, target=args.target)
    else:
        result = prepare_parameters(args.ligand, args.output, offxml=args.offxml, target=args.target)
    print(f"Ligand force field ready: {result}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
