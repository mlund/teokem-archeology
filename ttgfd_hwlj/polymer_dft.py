#!/usr/bin/env python3
"""
Run script for ttgfd_hwlj - Generalized Flory-Dimer Theory for Polymer Solutions
Generates input.tsph parameter file and executes the simulation.
"""

import argparse
import subprocess
import sys
from pathlib import Path


def write_input_file(args, output_path="input.tsph"):
    """Write the input.tsph parameter file."""
    kread_val = 1 if args.kread else 0

    # Write epsilon file
    with open("epfil", 'w') as f:
        f.write(f"{args.epslj}\n")

    with open(output_path, 'w') as f:
        f.write(f"  {args.bdm:22.15E}\n")      # Monomer bulk density
        f.write(f"  {args.bdtot:22.15f}\n")    # Total bulk density
        f.write(f"  {args.nmon:15d}\n")        # Number of monomers per polymer chain
        f.write(f"  {args.dz:22.15f}\n")       # Grid spacing in z direction
        f.write(f"  {args.drho:22.15f}\n")     # Grid spacing in radial direction
        f.write(f"  {args.dphi:22.15E}\n")     # Angular grid spacing (in units of pi)
        f.write(f"  {args.Rcoll:22.15f}\n")    # Colloid radius
        f.write(f"  {args.zc1:22.15f}\n")      # Position of first colloid center
        f.write(f"  {args.collsep:22.15f}\n")  # Separation between colloid centers
        f.write(f"  {args.Rcyl:22.15f}\n")     # Cylinder radius (system boundary)
        f.write(f"  {args.ioimaxm:10d}\n")     # Maximum number of iterations
        f.write(f"  {args.dmm:22.15f}       {args.dms:22.15f}\n")  # Density mixing params
        f.write(f"  {kread_val:15d}\n")        # Read initial guess from file (0=no, 1=yes)
        f.write(f"  {args.bl:22.15f}\n")       # Bond length
        f.write(f"  {args.dhs:22.15f}\n")      # Hard sphere diameter (monomer)
        f.write(f"  {args.dpphi:22.15E}\n")    # Angular grid spacing for potential calculation
        f.write("\n")


def main():
    parser = argparse.ArgumentParser(
        description='Generate input file and run ttgfd_hwlj polymer solution simulation',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Physical parameters
    parser.add_argument('--bdm', type=float, default=1.0e-2,
                        help='Monomer bulk density')
    parser.add_argument('--bdtot', type=float, default=1.0,
                        help='Total bulk density')
    parser.add_argument('--nmon', type=int, default=151,
                        help='Number of monomers per polymer chain')

    # Grid parameters
    parser.add_argument('--dz', type=float, default=0.25,
                        help='Grid spacing in z direction')
    parser.add_argument('--drho', type=float, default=0.25,
                        help='Grid spacing in radial direction')
    parser.add_argument('--dphi', type=float, default=5.0e-2,
                        help='Angular grid spacing (in units of pi)')
    parser.add_argument('--dpphi', type=float, default=2.5e-3,
                        help='Angular grid spacing for potential calculation')

    # Geometry parameters
    parser.add_argument('--Rcoll', type=float, default=10.0,
                        help='Colloid radius')
    parser.add_argument('--zc1', type=float, default=35.0,
                        help='Position of first colloid center')
    parser.add_argument('--collsep', type=float, default=20.0,
                        help='Separation between colloid centers')
    parser.add_argument('--Rcyl', type=float, default=35.0,
                        help='Cylinder radius (system boundary)')

    # Molecular parameters
    parser.add_argument('--bl', type=float, default=1.0,
                        help='Bond length')
    parser.add_argument('--dhs', type=float, default=1.0,
                        help='Hard sphere diameter (monomer)')
    parser.add_argument('--epslj', type=float, default=1.0,
                        help='Lennard-Jones epsilon parameter (energy scale)')

    # Numerical parameters
    parser.add_argument('--ioimaxm', type=int, default=5000,
                        help='Maximum number of iterations')
    parser.add_argument('--dmm', type=float, default=0.9,
                        help='Density mixing parameter (monomer)')
    parser.add_argument('--dms', type=float, default=0.5,
                        help='Density mixing parameter (solvent)')
    parser.add_argument('--kread', action='store_true',
                        help='Read initial guess from file fcdfil')

    # Script options
    parser.add_argument('--input-file', type=str, default='input.tsph',
                        help='Output path for input parameter file')
    parser.add_argument('--executable', type=str, default='./ttgfd_hwlj',
                        help='Path to ttgfd_hwlj executable')
    parser.add_argument('--dry-run', action='store_true',
                        help='Only generate input file, do not run simulation')

    args = parser.parse_args()

    # Write input file
    print(f"Generating {args.input_file}...")
    write_input_file(args, args.input_file)
    print(f"✓ Input file generated")

    if args.dry_run:
        print("Dry run mode: skipping execution")
        return 0

    # Check if executable exists
    exe_path = Path(args.executable)
    if not exe_path.exists():
        print(f"Error: Executable not found: {args.executable}", file=sys.stderr)
        return 1

    # Run simulation
    print(f"\nRunning {args.executable}...")
    print("=" * 80)
    try:
        result = subprocess.run([args.executable], check=True)
        print("=" * 80)
        print("✓ Simulation completed successfully")
        return result.returncode
    except subprocess.CalledProcessError as e:
        print("=" * 80)
        print(f"✗ Simulation failed with exit code {e.returncode}", file=sys.stderr)
        return e.returncode
    except KeyboardInterrupt:
        print("\n✗ Simulation interrupted by user", file=sys.stderr)
        return 130


if __name__ == '__main__':
    sys.exit(main())
