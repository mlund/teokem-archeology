#!/usr/bin/env python3
"""
Unit test comparison between legacy reference and optimized output.

Compares density profiles (fort.83, fort.85) and key parameters from stdout.
"""

import numpy as np
import re
import sys

def extract_params(filename):
    """Extract key parameters from stdout file."""
    with open(filename) as f:
        content = f.read()

    params = {}
    patterns = {
        'polymer_chempot': r'polymer chemical pot\. \(betamu\) =\s+([\d.E+-]+)',
        'bulk_pressure': r'total bulk pressure =\s+([\d.E+-]+)',
        'bFex': r'bFex =\s+([\d.E+-]+)',
        'bebelam': r'bebelam,behbclam =\s+([\d.E+-]+)\s+([\d.E+-]+)',
        'iterations': r'ddmax,niter =\s+([\d.E+-]+)\s+(\d+)',
        'aW': r'aW =\s+([-\d.E+-]+)',
        'bW': r'bW =\s+([-\d.E+-]+)',
        'rcliffF': r'rcliffF =\s+([-\d.E+-]+)',
        'ctF': r'ctF =\s+([-\d.E+-]+)',
    }

    for key, pattern in patterns.items():
        match = re.search(pattern, content)
        if match:
            if key == 'iterations':
                all_matches = re.findall(pattern, content)
                params[key] = (float(all_matches[-1][0]), int(all_matches[-1][1]))
            elif key == 'bebelam':
                params['bebelam'] = float(match.group(1))
                params['behbclam'] = float(match.group(2))
            else:
                params[key] = float(match.group(1))

    return params

def compare_profiles(legacy_file, opt_file, profile_name):
    """Compare density profile files."""
    legacy = np.loadtxt(legacy_file)
    opt = np.loadtxt(opt_file)

    print(f"\n=== {profile_name} ===")
    print(f"Shape: legacy={legacy.shape}, optimized={opt.shape}")

    if legacy.shape != opt.shape:
        print("ERROR: Different shapes!")
        return False

    # Compare each column
    for col in range(1, legacy.shape[1]):  # Skip first column (coordinate)
        diff = np.abs(legacy[:, col] - opt[:, col])
        nonzero = legacy[:, col] > 1e-10

        if not np.any(nonzero):
            print(f"Column {col}: All zeros in legacy")
            continue

        rel_diff = diff[nonzero] / legacy[nonzero, col] * 100

        print(f"Column {col}:")
        print(f"  Max absolute diff: {np.max(diff):.6e}")
        print(f"  Mean absolute diff: {np.mean(diff[nonzero]):.6e}")
        print(f"  Max relative diff: {np.max(rel_diff):.3f}%")
        print(f"  Mean relative diff: {np.mean(rel_diff):.3f}%")

    return True

def main():
    # Compare stdout parameters
    print("=" * 70)
    print("UNIT TEST: Optimized vs Legacy Reference")
    print("=" * 70)

    try:
        legacy_params = extract_params('test_legacy/stdout')
        opt_params = extract_params('stdout_optimized.txt')
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        return 1

    print(f"\n{'Parameter':<25} {'Legacy':>20} {'Optimized':>20} {'Diff %':>12}")
    print("=" * 80)

    test_passed = True
    tolerance = {
        'polymer_chempot': 0.01,  # 0.01% tolerance for bulk properties
        'bulk_pressure': 0.01,
        'bFex': 0.01,
        'bebelam': 0.01,
        'behbclam': 0.01,
        'aW': 5.0,  # 5% tolerance for interaction free energy
        'bW': 1.0,   # 1% tolerance
        'rcliffF': 1.0,  # 1% tolerance for forces
        'ctF': 1.0,
    }

    for key in sorted(legacy_params.keys()):
        if key == 'iterations':
            l_val, l_iter = legacy_params[key]
            o_val, o_iter = opt_params[key]
            print(f"{'Convergence (ddmax)':<25} {l_val:>20.6e} {o_val:>20.6e}")
            print(f"{'Iterations':<25} {l_iter:>20} {o_iter:>20}")
        else:
            l_val = legacy_params[key]
            o_val = opt_params[key]

            if abs(l_val) > 1e-10:
                diff_rel = abs((o_val - l_val) / l_val * 100)
                status = "PASS" if diff_rel < tolerance.get(key, 1.0) else "FAIL"
                print(f"{key:<25} {l_val:>20.12e} {o_val:>20.12e} {diff_rel:>11.2f}% {status}")
                if status == "FAIL":
                    test_passed = False
            else:
                print(f"{key:<25} {l_val:>20.12e} {o_val:>20.12e}")

    # Compare density profiles
    compare_profiles('test_legacy/fort.83', 'fort.83', 'fort.83 (radial profile at zc1)')
    compare_profiles('test_legacy/fort.85', 'fort.85', 'fort.85 (axial profile at rho=0)')

    print("\n" + "=" * 70)
    if test_passed:
        print("RESULT: PASSED")
        return 0
    else:
        print("RESULT: FAILED - Significant differences detected")
        return 1

if __name__ == '__main__':
    sys.exit(main())
