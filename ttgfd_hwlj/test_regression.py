#!/usr/bin/env python3
"""
Regression test for optimized ttgfd_hwlj.

Tests:
1. Fresh start (kread=0) produces expected forces
2. Restart (kread=1) converges from own fcdfil
3. Self-consistency between fresh and restart
4. Density profile consistency between fresh and restart
"""

import numpy as np
import re
import sys
import subprocess
import os
import shutil

def extract_param(content, pattern):
    """Extract a single parameter from output."""
    match = re.search(pattern, content)
    return float(match.group(1)) if match else None

def extract_final_iteration(content):
    """Extract final iteration data."""
    pattern = r'ddmax,niter =\s+([\d.E+-]+)\s+(\d+)'
    matches = re.findall(pattern, content)
    if matches:
        ddmax, niter = matches[-1]
        return float(ddmax), int(niter)
    return None, None

def compare_density_profile(file1, file2, name):
    """Compare two density profile files."""
    try:
        data1 = np.loadtxt(file1)
        data2 = np.loadtxt(file2)

        if data1.shape != data2.shape:
            return False, float('inf'), f"Shape mismatch: {data1.shape} vs {data2.shape}"

        # Compare all columns except the first (coordinate column)
        if len(data1.shape) == 1:
            # Single column file
            diff = np.abs(data1 - data2)
            max_diff = np.max(diff)
            max_rel_diff = max_diff / (np.max(np.abs(data1)) + 1e-16) * 100
        else:
            # Multi-column file - skip first column (coordinates)
            data_cols1 = data1[:, 1:]
            data_cols2 = data2[:, 1:]

            diff = np.abs(data_cols1 - data_cols2)
            max_diff = np.max(diff)
            max_rel_diff = max_diff / (np.max(np.abs(data_cols1)) + 1e-16) * 100

        passed = max_rel_diff < 1.0  # Pass if < 1% difference
        return passed, max_rel_diff, None

    except Exception as e:
        return False, float('inf'), str(e)

def run_test(kread, expected_force=None, max_iters=200):
    """Run the test with given kread value."""
    # Create input file
    input_content = f""" 0.10000000000000001     
  1.0000000000000000     
        151
 0.25000000000000000     
 0.25000000000000000     
  5.0000000000000003E-002
  10.000000000000000     
  30.000000000000000     
  44.000000000000000     
  30.000000000000000     
       5000
 0.90000000000000002       0.50000000000000000     
          {kread}
  1.0000000000000000     
  1.0000000000000000     
  2.5000000000000001E-003

"""
    with open('input.tsph', 'w') as f:
        f.write(input_content)
    
    # Run the program
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = '4'
    
    try:
        result = subprocess.run(
            ['./ttgfd_hwlj'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=300,
            env=env,
            text=True
        )
        
        if result.returncode != 0:
            print(f"ERROR: Program exited with code {result.returncode}")
            print(result.stderr)
            return None
            
        output = result.stdout
        
        # Extract parameters
        rcliffF = extract_param(output, r'rcliffF =\s+([-\d.E+-]+)')
        ctF = extract_param(output, r'ctF =\s+([-\d.E+-]+)')
        ch2 = extract_param(output, r'ch2 =\s+([-\d.E+-]+)')
        aW = extract_param(output, r'aW =\s+([-\d.E+-]+)')
        bW = extract_param(output, r'bW =\s+([-\d.E+-]+)')
        ddmax, niter = extract_final_iteration(output)

        # Save density profile files with kread suffix
        profile_files = ['fort.78', 'fort.83', 'fort.85', 'fort.87', 'fort.89']
        for f in profile_files:
            if os.path.exists(f):
                shutil.copy(f, f'{f}.kread{kread}')

        return {
            'rcliffF': rcliffF,
            'ctF': ctF,
            'ch2': ch2,
            'aW': aW,
            'bW': bW,
            'ddmax': ddmax,
            'iterations': niter,
            'output': output
        }
        
    except subprocess.TimeoutExpired:
        print(f"ERROR: Test with kread={kread} timed out (>300s)")
        return None
    except Exception as e:
        print(f"ERROR: {e}")
        return None

def main():
    print("=" * 70)
    print("REGRESSION TEST: ttgfd_hwlj Optimized Version")
    print("=" * 70)
    
    # Ensure epfil exists
    if not os.path.exists('test_legacy/epfil'):
        print("ERROR: test_legacy/epfil not found")
        return 1
    
    subprocess.run(['cp', 'test_legacy/epfil', '.'], check=True)
    
    all_passed = True
    
    # Test 1: Fresh start (kread=0)
    print("\nTEST 1: Fresh start (kread=0)")
    print("-" * 70)
    
    if os.path.exists('fcdfil'):
        os.remove('fcdfil')
    
    result_kread0 = run_test(kread=0)
    
    if result_kread0 is None:
        print("FAIL: Test did not complete")
        all_passed = False
    else:
        print(f"  Iterations:  {result_kread0['iterations']}")
        print(f"  rcliffF:     {result_kread0['rcliffF']:.15e}")
        print(f"  ctF:         {result_kread0['ctF']:.15e}")
        print(f"  aW:          {result_kread0['aW']:.15e}")
        print(f"  bW:          {result_kread0['bW']:.15e}")
        print(f"  Final ddmax: {result_kread0['ddmax']:.6e}")
        
        # Check convergence (with smart adaptive mixing, expect ~71 iterations)
        if result_kread0['iterations'] > 100:
            print(f"  WARNING: Took {result_kread0['iterations']} iterations (expected ~71)")
        
        # Expected values from verified kread=0 run with smart adaptive mixing
        expected_rcliffF = -2.3072255938144659E-005
        expected_aW = -50.812621716068243
        expected_bW = -5125.2051000404972

        # Allow tiny numerical differences from compiler optimizations and adaptive mixing
        tol_force = 1e-9  # Absolute tolerance for forces (adaptive mixing path dependent)
        tol_energy = 1e-5  # Absolute tolerance for energies (allows ~1e-6 variation)
        
        force_diff = abs(result_kread0['rcliffF'] - expected_rcliffF)
        aW_diff = abs(result_kread0['aW'] - expected_aW)
        bW_diff = abs(result_kread0['bW'] - expected_bW)
        
        if force_diff < tol_force and aW_diff < tol_energy and bW_diff < tol_energy:
            print("  PASS: Forces and energies match expected values ✓")
        else:
            print(f"  FAIL: Differences detected:")
            print(f"    rcliffF diff: {force_diff:.6e} (tol: {tol_force:.6e})")
            print(f"    aW diff:      {aW_diff:.6e} (tol: {tol_energy:.6e})")
            print(f"    bW diff:      {bW_diff:.6e} (tol: {tol_energy:.6e})")
            all_passed = False
    
    # Test 2: Restart (kread=1) from own fcdfil
    print("\nTEST 2: Restart (kread=1) from own fcdfil")
    print("-" * 70)
    
    if not os.path.exists('fcdfil'):
        print("SKIP: No fcdfil from Test 1")
        all_passed = False
    else:
        result_kread1 = run_test(kread=1)
        
        if result_kread1 is None:
            print("FAIL: Test did not complete")
            all_passed = False
        else:
            print(f"  Iterations:  {result_kread1['iterations']}")
            print(f"  rcliffF:     {result_kread1['rcliffF']:.15e}")
            print(f"  ctF:         {result_kread1['ctF']:.15e}")
            print(f"  aW:          {result_kread1['aW']:.15e}")
            print(f"  bW:          {result_kread1['bW']:.15e}")
            print(f"  Final ddmax: {result_kread1['ddmax']:.6e}")
            
            # Check convergence (should converge, not oscillate)
            if result_kread1['iterations'] > 200:
                print(f"  FAIL: Too many iterations ({result_kread1['iterations']})")
                print("        May be oscillating instead of converging")
                all_passed = False
            else:
                print(f"  PASS: Converged in {result_kread1['iterations']} iterations ✓")
            
            # Test 3: Self-consistency
            print("\nTEST 3: Self-consistency check")
            print("-" * 70)
            
            if result_kread0 and result_kread1:
                force_diff = abs(result_kread1['rcliffF'] - result_kread0['rcliffF'])
                rel_diff = force_diff / abs(result_kread0['rcliffF']) * 100
                
                print(f"  kread=0 force: {result_kread0['rcliffF']:.15e}")
                print(f"  kread=1 force: {result_kread1['rcliffF']:.15e}")
                print(f"  Difference:    {force_diff:.15e}")
                print(f"  Relative:      {rel_diff:.4f}%")
                
                if rel_diff < 0.01:  # Less than 0.01% difference
                    print("  PASS: Self-consistent ✓")
                else:
                    print(f"  FAIL: {rel_diff:.4f}% difference (expected < 0.01%)")
                    all_passed = False

            # Test 4: Density profile consistency
            print("\nTEST 4: Density profile consistency")
            print("-" * 70)

            profile_files = [
                ('fort.78', 'Integrated density per z-slice'),
                ('fort.83', 'Radial profile at z=zc1'),
                ('fort.85', 'Axial profile at rho=0'),
                ('fort.87', 'Chain propagators at z=zc1'),
                ('fort.89', 'Propagators along centerline')
            ]

            profiles_passed = True
            for fname, description in profile_files:
                file_kread0 = f'{fname}.kread0'
                file_kread1 = f'{fname}.kread1'

                if not os.path.exists(file_kread0) or not os.path.exists(file_kread1):
                    print(f"  {fname}: SKIP (files not found)")
                    continue

                passed, max_rel_diff, error = compare_density_profile(file_kread0, file_kread1, fname)

                if error:
                    print(f"  {fname}: ERROR - {error}")
                    profiles_passed = False
                    all_passed = False
                elif passed:
                    print(f"  {fname}: PASS (max rel diff = {max_rel_diff:.4f}%) ✓")
                else:
                    print(f"  {fname}: FAIL (max rel diff = {max_rel_diff:.4f}% > 1.0%)")
                    profiles_passed = False
                    all_passed = False

            if profiles_passed:
                print("\n  RESULT: All density profiles consistent ✓")

    # Final result
    print("\n" + "=" * 70)
    if all_passed:
        print("RESULT: ALL TESTS PASSED ✓")
        print("=" * 70)
        return 0
    else:
        print("RESULT: SOME TESTS FAILED ✗")
        print("=" * 70)
        return 1

if __name__ == '__main__':
    sys.exit(main())
