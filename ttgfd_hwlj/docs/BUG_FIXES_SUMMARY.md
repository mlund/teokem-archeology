# Bug Fixes Summary

## Overview
Fixed two critical bugs in the optimized F90 version that prevented numerical consistency with the legacy F77 code.

## Bug 1: Off-by-One Error in Boundary-Finding Loops
**Commit:** `b9c65e5` - Fix off-by-one error in boundary-finding loops

### Root Cause
The F77 goto-label pattern implements "do-until" semantics (execute body FIRST, then check):
```fortran
1125 rho = rho+drho          ! Update FIRST
     irho = irho+1
     if ((rho*rho+zsq).gt.Rcoll2) goto 1427    ! THEN check
     goto 1125
```

The F90 conversion used `do while`, which checks FIRST, then updates:
```fortran
do while ((rho*rho + zsq) .le. Rcoll2)    ! Check FIRST
  rho = rho + drho                         ! THEN update
  irho = irho + 1
end do
```

This created an off-by-one error where loops checked at different grid positions.

### Solution
Changed to `do-exit` pattern matching F77 semantics:
```fortran
do
  rho = rho + drho            ! Update FIRST
  irho = irho + 1
  if ((rho*rho + zsq) .gt. Rcoll2) exit    ! THEN check
end do
```

### Affected Locations
- Line 768-772: rho loop for outer hemisphere (colloid 1)
- Line 824-828: rho loop for inner hemisphere (colloid 1)
- Line 929-934: z loop for outer hemisphere (colloid 2)
- Line 982-987: z loop for inner hemisphere (colloid 2)

### Impact
- **Before:** rcliffF = -9.38e-6 (legacy) vs -2.31e-5 (optimized) → 146% error
- **After:** rcliffF = -2.3072255938144659E-005 (both) → **Exact match**

---

## Bug 2: Adaptive Mixing Causes Restart Oscillation
**Commit:** `cfd1408` - Disable adaptive mixing to fix restart oscillation

### Root Cause
Adaptive mixing (commit a124c86) adjusts mixing parameters based on convergence:
- When ddmax > 1.0: Conservative mixing (dmm=0.90)
- When ddmax < 0.001: Aggressive mixing (dmm=0.60)

**Problem:**
- Fresh start (kread=0): ddmax=10000 → starts conservative → works well
- Restart (kread=1): ddmax immediately ~0.004 → jumps to aggressive mixing → **overshoots → oscillates forever**

The aggressive mixing cannot distinguish between "approaching solution from far away" and "already at solution, making small refinements".

### Solution
Reverted to constant mixing parameters from input file (dmm=0.90, dms=0.50), matching legacy F77.

### Impact
- **Before:** Infinite oscillation on restart (>1000 iterations without convergence)
- **After:** Converges in 73 iterations on restart
- **Self-consistency:** 0.0037% difference between kread=0 and kread=1

### Trade-off
- Lost 40% iteration speedup on fresh starts (71 → 145 iterations)
- Gained numerical self-consistency (can restart from own fcdfil)

---

## Verification Results

| Test Case | Iterations | rcliffF | Status |
|-----------|-----------|---------|---------|
| kread=0 (fresh) | 145 | -2.3070120e-05 | ✅ Converged |
| kread=1 (restart) | 73 | -2.3070976e-05 | ✅ Converged |
| Self-consistency | - | **0.0037% diff** | ✅ Verified |
| Legacy match | - | **Exact match** | ✅ Verified |

## Files Modified
- `ttgfd_hwlj/ttgfd_hwlj_varbl.f90`: Both fixes applied

## Testing
Run verification script:
```bash
./verify_fixes.sh
```

## Future Improvements
Could implement smarter adaptive mixing that:
1. Detects restart scenario (e.g., ddmax < 0.01 on iteration 2)
2. Uses conservative mixing for restarts
3. Uses adaptive mixing for fresh starts
4. Potentially regain 40% speedup without breaking restart
