# Remaining Vectorization Opportunities

**Date**: 2025-12-31
**Analysis**: Comprehensive scan for missed SIMD vectorization due to `cycle` statements and control flow

---

## Summary

After optimizing the main chain propagation loop (41% of runtime), there are still a few remaining `cycle` statements that prevent vectorization in minor subroutines. However, **the impact is negligible** since these account for less than 1% of total runtime.

**Performance Profile Context**:
- Chain propagation (already optimized): 41.1% ✅
- Thread synchronization: 37.1% (reduced via parallel region merging) ✅
- EBDU: 15.5% (already vectorizing) ✅
- **EBLMNEW**: 0.2% ⚠️ Has cycle statements
- **CDCALC**: 0.2% ⚠️ Indirect indexing prevents vectorization
- Other: <1%

**Conclusion**: Remaining optimizations would yield <0.5% speedup. Not worth the effort.

---

## Detailed Analysis

### 1. EBLMNEW Subroutine - Angular Loop (0.2% runtime)

**Location**: `ttgfd_hwlj_varbl.f90:1458-1473`

**Status**: ❌ NOT vectorizing due to `cycle` statements

**Code**:
```fortran
!$omp simd reduction(+:phisum)
do iphi = 1, nphi
  rho2 = rhomax2 - fphi*cos_phi(iphi)
  ! Skip if integration point is inside first colloid
  rsq = rho2 + zpcsq
  if (rsq .lt. Rcoll2) cycle     ! ← Line 1464: Prevents vectorization
  ! Skip if integration point is inside second colloid
  rsq = rho2 + zpc2sq
  if (rsq .lt. Rcoll2) cycle     ! ← Line 1467: Prevents vectorization
  rho = dsqrt(rho2)
  irho = int(rho*rdrho) + 1
  phisum = convp(irho, jz) + phisum
end do
!$omp end simd
```

**Compiler Report**:
```
Line 1458: missed: couldn't vectorize loop
Line 1458: missed: not vectorized: unsupported control flow in loop
```

**Fix**: Apply same mask-based approach as in chain propagation:
```fortran
!$omp simd reduction(+:phisum)
do iphi = 1, nphi
  rho2 = rhomax2 - fphi*cos_phi(iphi)
  rsq1 = rho2 + zpcsq
  rsq2 = rho2 + zpc2sq
  valid = merge(1.0d0, 0.0d0, rsq1 >= Rcoll2 .and. rsq2 >= Rcoll2)
  rho = dsqrt(rho2)
  irho = int(rho*rdrho) + 1
  phisum = phisum + valid * convp(irho, jz)
end do
!$omp end simd
```

**Expected speedup**: ~0.1% (EBLMNEW is only 0.2% of runtime, mask overhead ~50%)

**Recommendation**: ⚠️ **Skip** - not worth the effort for 0.1% gain

---

### 2. EBLMNEW Subroutine - Outer Loops

**Location**: `ttgfd_hwlj_varbl.f90:1413-1430`

**Status**: ❌ NOT vectorizing (but not performance-critical)

**Code**:
```fortran
do kz = irho0min, mxrho
  rho0 = rho0 + drho
  rho02 = rho0**2

  ! Skip points inside first colloid
  rt2 = rho02 + (z - zc1)**2
  if (rt2 .lt. Rcoll2) then
    ehbclam(kz, iz) = 0.d0
    ebelam(kz, iz) = 0.d0
    cycle                        ! ← Line 1422
  end if
  ! Skip points inside second colloid
  rt2 = rho02 + (z - zc2)**2
  if (rt2 .lt. Rcoll2) then
    ehbclam(kz, iz) = 0.d0
    ebelam(kz, iz) = 0.d0
    cycle                        ! ← Line 1429
  end if
  ...
```

**Compiler Report**:
```
Line 1413: missed: couldn't vectorize loop
Line 1413: missed: not vectorized: unsupported control flow in loop
```

**Fix**: Same mask-based approach
**Expected speedup**: <0.1%
**Recommendation**: ⚠️ **Skip**

---

### 3. Chain Propagation - Outer Loop (41% runtime, already optimized)

**Location**: `ttgfd_hwlj_varbl.f90:476-493`

**Status**: ⚠️ Outer loop has `cycle`, but **innermost loop already vectorizing** ✅

**Code**:
```fortran
do kz = irho0min, mxrho - ibl          ! ← Outer loop (NOT the critical one)
  rho0 = rho0 + drho
  rho02 = rho0**2
  rt2 = rho02 + (z - zc1)**2
  if (rt2 .lt. Rcoll2) then
    c(kz, iz, imon) = 0.d0
    cB(kz, iz) = 0.d0
    cycle                                 ! ← Line 485
  end if
  rt2 = rho02 + (z - zc2)**2
  if (rt2 .lt. Rcoll2) then
    c(kz, iz, imon) = 0.d0
    cB(kz, iz) = 0.d0
    cycle                                 ! ← Line 492
  end if

  ! Inner loops below (jz, then iphi)
  do jz = jstart, iz + ibl
    ...
    !$omp simd reduction(+:phisum)      ! ← This loop IS vectorizing ✅
    do iphi = 1, nphi
      ! Already optimized with mask-based approach
```

**Compiler Report**:
```
Line 476: missed: couldn't vectorize loop (expected - outer loop)
Line 517: optimized: loop vectorized using 16 byte vectors  ✅
```

**Analysis**:
- The cycle statements are in the **middle loop (kz)**, not the **innermost loop (iphi)**
- The iphi loop (line 517) is the performance-critical one and IS vectorizing
- Vectorizing the kz loop would provide minimal benefit (it's the middle loop, not tightly nested)

**Recommendation**: ✅ **Already optimal** - no action needed

---

### 4. CDCALC Subroutine - Angular Loop (0.2% runtime)

**Location**: `ttgfd_hwlj_varbl.f90:1300-1307`

**Status**: ❌ NOT vectorizing due to **indirect array indexing** (not cycle statements)

**Code**:
```fortran
!$omp simd reduction(+:phisum)
do iphi = 1, nphi
  rho2 = rhomax2 - fphi*cos_phi(iphi)
  rho = dsqrt(rho2)
  irho = int(rho*rdrho) + 1          ! ← Indirect index
  phisum = fdmon(irho, jz) + phisum  ! ← Gather operation
end do
!$omp end simd
```

**Compiler Report**:
```
Line 1300: missed: couldn't vectorize loop
Line 1300: missed: may need non-SLP handling
```

**Analysis**:
- The issue is NOT `cycle` statements (there are none)
- The issue is **gather operation**: `fdmon(irho, jz)` where `irho` varies unpredictably
- Modern CPUs support vector gather, but gfortran may not auto-vectorize this pattern
- CDCALC is only 0.2% of runtime (52 samples out of 28,500)

**Potential Fix**: Manual SIMD with gather intrinsics (very complex)
**Expected speedup**: <0.1%
**Recommendation**: ⚠️ **Skip** - gather overhead may negate any benefit

---

## Summary of Remaining cycle Statements

| Location | Line | Subroutine | % Runtime | Impact if Fixed | Recommendation |
|----------|------|------------|-----------|-----------------|----------------|
| iphi loop | 1464, 1467 | EBLMNEW | 0.2% | ~0.1% | Skip |
| kz loop | 1422, 1429 | EBLMNEW | 0.2% | <0.1% | Skip |
| kz loop | 485, 492 | Chain prop | 41% | None* | Already optimal |

*The cycle statements are in the middle loop; the critical innermost loop (iphi) is already vectorizing

---

## Performance Impact Analysis

### Current Status After Optimizations

| Optimization | Speedup | Cumulative |
|--------------|---------|------------|
| Baseline | - | 14.55s |
| Merge parallel regions | 8.0% | 13.39s |
| Eliminate cycle in chain prop iphi loop | 5.6% | 12.69s |
| **Total** | **14.7%** | **12.69s** |

### If We Fixed Remaining Issues

| Additional Optimization | Est. Speedup | New Time | Effort |
|-------------------------|--------------|----------|--------|
| Fix EBLMNEW iphi cycle | 0.1% | 12.68s | Low |
| Fix EBLMNEW kz cycle | <0.1% | 12.68s | Low |
| Fix CDCALC gather | <0.1% | 12.68s | Very High |
| **Total additional** | **~0.2%** | **~12.67s** | Medium-High |

**Conclusion**: Fixing remaining vectorization issues would gain <0.2% speedup for medium-high effort. **Not recommended.**

---

## Recommendations

### ✅ Keep Current State
The major vectorization opportunities have been captured:
1. Chain propagation innermost loop (41% of runtime): ✅ Vectorizing
2. EBDU LJ convolution (15.5% of runtime): ✅ Vectorizing
3. hvec tabulation (3.1% of runtime): ✅ Vectorizing

**Total vectorized code: >60% of runtime**

### ⚠️ Skip Remaining Optimizations
- EBLMNEW and CDCALC together account for 0.4% of runtime
- Fixing them would yield <0.2% total speedup
- Better to focus on:
  - **Algorithmic improvements** (adaptive mixing, better convergence)
  - **Cache optimization** (loop tiling, blocking)
  - **AVX-512** on capable hardware (30-80% speedup potential)

### 📊 Alternative High-Impact Optimizations

If seeking further performance, consider these instead:

1. **AVX-512 (Intel/AMD only)**: 1.3-1.8x speedup on 512-bit SIMD
   - Already enabled with `-march=native`
   - Test on Intel Xeon or AMD Zen4+

2. **Cache Optimization**: 5-15% potential
   - Loop tiling for better L1/L2 cache usage
   - Array layout optimization (AoS → SoA)

3. **GPU Port**: 5-20x potential (for large systems)
   - OpenACC or CUDA Fortran
   - Most effective for systems with >10M grid points

---

## Conclusion

**Current state: Excellent** ✅
- All major hotspots (>60% of runtime) are vectorizing
- 14.7% speedup achieved from vectorization + parallel region merging

**Remaining opportunities: Negligible** ⚠️
- <0.2% total speedup from fixing remaining `cycle` statements
- Not worth the development/testing effort

**Recommendation**: Declare vectorization optimization **complete** and focus on other performance avenues (cache, AVX-512, or algorithmic improvements).
