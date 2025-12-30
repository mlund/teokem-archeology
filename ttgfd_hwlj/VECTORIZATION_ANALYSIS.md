# Vectorization Analysis - ttgfd_hwlj

## Summary

Analyzed compiler auto-vectorization on the 3 main performance hotspots:
1. **Chain propagation loop** (41% CPU - lines 470-574)
2. **EBDU subroutine** (15.5% CPU - lines 1537-1570)  
3. **hvec computation** (3.1% CPU - lines 357-387)

Compiler: gfortran with `-O3 -ftree-vectorize -march=native -ffast-math`

## Hotspot 1: Chain Propagation Loop (41% CPU)

**Location**: `ttgfd_hwlj_varbl.f90:470-574` (merged parallel region)

### Vectorization Status:

```fortran
! Line 470: Outer loop over iz (z-grid)
!$omp do schedule(static)
do iz = istp1 + ibl, imitt
  ...
  ! Line 476: Middle loop over kz (rho-grid)  
  do kz = irho0min, mxrho - ibl
    ...
    ! Line 499: Bond propagation loop over jz
    do jz = jstart, iz + ibl
      ...
      ! Line 508: Angular integration (CRITICAL INNERMOST LOOP)
!$omp simd reduction(+:phisum)
      do iphi = 1, nphi
        rho2 = rhoz2 - fphi*cos_phi(iphi)
        rsq = rho2 + zpcsq
        if (rsq .lt. Rcoll2) cycle  ! Inside first colloid
        rsq = rho2 + zpc2sq
        if (rsq .lt. Rcoll2) cycle  ! Inside second colloid
        rho = dsqrt(rho2)
        irho = int(rho*rdrho) + 1
        phisum = cA(irho, jz) + phisum
      end do
!$omp end simd
```

**Compiler Report:**
```
Line 470: missed: not vectorized: loop nest containing two or more consecutive inner loops
Line 476: missed: not vectorized: loop nest containing two or more consecutive inner loops  
Line 499: missed: not vectorized: loop nest containing two or more consecutive inner loops
Line 508: optimized: loop vectorized using 16 byte vectors
Line 511: missed: not vectorized: unsupported control flow in loop
```

### Analysis:

✅ **SUCCESS**: Line 508 (angular integration) DOES vectorize despite the `!$omp simd` directive
- Uses 16-byte vectors (2x double precision values per SIMD operation)
- This is the innermost computational kernel with ~400 iterations

❌ **FAILURE**: Line 511 reports "unsupported control flow" but line 508 reports success
- **Root cause**: The `cycle` statements (conditional exits) prevent full vectorization
- Compiler creates a masked/predicated version that partially works

❌ **MISSED**: Outer loops (iz, kz, jz) cannot vectorize due to nested loop structure

### Performance Impact:

The angular integration loop (line 508) is the **most compute-intensive** part:
- 150 monomers × ~45 z-grid points × ~32 rho-grid points × ~2 jz-grid points × **400 phi angles**
- Total: ~17 million angular integrations per iteration × 300 iterations = **5 billion operations**

Vectorizing this loop with 2-wide SIMD provides **up to 2x speedup on this kernel**.

### Optimization Opportunities:

🔴 **HIGH PRIORITY**: Eliminate conditional exits (`cycle` statements) to enable full vectorization
```fortran
! Current (prevents vectorization):
if (rsq .lt. Rcoll2) cycle

! Better approach:
mask = merge(1.0d0, 0.0d0, rsq >= Rcoll2)  
phisum = phisum + mask * cA(irho, jz)
```

🟡 **MEDIUM**: The `cycle` statements cause branch mispredictions and break SIMD pipelines

---

## Hotspot 2: EBDU Subroutine (15.5% CPU)

**Location**: `ttgfd_hwlj_varbl.f90:1537-1570`

### Vectorization Status:

```fortran
! Line 1537: Outer loop over z-grid
!$omp parallel do private(...)
do iz = istp1, imitt
  ...
  ! Line 1540: Loop over rho-grid
  do krho = 1, mxrho
    ...
    ! Line 1550: Angular integration (CRITICAL INNERMOST LOOP - transposed access)
!$omp simd reduction(+:sumrho)
    do kprho = 1, mxrho
      sumrho = sumrho + (fdmon(kprho, ipz) - bdm)*hvec(kprho, krho, itdz)
    end do
!$omp end simd
```

**Compiler Report:**
```
Line 1540: missed: not vectorized: loop nest containing two or more consecutive inner loops
Line 1550: optimized: loop vectorized using 16 byte vectors
Line 1562: optimized: loop vectorized using 16 byte vectors  
```

### Analysis:

✅ **SUCCESS**: Lines 1550 and 1562 vectorize successfully
- Uses 16-byte vectors (2x double precision)
- **Stride-1 memory access** thanks to hvec transposition optimization
- No conditional exits, clean reduction operation

✅ **ALREADY OPTIMIZED**: hvec transposition provides 15% speedup (from previous work)
- Original: `hvec(itdz, krho, kprho)` - large stride
- Optimized: `hvec(kprho, krho, itdz)` - stride-1 innermost access

### Performance Impact:

This subroutine is **well-optimized** and vectorizing correctly. The 15% speedup from transposition + vectorization is already captured.

---

## Hotspot 3: hvec Computation (3.1% CPU)

**Location**: `ttgfd_hwlj_varbl.f90:357-387`

### Vectorization Status:

```fortran
! Line 359: Outer loop over tdz
!$omp parallel do private(...)
do itdz = 0, nfack
  ...
  ! Line 362: Loop over kprho
  do kprho = 0, mxrho
    ...
    ! Line 367: Loop over krho  
    do krho = kprho, mxrho
      ...
      ! Line 374: Angular integration (CRITICAL INNERMOST LOOP)
!$omp simd reduction(+:pint)
      do iphi = 1, nphi
        s2 = trhosq - trmix*cos_phi(iphi)
        useful = s2 + tdzsq
        use1 = useful**(-3.d0)
        pint = pint + (use1*use1 - use1)
      end do
!$omp end simd
```

**Compiler Report:**
```
Line 359: missed: not vectorized: loop nest containing two or more consecutive inner loops
Line 362: missed: not vectorized: loop nest containing two or more consecutive inner loops
Line 367: missed: not vectorized: loop nest containing two or more consecutive inner loops  
Line 374: optimized: loop vectorized using 16 byte vectors
Line 376: optimized: loop vectorized using 16 byte vectors
```

### Analysis:

✅ **SUCCESS**: Lines 374 and 376 vectorize successfully
- Uses 16-byte vectors (2x double precision)
- No conditional exits, clean arithmetic

**Only executes once per run** (not in iteration loop), so 3.1% overhead is acceptable.

---

## Overall Vectorization Summary

| Hotspot | CPU % | Innermost Loop | Vectorized? | Speedup Potential |
|---------|-------|----------------|-------------|-------------------|
| Chain propagation | 41% | Line 508 (phisum) | ✅ Partial (masked) | 🔴 **HIGH**: Fix cycle statements |
| EBDU | 15.5% | Line 1550 (sumrho) | ✅ Full | Already optimized |
| hvec | 3.1% | Line 374 (pint) | ✅ Full | Already optimized |

## Key Findings

### ✅ Good News:
1. **All critical innermost loops are vectorizing** (lines 508, 1550, 1562, 374, 376)
2. Compiler uses 16-byte vectors (2-wide SIMD for double precision)
3. `!$omp simd` directives are working correctly

### ⚠️ Problem:
**Chain propagation angular loop** (line 508-520) has **conditional exits that prevent full vectorization**:
```fortran
if (rsq .lt. Rcoll2) cycle  ! Breaks SIMD pipeline
```

The compiler creates a "masked" version that partially vectorizes but:
- Generates predicated instructions (slower than pure SIMD)
- May cause branch mispredictions
- Doesn't achieve full 2x speedup potential

## Recommended Optimizations

### Priority 1: Eliminate Conditional Exits in Chain Propagation Loop

**Current code (line 508-520):**
```fortran
!$omp simd reduction(+:phisum)
do iphi = 1, nphi
  rho2 = rhoz2 - fphi*cos_phi(iphi)
  rsq = rho2 + zpcsq
  if (rsq .lt. Rcoll2) cycle  ! BREAKS VECTORIZATION
  rsq = rho2 + zpc2sq
  if (rsq .lt. Rcoll2) cycle  ! BREAKS VECTORIZATION
  rho = dsqrt(rho2)
  irho = int(rho*rdrho) + 1
  phisum = cA(irho, jz) + phisum
end do
!$omp end simd
```

**Optimized approach (mask-based):**
```fortran
!$omp simd reduction(+:phisum)
do iphi = 1, nphi
  rho2 = rhoz2 - fphi*cos_phi(iphi)
  rsq1 = rho2 + zpcsq
  rsq2 = rho2 + zpc2sq
  
  ! Create mask: 1.0 if valid, 0.0 if inside colloid
  valid = merge(1.0d0, 0.0d0, rsq1 >= Rcoll2 .and. rsq2 >= Rcoll2)
  
  rho = dsqrt(rho2)
  irho = int(rho*rdrho) + 1
  phisum = phisum + valid * cA(irho, jz)
end do
!$omp end simd
```

**Expected Impact:**
- Remove branch mispredictions (~5-10% speedup)
- Enable full SIMD utilization (additional 1.2-1.5x on this loop)
- Combined: **1.3-1.6x speedup on 41% hotspot = ~12-25% total speedup**

### Priority 2: Check Array Bounds for Autovectorization Hints

Add `!GCC$ ivdep` or restrict pointers to help compiler prove no aliasing:
```fortran
!GCC$ ivdep  ! Ignore vector dependencies
!$omp simd reduction(+:phisum)
do iphi = 1, nphi
  ...
end do
```

### Priority 3: Experiment with AVX-512 (if available)

Current: 16-byte vectors (2x double = SSE/AVX2)
Potential: 64-byte vectors (8x double = AVX-512)

Add compiler flag: `-mavx512f` (if CPU supports it)

---

## Testing Recommendations

1. **Verify current vectorization** with Intel Advisor or similar:
   ```bash
   gfortran -O3 -ftree-vectorize -march=native -fopt-info-vec-optimized ttgfd_hwlj_varbl.f90
   ```

2. **Measure baseline performance** of vectorized code

3. **Implement Priority 1 optimization** (eliminate cycle statements)

4. **Benchmark improvement** - expect 10-25% total speedup

5. **Profile with hardware counters** to confirm SIMD utilization:
   ```bash
   perf stat -e fp_arith_inst_retired.scalar_double,fp_arith_inst_retired.128b_packed_double ./ttgfd_hwlj
   ```

