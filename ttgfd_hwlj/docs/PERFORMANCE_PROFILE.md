# Performance Profile Analysis - ttgfd_hwlj

**Profile Date**: 2025-12-30
**Tool**: macOS `sample` (10-second sampling at 1ms interval)
**Parameters**: Fast test parameters (bdm=0.1, epslj=0.5, dz=1.0, drho=0.3)
**Threads**: OMP_NUM_THREADS=4

## Executive Summary

Total samples collected: ~28,500 over 10 seconds

The **single biggest performance bottleneck** is OpenMP thread synchronization overhead (37%), which exceeds the time spent in any single computational kernel. This presents the highest-priority optimization opportunity.

## Top Hotspots (by CPU time)

| Rank | Function | Samples | % CPU | Description |
|------|----------|---------|-------|-------------|
| 1 | `MAIN__._omp_fn.1` | 11,710 | 41.1% | Polymer chain propagation loop |
| 2 | `__psynch_cvwait` | 10,575 | 37.1% | **Thread synchronization overhead** |
| 3 | `ebdu_._omp_fn.1` | 4,429 | 15.5% | External potential (EBDU) calculation |
| 4 | `MAIN__._omp_fn.0` | 872 | 3.1% | Lennard-Jones hvec tabulation |
| 5 | `eblmnew_._omp_fn.1` | 67 | 0.2% | Chain end Boltzmann factors |
| 6 | `cdcalc_._omp_fn.0` | 52 | 0.2% | Contact density calculation |
| 7 | `MAIN__._omp_fn.3` | 18 | 0.1% | Density field mixing |
| - | Other | 777 | 2.7% | System calls, memory operations, etc. |

## Detailed Analysis

### 1. Chain Propagation Loop (41.1% - MAIN__._omp_fn.1)

**Location**: `ttgfd_hwlj_varbl.f90:1020-1154`

This is the primary computational kernel, performing polymer chain propagation through the density field. It includes:
- Triple nested loops over monomers, z-positions, and radial positions
- Bond propagation with Gaussian convolution
- Monomer-monomer excluded volume interactions

**Current Optimizations**:
- ✅ OpenMP parallelization over monomer index
- ✅ SIMD directives on inner loops
- ✅ Precomputed bond propagation factors

**Potential Optimizations**:
- Verify vectorization quality with compiler reports
- Consider loop tiling/blocking for cache locality
- Profile for branch mispredictions

### 2. Thread Synchronization Overhead (37.1% - Critical!)

**Function**: `__psynch_cvwait` (pthread condition variable waits)

This represents pure overhead from OpenMP barriers and thread synchronization.

**Root Cause**: Multiple separate parallel regions per iteration
- hvec computation (once per run)
- Chain propagation loop (MAIN__._omp_fn.1)
- EBDU calculation
- AVEC calculation
- CDCALC calculation
- EBLMNEW calculation
- Density mixing loop

Each `!$omp parallel` directive creates a barrier where all threads must synchronize.

**Impact**: With 4 threads, overhead from 6-7 parallel regions per iteration adds up to 37% of runtime.

### 3. EBDU - External Potential (15.5% - ebdu_._omp_fn.1)

**Location**: `ttgfd_hwlj_varbl.f90:1505-1562`

Calculates the external potential from Lennard-Jones interactions with the monomer density field.

**Current Optimizations**:
- ✅ hvec array transposition for cache-friendly access (+15% speedup)
- ✅ OpenMP parallelization over radial/z grid
- ✅ Stride-1 memory access in innermost loop

This subroutine is already well-optimized from the hvec transposition work.

### 4. hvec Computation (3.1% - MAIN__._omp_fn.0)

**Location**: `ttgfd_hwlj_varbl.f90:354-387`

Pre-computes the Lennard-Jones interaction potential lookup table.

**Current Optimizations**:
- ✅ Parallelized with OpenMP
- ✅ SIMD directives on angular integration
- ✅ Cosine lookup table to avoid repeated trig calls

Only executed once per run, so 3.1% overhead is acceptable.

## Optimization Recommendations

### Priority 1: Reduce Thread Synchronization (HIGHEST IMPACT)

**Potential Speedup**: 1.3-1.5x (by recovering 20-30% of the 37% overhead)

**Strategy**: Merge parallel regions to reduce barriers

**Specific Actions**:

1. **Merge subroutine calls into main parallel region**:
   ```fortran
   !$omp parallel private(...)
   do while (iteration_loop)
     ! All field calculations within same parallel region
     call CDCALC_PARALLEL  ! No barrier
     call AVEC_PARALLEL    ! No barrier
     call EBLMNEW_PARALLEL ! No barrier
     call EBDU_PARALLEL    ! No barrier
     ! Chain propagation
     ! Density mixing
   end do
   !$omp end parallel
   ```

2. **Use `!$omp single` for serial sections** instead of ending/restarting parallel regions

3. **Consider taskloop** for independent operations that don't need barriers

**Effort**: Medium (requires refactoring parallel structure)

### Priority 2: Optimize Chain Propagation Loop

**Potential Speedup**: 1.1-1.2x

**Actions**:
1. Generate compiler vectorization reports: `gfortran -fopt-info-vec-all`
2. Check if inner loops are vectorizing properly
3. Consider loop interchange if innermost loop has poor stride
4. Profile for branch mispredictions

**Effort**: Medium

### Priority 3: Cache Optimization via Tiling

**Potential Speedup**: 1.05-1.15x

**Actions**:
1. Profile L1/L2/L3 cache miss rates
2. Implement loop tiling on large array sweeps
3. Optimize array layout for access patterns

**Effort**: High (requires careful profiling and tuning)

### Priority 4: Load Balancing

**Actions**:
1. Use `schedule(dynamic)` instead of `schedule(static)` for unbalanced loops
2. Profile thread utilization to identify load imbalance

**Effort**: Low

## Combined Performance Potential

| Optimization | Effort | Speedup | Cumulative |
|--------------|--------|---------|------------|
| Baseline (current) | - | 1.0x | 1.0x |
| Merge parallel regions | Medium | 1.35x | 1.35x |
| Optimize chain propagation | Medium | 1.15x | 1.55x |
| Cache tiling | High | 1.1x | 1.71x |

**Total potential speedup: ~1.5-1.7x** from parallelization improvements alone.

## Next Steps

1. **Immediate**: Implement parallel region merging (Priority 1)
2. **Short-term**: Generate vectorization reports and optimize loops (Priority 2)
3. **Long-term**: Profile cache behavior and implement tiling (Priority 3)

## Tools for Further Analysis

- **Xcode Instruments**: Time Profiler, System Trace for detailed thread behavior
- **perf** (Linux): Hardware counter analysis for cache misses, branch mispredictions
- **Compiler reports**: `gfortran -fopt-info-vec-all` for vectorization analysis
- **VTune** (Intel): Comprehensive microarchitecture analysis
