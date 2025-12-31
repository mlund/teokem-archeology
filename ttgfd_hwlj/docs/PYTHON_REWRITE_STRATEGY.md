# Python Rewrite Strategy for ttgfd_hwlj

**Date**: 2025-12-31
**Current Version**: Fortran 90, optimized with 14.7% speedup over baseline
**Target**: Python 3.10+ with numpy/scipy

---

## Executive Summary

This document outlines a strategy for rewriting the ttgfd_hwlj polymer density functional theory code from Fortran 90 to Python using numpy/scipy. The Fortran version is highly optimized (SIMD vectorization, OpenMP parallelization) and achieves 12.69s runtime on test cases. A naive Python port would be **10-100x slower**. This strategy focuses on:

1. **Hybrid approach**: Keep performance-critical kernels in compiled code (Fortran/C/Cython)
2. **NumPy vectorization**: Replace Fortran array operations with numpy equivalents
3. **Numba JIT**: Compile performance-critical loops to native code
4. **Incremental migration**: Start with I/O and control logic, gradually replace compute kernels

**Recommended outcome**: Python interface with compiled performance kernels, targeting **2-5x slowdown** compared to optimized Fortran (acceptable for research code with better maintainability).

---

## Current Fortran Code Structure

### Code Organization

```
ttgfd_hwlj_varbl.f90 (2171 lines)
├── Main program (platem)
│   ├── Initialization and I/O
│   ├── Self-consistent field iteration loop
│   └── Convergence checking and output
└── Subroutines (called within SCF loop)
    ├── CDCALC   - Contact density calculation
    ├── AVEC     - Monomer density averaging
    ├── EBLMNEW  - Chain end Boltzmann factors
    └── EBDU     - External potential (Lennard-Jones)
```

### Performance Characteristics (from profiling)

| Component | % Runtime | Lines | Optimization Status |
|-----------|-----------|-------|---------------------|
| Chain propagation | 41% | ~150 | ✅ SIMD + OpenMP |
| Thread barriers | 37% | - | ✅ Reduced via merging |
| EBDU (LJ potential) | 15.5% | ~60 | ✅ Cache-optimized |
| hvec tabulation | 3.1% | ~40 | ✅ SIMD + OpenMP |
| Other subroutines | 3.3% | ~100 | ✅ OpenMP |

### Key Algorithms

1. **Self-Consistent Field (SCF) Iteration**
   - Picard iteration with density mixing (adaptive + fixed mixing)
   - Convergence tolerance: 1e-8
   - Typical iterations: 500-2000

2. **Chain Propagation** (41% of runtime)
   - Forward/backward recursive propagation through density field
   - Gaussian bond convolution with precomputed kernels
   - Triple nested loops: monomers × z-positions × radial positions
   - **Critical for performance**

3. **Lennard-Jones Convolution** (15.5% of runtime)
   - Precomputed hvec lookup table (angular integration)
   - 3D convolution of density with LJ potential
   - Memory-intensive (6GB COMMON blocks)

4. **Excluded Volume Interactions**
   - Hard-sphere potential with Heaviside function
   - Contact density calculations
   - Spherical colloid boundary conditions

---

## Python Architecture Recommendations

### Option 1: Hybrid Python/Fortran (Recommended)

**Architecture**:
```
Python Interface (run.py extended)
├── Simulation setup and I/O (pure Python)
├── SCF convergence control (Python)
└── Performance kernels (compiled)
    ├── Chain propagation (keep Fortran or rewrite in Cython/C++)
    ├── EBDU convolution (keep Fortran or use scipy.ndimage.convolve)
    └── Field calculations (numpy or Fortran)
```

**Pros**:
- **Best performance**: 1.5-3x slowdown vs pure Fortran
- Gradual migration path
- Python interface for experiments
- Keep existing optimizations

**Cons**:
- Mixed-language complexity
- Build system overhead (f2py, Cython compilation)
- Debugging across language boundaries

**Implementation**:
- Use `f2py` to wrap existing Fortran subroutines
- Replace I/O and control flow with Python
- Gradually replace Fortran kernels with numpy/numba

### Option 2: Pure Python with Numba (Good for prototyping)

**Architecture**:
```
Pure Python with numpy/numba
├── Main loop and I/O (Python)
├── Array operations (numpy)
└── Tight loops (numba @jit decorated)
```

**Pros**:
- Single language
- Easy to modify and experiment
- Good IDE/debugging support
- Numba can achieve near-C performance on suitable code

**Cons**:
- **Performance**: 2-10x slower than Fortran even with Numba
- Numba compilation overhead on first run
- Not all Fortran patterns map well to Numba
- Memory usage may be higher

**Implementation**:
- Translate Fortran to Python with numpy arrays
- Add `@numba.jit(nopython=True)` to hot loops
- Use `numba.prange` for parallelization

### Option 3: Pure NumPy/SciPy (Not Recommended)

**Architecture**:
```
Pure numpy/scipy vectorized operations
```

**Pros**:
- Simplest code
- Maximum Pythonic style
- No compilation step

**Cons**:
- **Severe performance degradation**: 20-100x slower
- Memory overhead from temporary arrays
- Some algorithms don't vectorize well (SCF iteration)

**Only viable if**: Performance is not critical and code runs infrequently.

---

## Fortran-to-Python Pattern Mapping

### 1. Array Operations

**Fortran**:
```fortran
! Allocate and initialize
allocatable :: rho(:,:), c(:,:,:)
allocate(rho(nrho, nz))
rho = 0.0d0

! Element-wise operations
rho(i,j) = rho(i,j) + delta * exp(-factor)

! Array slicing
rhoz2 = rho0**2 + z**2
```

**Python (NumPy)**:
```python
# Allocate and initialize
rho = np.zeros((nrho, nz), dtype=np.float64)

# Element-wise operations (vectorized)
rho[i,j] += delta * np.exp(-factor)

# Array operations (fully vectorized)
rhoz2 = rho0**2 + z**2
```

**Performance**: ✅ Equivalent (both use SIMD)

### 2. Nested Loops (Chain Propagation)

**Fortran** (ttgfd_hwlj_varbl.f90:511-521):
```fortran
!$omp parallel do private(iz, jz, irho, ...)
do imon = 1, nmon
  do jz = 1, nz
    do irho = 1, nrho
      !$omp simd reduction(+:phisum)
      do iphi = 1, nphi
        rho2 = rhoz2 - fphi*cos_phi(iphi)
        rsq1 = rho2 + zpcsq
        rsq2 = rho2 + zpc2sq
        valid = merge(1.0d0, 0.0d0, rsq1 >= Rcoll2 .and. rsq2 >= Rcoll2)
        rho = dsqrt(rho2)
        irho_idx = int(rho*rdrho) + 1
        phisum = phisum + valid * cA(irho_idx, jz)
      end do
    end do
  end do
end do
```

**Python Option A (Numba)**:
```python
@numba.jit(nopython=True, parallel=True)
def chain_propagation(nmon, nz, nrho, nphi, cos_phi, cA, ...):
    for imon in numba.prange(nmon):
        for jz in range(nz):
            for irho in range(nrho):
                phisum = 0.0
                for iphi in range(nphi):
                    rho2 = rhoz2 - fphi * cos_phi[iphi]
                    rsq1 = rho2 + zpcsq
                    rsq2 = rho2 + zpc2sq
                    valid = 1.0 if (rsq1 >= Rcoll2 and rsq2 >= Rcoll2) else 0.0
                    rho = np.sqrt(rho2)
                    irho_idx = int(rho * rdrho)
                    phisum += valid * cA[irho_idx, jz]
                # ... use phisum
```

**Performance**: ⚠️ Numba: 1.5-3x slower than Fortran SIMD

**Python Option B (Keep Fortran via f2py)**:
```python
# Interface to Fortran subroutine
import fortran_kernels  # compiled with f2py

result = fortran_kernels.chain_propagation(
    nmon, nz, nrho, nphi, cos_phi, cA, ...
)
```

**Performance**: ✅ Same as Fortran

### 3. OpenMP Parallelization

**Fortran**:
```fortran
!$omp parallel do private(i, j) shared(array)
do i = 1, n
  array(i) = compute(i)
end do
!$omp end parallel do
```

**Python (Numba)**:
```python
@numba.jit(nopython=True, parallel=True)
def compute_parallel(n, array):
    for i in numba.prange(n):  # Parallel range
        array[i] = compute(i)
```

**Python (multiprocessing)**:
```python
from multiprocessing import Pool

def compute_wrapper(i):
    return compute(i)

with Pool() as pool:
    results = pool.map(compute_wrapper, range(n))
```

**Performance**:
- Numba `prange`: ✅ Good (low overhead)
- multiprocessing: ⚠️ High overhead for fine-grained parallelism

### 4. Convolution Operations

**Fortran** (EBDU subroutine):
```fortran
do irho = 1, nrho
  do jz = 1, nz
    do krho = 1, nrho
      do kz = 1, nz
        eb(irho,jz) = eb(irho,jz) + hvec(idist) * rho(krho,kz)
      end do
    end do
  end do
end do
```

**Python (SciPy)**:
```python
from scipy.ndimage import convolve

# If hvec can be expressed as convolution kernel
eb = convolve(rho, hvec_kernel, mode='constant')
```

**Python (NumPy - manual)**:
```python
# For irregular kernels (e.g., LJ potential on cylindrical grid)
for irho in range(nrho):
    for jz in range(nz):
        eb[irho, jz] = np.sum(hvec[distances] * rho)
```

**Performance**:
- `scipy.ndimage.convolve`: ✅ Highly optimized (C backend)
- Manual loops: ❌ 10-50x slower unless Numba JIT

### 5. File I/O

**Fortran**:
```fortran
open(unit=38, file='fcdfil', form='formatted')
write(38, '(2E20.12)') rho_array(i), z_coord(i)
close(38)
```

**Python**:
```python
import numpy as np

# NumPy text format (slower but human-readable)
np.savetxt('fcdfil', np.column_stack([rho_array, z_coord]), fmt='%.12E')

# NumPy binary (faster, not human-readable)
np.save('fcdfil.npy', {'rho': rho_array, 'z': z_coord})

# HDF5 (best for large datasets)
import h5py
with h5py.File('fcdfil.h5', 'w') as f:
    f.create_dataset('rho', data=rho_array)
    f.create_dataset('z', data=z_coord)
```

**Performance**: Binary formats 10-100x faster for large arrays

---

## Migration Strategy: Phased Approach

### Phase 1: Python Interface (1-2 weeks)

**Goal**: Replace run.py and I/O with full Python interface, keep Fortran compute kernels

**Tasks**:
1. Wrap Fortran subroutines with `f2py`:
   ```bash
   f2py -c ttgfd_hwlj_varbl.f90 -m fortran_kernels
   ```

2. Create Python class for simulation:
   ```python
   class PolymerDFT:
       def __init__(self, params):
           self.params = params
           self.rho = None  # Density fields

       def run_scf(self):
           while not converged:
               fortran_kernels.chain_propagation(...)
               fortran_kernels.ebdu(...)
               # Convergence check in Python
   ```

3. Replace file I/O with NumPy:
   - Read/write arrays with `np.load/save`
   - Use HDF5 for large datasets

**Expected Performance**: ✅ ~1.05x slowdown (minimal overhead)

### Phase 2: NumPy Field Calculations (2-3 weeks)

**Goal**: Replace simple subroutines (CDCALC, AVEC, EBLMNEW) with NumPy

**Tasks**:
1. Translate CDCALC (contact density):
   ```python
   def cdcalc(rho, dhs, nrho, nz):
       # Vectorized NumPy operations
       cd = np.zeros_like(rho)
       # ... translate Fortran logic
       return cd
   ```

2. Translate AVEC (averaging):
   ```python
   def avec(c, weights):
       return np.average(c, axis=2, weights=weights)
   ```

3. Keep performance-critical kernels in Fortran initially

**Expected Performance**: ⚠️ 1.1-1.3x slowdown (minor overhead from Python/NumPy)

### Phase 3: Numba JIT for Chain Propagation (3-4 weeks)

**Goal**: Rewrite chain propagation in Numba (most critical kernel)

**Tasks**:
1. Translate chain propagation loop structure to Python
2. Add `@numba.jit(nopython=True, parallel=True)` decorator
3. Optimize array access patterns for Numba
4. Benchmark against Fortran version

**Expected Performance**: ⚠️ 1.5-2.5x slowdown compared to optimized Fortran

**Fallback**: If Numba performance is inadequate, rewrite in Cython:
```cython
# chain_prop.pyx
cimport numpy as np
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
def chain_propagation(double[:,:] rho, ...):
    cdef int i, j, k
    cdef double phisum
    # C-level performance with Python interface
```

### Phase 4: EBDU Convolution (2-3 weeks)

**Goal**: Replace EBDU (Lennard-Jones potential) with SciPy or custom convolution

**Options**:
1. **SciPy convolution** (if kernel is regular):
   ```python
   from scipy.ndimage import convolve
   eb = convolve(rho, hvec_kernel, mode='constant')
   ```

2. **Custom Numba kernel** (for irregular cylindrical geometry):
   ```python
   @numba.jit(nopython=True, parallel=True)
   def ebdu_custom(rho, hvec, distances, ...):
       # Optimized distance-based lookup
   ```

3. **Keep Fortran** (if performance is critical)

**Expected Performance**:
- SciPy: ✅ 1.1-1.5x slowdown
- Numba: ⚠️ 1.5-2.5x slowdown
- Fortran (f2py): ✅ 1.0x (same)

### Phase 5: Optimization and Validation (2-4 weeks)

**Goal**: Optimize Python version, comprehensive validation

**Tasks**:
1. Profile Python code with `cProfile` and `line_profiler`
2. Identify bottlenecks and optimize
3. Memory profiling with `memory_profiler`
4. Comprehensive regression testing against Fortran
5. Document performance trade-offs

**Target**: 2-5x slowdown vs optimized Fortran is acceptable for research code

---

## Performance Expectations

### Realistic Timings (test case: 12.69s in Fortran)

| Implementation | Expected Time | Slowdown | Effort |
|----------------|---------------|----------|--------|
| Current Fortran (optimized) | 12.69s | 1.0x | ✅ Done |
| Phase 1 (Python + f2py) | 13-14s | 1.05-1.1x | Low |
| Phase 2 (+ NumPy fields) | 14-16s | 1.1-1.3x | Low |
| Phase 3 (+ Numba chain prop) | 20-30s | 1.5-2.5x | Medium |
| Phase 4 (+ Numba/SciPy EBDU) | 25-40s | 2.0-3.0x | Medium |
| Pure NumPy (no Numba/f2py) | 100-500s | 10-40x | Low |

### Memory Usage

**Fortran**: ~6GB for large systems (COMMON blocks)
**Python**:
- NumPy arrays: Same as Fortran (efficient)
- Overhead: +10-20% from Python objects
- **Risk**: Temporary arrays in vectorized operations can double memory usage

**Mitigation**:
- Use in-place operations: `+=`, `*=`
- Pre-allocate arrays, avoid creating temporaries
- Use `np.empty()` instead of `np.zeros()` where safe

---

## Recommended Technology Stack

### Core Scientific Stack
```
python >= 3.10
numpy >= 1.24
scipy >= 1.10
numba >= 0.57    # JIT compilation for performance
h5py >= 3.8      # Efficient I/O
```

### Optional Performance Enhancements
```
cython >= 3.0    # If Numba insufficient
mpi4py           # For distributed memory parallelism (large systems)
pybind11         # C++ alternative to f2py
```

### Development Tools
```
pytest           # Testing
pytest-benchmark # Performance regression testing
line_profiler    # Line-by-line profiling
memory_profiler  # Memory usage tracking
black            # Code formatting
mypy             # Static type checking
```

---

## Example: Minimal Python Interface (Phase 1)

### Step 1: Compile Fortran kernels with f2py

```bash
# Create f2py signature file
f2py ttgfd_hwlj_varbl.f90 -m fortran_kernels -h fortran_kernels.pyf

# Edit .pyf file to expose only needed subroutines
# Compile
f2py -c fortran_kernels.pyf ttgfd_hwlj_varbl.f90
```

### Step 2: Python wrapper class

```python
# polymer_dft.py
import numpy as np
from dataclasses import dataclass
from pathlib import Path
import fortran_kernels  # Compiled Fortran

@dataclass
class SimulationParams:
    """Simulation parameters"""
    bdm: float = 0.01           # Monomer bulk density
    nmon: int = 151             # Polymer chain length
    dz: float = 0.25            # Grid spacing
    drho: float = 0.25
    Rcoll: float = 10.0         # Colloid radius
    epslj: float = 1.0          # LJ energy
    ioimaxm: int = 5000         # Max iterations
    dmm: float = 0.9            # Mixing parameter
    # ... other parameters

class PolymerDFT:
    """Python interface to Fortran polymer DFT code"""

    def __init__(self, params: SimulationParams):
        self.params = params
        self.iteration = 0
        self.converged = False

    def initialize_fields(self):
        """Initialize density fields"""
        nrho = int(self.params.Rcyl / self.params.drho) + 1
        nz = int(2 * self.params.zc1 / self.params.dz) + 1

        self.rho_monomer = np.full((nrho, nz), self.params.bdm)
        self.rho_solvent = np.full((nrho, nz), self.params.bds)

    def run_scf_iteration(self):
        """Single SCF iteration using Fortran kernels"""
        # Call Fortran subroutines via f2py
        rho_new = fortran_kernels.chain_propagation(
            self.rho_monomer,
            self.params.nmon,
            # ... pass parameters as arrays
        )

        # Convergence check in Python
        ddmax = np.max(np.abs(rho_new - self.rho_monomer))

        # Mixing
        self.rho_monomer = (self.params.dmm * self.rho_monomer +
                           (1 - self.params.dmm) * rho_new)

        return ddmax

    def run(self):
        """Run full SCF calculation"""
        self.initialize_fields()

        for self.iteration in range(self.params.ioimaxm):
            ddmax = self.run_scf_iteration()

            if ddmax < 1e-8:
                self.converged = True
                break

        return self.converged

    def save_results(self, filename: str):
        """Save results to HDF5"""
        import h5py
        with h5py.File(filename, 'w') as f:
            f.create_dataset('rho_monomer', data=self.rho_monomer)
            f.create_dataset('rho_solvent', data=self.rho_solvent)
            f.attrs['converged'] = self.converged
            f.attrs['iterations'] = self.iteration

# Usage
if __name__ == '__main__':
    params = SimulationParams(
        bdm=0.01,
        epslj=0.5,
        nmon=151
    )

    sim = PolymerDFT(params)
    sim.run()
    sim.save_results('output.h5')
```

---

## Testing Strategy

### Regression Tests Against Fortran

```python
# test_python_vs_fortran.py
import numpy as np
import pytest
from polymer_dft import PolymerDFT, SimulationParams

def test_phase1_matches_fortran():
    """Python interface should match Fortran output exactly"""
    params = SimulationParams(bdm=0.1, epslj=0.5, nmon=151)

    # Run Python version (Phase 1: calls Fortran via f2py)
    sim_py = PolymerDFT(params)
    sim_py.run()

    # Load Fortran reference output
    rho_fortran = np.loadtxt('reference_output_fortran.dat')

    # Should be identical (same numerical precision)
    np.testing.assert_allclose(
        sim_py.rho_monomer.ravel(),
        rho_fortran,
        rtol=1e-12
    )

def test_phase3_numba_matches_fortran():
    """Numba version should match within numerical tolerance"""
    # After replacing kernels with Numba
    # Allow small differences due to floating-point order
    np.testing.assert_allclose(
        result_numba,
        result_fortran,
        rtol=1e-10,
        atol=1e-12
    )

@pytest.mark.benchmark
def test_performance_regression():
    """Ensure Python version stays within 5x of Fortran"""
    import time

    params = SimulationParams(bdm=0.1, epslj=0.5)

    start = time.time()
    sim = PolymerDFT(params)
    sim.run()
    python_time = time.time() - start

    # Fortran baseline: 12.69s
    fortran_time = 12.69

    slowdown = python_time / fortran_time
    assert slowdown < 5.0, f"Too slow: {slowdown:.1f}x vs Fortran"
```

---

## Trade-offs and Recommendations

### When to Use Python Version

✅ **Good use cases**:
- **Research and exploration**: Testing new theories, parameters
- **Visualization and analysis**: Integrate with matplotlib, pandas
- **Workflow integration**: Combine with other Python tools
- **Teaching**: More accessible than Fortran for students
- **Prototyping**: Faster development cycle

❌ **Avoid Python for**:
- **Production runs**: Large-scale simulations (use Fortran)
- **HPC clusters**: Fortran better for MPI + GPU
- **Time-critical applications**: Real-time processing

### Final Recommendation: Hybrid Approach

**Best of both worlds**:
1. **Keep Fortran for performance**: Compile as shared library
2. **Python for interface**: I/O, visualization, experimentation
3. **Incremental migration**: Replace Fortran kernels only when Python equivalent is validated

**Implementation timeline**:
- **Phase 1** (Python + f2py): 1-2 weeks → **Immediate value**
- **Phase 2-3** (NumPy + Numba): 2-3 months → **Most code in Python**
- **Phase 4** (Full Python): 3-6 months → **Complete rewrite**

**Suggested approach**:
- Start with Phase 1 to get Python interface working
- Evaluate if Phases 2-4 are worth the performance trade-off
- Consider keeping performance-critical kernels (chain propagation, EBDU) in Fortran long-term

---

## Additional Resources

### f2py Documentation
- https://numpy.org/doc/stable/f2py/
- Wrapping Fortran for Python: https://www.numfys.net/howto/F2PY/

### Numba Best Practices
- https://numba.readthedocs.io/en/stable/user/performance-tips.html
- Parallelization: https://numba.readthedocs.io/en/stable/user/parallel.html

### SciPy Convolution
- `scipy.ndimage.convolve`: https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.convolve.html
- FFT convolution: `scipy.signal.fftconvolve`

### Cython (if Numba insufficient)
- https://cython.readthedocs.io/
- Faster than Numba for complex control flow, comparable for tight loops

### Performance Profiling
- `cProfile`: Built-in Python profiler
- `line_profiler`: Line-by-line timing
- `py-spy`: Sampling profiler with flame graphs
- `scalene`: Modern profiler with memory tracking

---

## Conclusion

A full Python rewrite of ttgfd_hwlj is feasible but requires careful attention to performance. The recommended strategy is:

1. **Start small**: Phase 1 (Python + f2py) provides immediate benefits with minimal risk
2. **Validate thoroughly**: Compare every change against Fortran reference
3. **Accept trade-offs**: 2-5x slowdown is reasonable for better maintainability
4. **Stay pragmatic**: Keep Fortran kernels if Python replacements are too slow

The optimized Fortran version (12.69s) is excellent for production. The Python version excels at flexibility, integration, and rapid experimentation. Use the right tool for the job.
