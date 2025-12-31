# Fortran 90 Modernization Opportunities for ttgfd_hwlj

**Date**: 2025-12-31
**Current Status**: Fortran 90 with legacy F77 patterns
**Target**: Modern Fortran 2008/2018 best practices

---

## Executive Summary

The code is currently Fortran 90 with significant Fortran 77 legacy patterns. The most impactful modernizations would be:

1. **Replace COMMON blocks with modules** (HIGH impact, MEDIUM effort)
2. **Use ISO_FORTRAN_ENV types** (LOW impact, LOW effort)
3. **Modern file I/O with newunit** (LOW impact, LOW effort)
4. **Better code organization** (MEDIUM impact, HIGH effort)

**Recommended approach**: Incremental modernization, starting with low-risk improvements.

---

## Current Legacy Patterns

### 1. COMMON Blocks (Fortran 77 style)

**Location**: `t2.inc.f90` lines 17, 37, 40

**Current code**:
```fortran
! Large arrays in COMMON block (subroutines need access via COMMON)
DOUBLE PRECISION :: fdmon(0:maxrho, 0:maxel), ebelam(0:maxrho, 0:maxel)
DOUBLE PRECISION :: convp(0:maxrho, 0:maxel), hvec(0:maxrho, 0:maxrho, 0:maxel)
...
COMMON/VECT/fdmon, ebelam, convp, hvec, fem, ehbclam, cdmonm, ae1, ae2, edu, cos_phi, cos_pphi

DOUBLE PRECISION :: dz, scalem, emscale, rdz, Yfact, ...
COMMON/VAR/dz, scalem, emscale, rdz, Yfact, ...

INTEGER :: istart, istp1, islut, ism, nfack, imitt, nmon
COMMON/HELTAL/istart, istp1, islut, ism, nfack, ...
```

**Problems**:
- No encapsulation or access control
- Global mutable state makes code hard to reason about
- Error-prone: order in COMMON must match across all files
- Cannot use allocatable arrays in COMMON (fixed size only)
- Not thread-safe (though we use OpenMP, COMMON is shared)

**Modern alternative**: Module with derived types or module variables
```fortran
module polymer_dft_data
  use iso_fortran_env, only: real64, int32
  implicit none

  ! Grid parameters
  type :: grid_params_t
    real(real64) :: dz, drho, dphi
    real(real64) :: rdz, rdrho
    integer(int32) :: mxrho, nphi, nmon
  end type

  type(grid_params_t), save :: grid

  ! Field arrays (allocatable for flexibility)
  real(real64), allocatable, save :: fdmon(:,:), ebelam(:,:)
  real(real64), allocatable, save :: hvec(:,:,:)
  ...

contains
  subroutine initialize_arrays(nrho, nz)
    integer, intent(in) :: nrho, nz
    allocate(fdmon(0:nrho, 0:nz))
    allocate(ebelam(0:nrho, 0:nz))
    ...
  end subroutine

  subroutine cleanup_arrays()
    deallocate(fdmon, ebelam, hvec, ...)
  end subroutine
end module
```

**Benefits**:
- Explicit dependencies via `use polymer_dft_data`
- Allocatable arrays adapt to problem size
- Better documentation and organization
- Type safety with derived types
- Easier to parallelize (can use threadprivate if needed)

**Effort**: MEDIUM (need to refactor all subroutines to use module)
**Risk**: LOW (mechanical transformation, easy to test)

---

### 2. Include Files Instead of Modules

**Location**: `ttgfd_hwlj_varbl.f90` line 11

**Current code**:
```fortran
program platem
  implicit none
  include 't2.inc.f90'
```

**Problems**:
- Include is textual substitution (C preprocessor style)
- No namespace control
- Cannot selectively import
- Include files can't contain their own `implicit none`
- Makes dependencies unclear

**Modern alternative**:
```fortran
program platem
  use iso_fortran_env, only: real64, int32
  use polymer_dft_data
  use polymer_dft_constants, only: PI, TWOPI
  implicit none
```

**Benefits**:
- Explicit dependencies
- Can use `only:` for selective import
- Better IDE support (completion, go-to-definition)
- Module can enforce `implicit none`

**Effort**: LOW (convert t2.inc.f90 to module)
**Risk**: LOW

---

### 3. Old-Style Type Declarations

**Location**: Throughout code

**Current code**:
```fortran
DOUBLE PRECISION :: variable
INTEGER :: index
```

**Modern alternative**:
```fortran
use iso_fortran_env, only: real64, int32, int64
real(real64) :: variable
integer(int32) :: index
```

**Benefits**:
- Portable precision (8-byte real guaranteed)
- Explicit about integer sizes (important for large arrays)
- Standard across platforms
- Better interoperability with C

**Effort**: LOW (search and replace)
**Risk**: VERY LOW

---

### 4. Fixed-Size Arrays

**Location**: `t2.inc.f90` lines 2-4, 7-15

**Current code**:
```fortran
INTEGER, PARAMETER :: maxel = 1001
INTEGER, PARAMETER :: maxrho = 321
DOUBLE PRECISION :: fdmon(0:maxrho, 0:maxel)
```

**Problems**:
- Wastes memory for small problems
- Hard limit on problem size
- Cannot grow dynamically

**Modern alternative**:
```fortran
real(real64), allocatable :: fdmon(:,:)
! Later...
allocate(fdmon(0:nrho, 0:nz))
```

**Benefits**:
- Adapt to problem size
- Can run larger systems without recompiling
- Better memory usage

**Effort**: MEDIUM (need to add allocation/deallocation)
**Risk**: LOW (but need careful testing of bounds)

**Note**: Currently limited by ~6GB COMMON block size. Modules don't have this limit.

---

### 5. Hardcoded File Unit Numbers

**Location**: `ttgfd_hwlj_varbl.f90` lines 79-84

**Current code**:
```fortran
ifc = 38  ! Output: concentration profiles
ins = 49  ! Input: simulation parameters
iep = 50  ! Input: Lennard-Jones epsilon parameter

open (ifc, file='fcdfil', form='formatted')
open (ins, file='input.tsph', form='formatted')
```

**Problems**:
- Magic numbers (why 38, 49, 50?)
- Risk of conflicts with other libraries
- Manual tracking of units

**Modern alternative**:
```fortran
integer :: ifc, ins, iep

open (newunit=ifc, file='fcdfil', form='formatted', status='replace')
open (newunit=ins, file='input.tsph', form='formatted', status='old')
```

**Benefits**:
- Automatic unit allocation by runtime
- No conflicts
- Cleaner code

**Effort**: VERY LOW
**Risk**: VERY LOW

---

### 6. Mixed Case and Style

**Location**: Throughout

**Current code**:
```fortran
DOUBLE PRECISION :: fdmon(0:maxrho, 0:maxel)  ! Uppercase
double precision :: variable                  ! Lowercase
INTEGER, PARAMETER :: maxel = 1001            ! Uppercase
```

**Modern style**: Consistent lowercase (community standard)
```fortran
real(real64) :: fdmon(0:maxrho, 0:maxel)
integer, parameter :: maxel = 1001
```

**Effort**: VERY LOW (fprettify can do this)
**Risk**: VERY LOW

---

### 7. Subroutines Without Modules

**Location**: Lines 1204-1601 in `ttgfd_hwlj_varbl.f90`

**Current structure**:
```fortran
program platem
  ! Main program code
  ! ... 1200 lines ...
end program

subroutine CDFACT
  ! ...
end subroutine

subroutine CDCALC
  ! ...
end subroutine
```

**Problems**:
- Subroutines are global scope
- No encapsulation
- Long single file (1600 lines)
- Harder to maintain

**Modern alternative**:
```fortran
module polymer_dft_subroutines
  use polymer_dft_data
  implicit none
  private
  public :: cdfact, cdcalc, avec, eblmnew, ebdu

contains
  subroutine cdfact()
    ! Can access module variables directly
  end subroutine

  subroutine cdcalc()
    ! ...
  end subroutine
end module

program platem
  use polymer_dft_data
  use polymer_dft_subroutines
  implicit none
  ! Main program
end program
```

**Benefits**:
- Better organization
- Can split into multiple files
- Explicit dependencies
- Private helper functions possible

**Effort**: HIGH (major refactoring)
**Risk**: MEDIUM (need careful testing)

---

## Recommended Modernization Plan

### Phase 1: Low-Risk Quick Wins (1-2 days)

**Goal**: Modernize without changing logic

1. **Use ISO_FORTRAN_ENV types**
   ```fortran
   use iso_fortran_env, only: real64, int32
   real(real64) :: variable
   integer(int32) :: index
   ```
   - Effort: LOW (search/replace)
   - Risk: VERY LOW
   - Benefit: Portable, explicit precision

2. **Use newunit for file I/O**
   ```fortran
   open(newunit=ifc, file='fcdfil', ...)
   ```
   - Effort: VERY LOW
   - Risk: VERY LOW
   - Benefit: Cleaner, safer

3. **Consistent lowercase style**
   ```bash
   fprettify --case 1 1 1 1 *.f90
   ```
   - Effort: VERY LOW (automated)
   - Risk: VERY LOW
   - Benefit: Modern, readable

**Expected outcome**: Code looks more modern, no functional changes

---

### Phase 2: Convert COMMON Blocks to Module (1 week)

**Goal**: Eliminate COMMON blocks while maintaining compatibility

**Steps**:

1. Create `polymer_dft_data.f90`:
   ```fortran
   module polymer_dft_data
     use iso_fortran_env, only: real64, int32
     implicit none

     ! Constants
     integer, parameter :: maxel = 1001
     integer, parameter :: maxrho = 321

     ! Module variables (was COMMON/VECT/)
     real(real64), save :: fdmon(0:maxrho, 0:maxel)
     real(real64), save :: ebelam(0:maxrho, 0:maxel)
     ! ... rest of arrays

     ! Parameters (was COMMON/VAR/)
     real(real64), save :: dz, drho, dphi
     ! ... rest of variables

     ! Grid dimensions (was COMMON/HELTAL/)
     integer(int32), save :: mxrho, nphi, nmon
     ! ... rest of integers
   end module
   ```

2. Replace `include 't2.inc.f90'` with `use polymer_dft_data`

3. Remove `t2.inc.f90` include from all subroutines, add `use polymer_dft_data`

4. Test thoroughly with regression suite

**Effort**: MEDIUM
**Risk**: LOW (mechanical transformation)
**Benefit**: Foundation for further improvements

---

### Phase 3: Use Allocatable Arrays (2-3 weeks)

**Goal**: Dynamic sizing based on input

**Steps**:

1. Make arrays allocatable in module:
   ```fortran
   real(real64), allocatable, save :: fdmon(:,:)
   ```

2. Add initialization routine:
   ```fortran
   subroutine initialize_arrays(nrho, nz)
     allocate(fdmon(0:nrho, 0:nz))
     allocate(ebelam(0:nrho, 0:nz))
     ! ...
   end subroutine
   ```

3. Call from main program after reading grid parameters

4. Add cleanup routine

**Effort**: MEDIUM
**Risk**: MEDIUM (array bounds need careful checking)
**Benefit**:
- No artificial size limits
- Better memory usage
- Can run huge systems

---

### Phase 4: Organize into Modules (3-4 weeks) - OPTIONAL

**Goal**: Better code organization

**Proposed structure**:
```
polymer_dft_constants.f90  - Physical/math constants
polymer_dft_types.f90      - Derived types
polymer_dft_data.f90       - Module variables
polymer_dft_io.f90         - File I/O routines
polymer_dft_fields.f90     - Field calculation subroutines
polymer_dft_chain.f90      - Chain propagation
polymer_dft.f90            - Main program
```

**Effort**: HIGH
**Risk**: MEDIUM
**Benefit**: Much easier to maintain and extend

---

## Non-Recommended Changes

### 1. ❌ Object-Oriented Refactoring

**Why skip**:
- Fortran OOP is verbose
- No clear benefit for this numerical code
- Would slow down development
- OOP doesn't improve performance

### 2. ❌ Submodules

**Why skip**:
- Compiler support still spotty (especially gfortran)
- Adds complexity without clear benefit
- Modules are sufficient

### 3. ❌ Coarrays for Parallelism

**Why skip**:
- OpenMP works well already
- Coarrays require special runtime
- Not worth the porting effort

---

## Compatibility Considerations

### Fortran Standard Support

| Feature | Standard | gfortran | ifort | nvfortran |
|---------|----------|----------|-------|-----------|
| Modules | F90 | ✅ | ✅ | ✅ |
| ISO_FORTRAN_ENV | F2003 | ✅ | ✅ | ✅ |
| Allocatable in modules | F90 | ✅ | ✅ | ✅ |
| newunit= | F2008 | ✅ 4.6+ | ✅ | ✅ |
| Submodules | F2008 | ⚠️ 6.0+ | ✅ | ⚠️ |

**Recommendation**: Target Fortran 2008 (widely supported)

---

## Testing Strategy

For each modernization phase:

1. **Unit testing**: Test individual subroutines
2. **Regression testing**: Run all 5 existing regression tests
3. **Numerical validation**: Ensure bit-for-bit identical results
4. **Performance testing**: Verify no slowdown

**Key validation**:
```bash
# Before changes
./ttgfd_hwlj < input.tsph > output_before.txt

# After changes
./ttgfd_hwlj < input.tsph > output_after.txt

# Compare
diff output_before.txt output_after.txt  # Should be identical
```

---

## Estimated Effort Summary

| Phase | Effort | Risk | Benefit | Priority |
|-------|--------|------|---------|----------|
| Phase 1: Quick wins | 1-2 days | Very Low | Low | HIGH |
| Phase 2: Module conversion | 1 week | Low | High | HIGH |
| Phase 3: Allocatable arrays | 2-3 weeks | Medium | High | MEDIUM |
| Phase 4: Full reorganization | 3-4 weeks | Medium | Medium | LOW |

**Recommended**: Do Phases 1 and 2 (total ~1.5 weeks), defer Phase 3-4 unless needed.

---

## Example: Phase 1 Quick Wins

Here's what a minimally modernized version would look like:

**Before** (`ttgfd_hwlj_varbl.f90`):
```fortran
program platem
  implicit none
  include 't2.inc.f90'

  INTEGER, PARAMETER :: MAXMON = 151
  double precision, allocatable :: c(:, :, :)
  integer :: ifc, ins

  ifc = 38
  open (ifc, file='fcdfil', form='formatted')
```

**After** (Phase 1 only):
```fortran
program platem
  use iso_fortran_env, only: real64, int32
  implicit none
  include 't2.inc.f90'

  integer(int32), parameter :: maxmon = 151
  real(real64), allocatable :: c(:, :, :)
  integer(int32) :: ifc, ins

  open (newunit=ifc, file='fcdfil', form='formatted', status='replace')
```

**Changes**:
- Uses `iso_fortran_env` for portable types
- Lowercase keywords (modern style)
- `newunit=` for automatic file unit allocation
- Explicit status for file open

**Testing**: Should produce identical numerical results

---

## Conclusion

The code is well-written Fortran 90 but uses Fortran 77 patterns (COMMON blocks, includes).

**Immediate recommendation**:
1. Start with Phase 1 (quick wins) - minimal risk, modernizes appearance
2. Consider Phase 2 (module conversion) - unlocks future improvements

**Long-term**: Phases 3-4 only if planning significant new features or want to remove size limits.

The current code works well and is optimized. Modernization is about **maintainability and future-proofing**, not performance.
