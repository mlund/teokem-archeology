! ============================================================================
! MODULE: polymer_dft_data
! ============================================================================
! Contains all shared data for the polymer DFT calculation.
! This module replaces the old COMMON blocks and include file approach.
!
! Migrated from: t2.inc.f90
! Date: 2025-12-31
! Phase 3: Converted to allocatable arrays - 2026-01-01
! ============================================================================
module polymer_dft_data
  use iso_fortran_env, only: real64, int32
  implicit none

  ! ========================================================================
  ! Derived types for structured data organization
  ! ========================================================================

  ! Input parameters read directly from input files
  ! This encapsulates all user-specified parameters from input.tsph and epfil
  type :: input_params_t
    ! Bulk densities (from input.tsph)
    real(real64) :: bdm         ! Monomer bulk density
    real(real64) :: bdtot       ! Total bulk density

    ! Polymer parameters (from input.tsph)
    integer(int32) :: nmon      ! Number of monomers per polymer chain

    ! Grid spacing parameters (from input.tsph)
    real(real64) :: dz          ! Grid spacing in z direction
    real(real64) :: drho        ! Grid spacing in radial direction
    real(real64) :: dphi        ! Angular grid spacing (in units of pi)
    real(real64) :: dpphi       ! Angular grid spacing for potential calculation

    ! Geometry parameters (from input.tsph)
    real(real64) :: Rcoll       ! Colloid radius
    real(real64) :: zc1         ! Position of first colloid center
    real(real64) :: collsep     ! Separation between colloid centers
    real(real64) :: Rcyl        ! Cylinder radius (system boundary)

    ! Algorithm parameters (from input.tsph)
    integer(int32) :: ioimaxm   ! Maximum number of iterations
    real(real64) :: dmm         ! Density mixing parameter for monomer
    real(real64) :: dms         ! Density mixing parameter for solvent
    integer(int32) :: kread     ! Restart flag (0=fresh start, 1=restart from file)

    ! Molecular parameters (from input.tsph)
    real(real64) :: bl          ! Bond length
    real(real64) :: dhs         ! Hard sphere diameter (monomer)

    ! Lennard-Jones parameters (from epfil)
    real(real64) :: epslj       ! Lennard-Jones epsilon parameter
  end type input_params_t

  ! ========================================================================
  ! Module variables - structured input parameters
  ! ========================================================================
  type(input_params_t) :: input  ! User input parameters

  ! ========================================================================
  ! Module-level aliases for commonly-used input parameters
  ! (For backward compatibility with subroutines)
  ! These are set from input% after reading
  ! ========================================================================
  real(real64) :: dz, drho, dphi     ! Grid spacings
  real(real64) :: dhs, bl            ! Molecular parameters
  real(real64) :: bdm                ! Monomer bulk density
  real(real64) :: zc1, Rcoll         ! Geometry

  ! ========================================================================
  ! Derived parameters (computed from input, shared across subroutines)
  ! ========================================================================
  ! Grid spacing reciprocals
  real(real64) :: rdz, rdrho, rdphi

  ! Geometric derived parameters
  real(real64) :: dhs2, dhs3, rdhs3  ! Hard sphere diameter: squared, cubed, 1/cubed
  real(real64) :: bl2                ! Bond length squared
  real(real64) :: dzrfp              ! dz/(4*pi)
  real(real64) :: twopidz            ! 2*pi*dz

  ! ========================================================================
  ! Array dimensions - now allocatable (no fixed limits)
  ! ========================================================================
  integer, parameter :: maxphi = 5000  ! Maximum phi grid points for lookup tables

  ! ========================================================================
  ! Large arrays - now allocatable (formerly fixed-size in COMMON/VECT/)
  ! ========================================================================
  real(real64), allocatable :: fdmon(:, :)
  real(real64), allocatable :: ebelam(:, :)
  real(real64), allocatable :: convp(:, :)
  real(real64), allocatable :: hvec(:, :, :)
  real(real64), allocatable :: fem(:, :)
  real(real64), allocatable :: ehbclam(:, :)
  real(real64), allocatable :: cdmonm(:, :)
  real(real64), allocatable :: ae1(:, :)
  real(real64), allocatable :: ae2(:, :)
  real(real64), allocatable :: edu(:, :)

  ! Cosine lookup tables for performance optimization
  real(real64) :: cos_phi(maxphi)
  real(real64) :: cos_pphi(maxphi)

  ! ========================================================================
  ! Other module-level variables (computed values, not from input)
  ! ========================================================================
  real(real64) :: scalem, emscale
  real(real64) :: Yfact, rnmon, rrnmon
  real(real64) :: chemps
  real(real64) :: Rcoll2, behbclam, bebelam
  real(real64) :: zc2
  real(real64) :: cdnorm

  ! ========================================================================
  ! Integer grid dimensions (formerly in COMMON/HELTAL/)
  ! ========================================================================
  integer(int32) :: istart, istp1, islut, ism, nfack, imitt, nmon
  integer(int32) :: istp1s, isluts, mxrho, ksm, nphi
  integer(int32) :: ibl, kbl

  ! ========================================================================
  ! Physical and mathematical constants (compile-time parameters)
  ! ========================================================================

  ! Mathematical constants
  ! Note: Precision limited to double precision (15-17 significant digits)
  real(real64), parameter :: PI = 3.14159265358979d0
  real(real64), parameter :: TWOPI = 2.d0*PI
  real(real64), parameter :: FOURPI = 4.d0*PI

  ! Carnahan-Starling equation of state parameters for hard spheres
  real(real64), parameter :: A1_CS = 1.d0
  real(real64), parameter :: A2_CS = 2.45696d0
  real(real64), parameter :: B1_CS = 1.d0
  real(real64), parameter :: B2_CS = 4.10386d0

  ! Convergence tolerance
  real(real64), parameter :: CONV_TOL = 0.00001d0

  ! Derived constants (computed from fundamental constants)
  real(real64), parameter :: PIS = PI/6.d0
  real(real64), parameter :: C1 = -1.d0
  real(real64), parameter :: C2 = -3.75503d0
  real(real64), parameter :: AA1 = 2.d0*C1 - 2.d0*A1_CS - 4.d0
  real(real64), parameter :: AA2 = 2.d0*C2 - 2.d0*A2_CS - 4.d0
  real(real64), parameter :: BB1 = 3.d0 - B1_CS + A1_CS - 3.d0*C1
  real(real64), parameter :: BB2 = 3.d0 - B2_CS + A2_CS - 3.d0*C2
  real(real64), parameter :: Y = (9.82605d0 - 9.d0*PI*0.25d0)/(9.d0*PI*0.25d0 - 4.d0*PI/3.d0)

contains

  ! ==========================================================================
  ! SUBROUTINE: allocate_arrays
  ! ==========================================================================
  ! Allocates all dynamic arrays based on grid dimensions
  ! Must be called after mxrho, imitt, and nfack are calculated from input
  ! ==========================================================================
  subroutine allocate_arrays(nrho, nz, nz_hvec)
    integer(int32), intent(in) :: nrho     ! Maximum rho grid points (mxrho + kbl)
    integer(int32), intent(in) :: nz       ! Maximum z grid points (imitt + ibl)
    integer(int32), intent(in) :: nz_hvec  ! Maximum z for hvec (nfack - 1)

    ! Allocate 2D arrays with 0-based indexing
    allocate(fdmon(0:nrho, 0:nz))
    allocate(ebelam(0:nrho, 0:nz))
    allocate(convp(0:nrho, 0:nz))
    allocate(fem(0:nrho, 0:nz))
    allocate(ehbclam(0:nrho, 0:nz))
    allocate(cdmonm(0:nrho, 0:nz))
    allocate(ae1(0:nrho, 0:nz))
    allocate(ae2(0:nrho, 0:nz))
    allocate(edu(0:nrho, 0:nz))

    ! Allocate 3D array (hvec has different z-dimension for LJ potential table)
    allocate(hvec(0:nrho, 0:nrho, 0:nz_hvec))

    ! Initialize arrays to zero
    fdmon = 0.d0
    ebelam = 0.d0
    convp = 0.d0
    hvec = 0.d0
    fem = 0.d0
    ehbclam = 0.d0
    cdmonm = 0.d0
    ae1 = 0.d0
    ae2 = 0.d0
    edu = 0.d0

  end subroutine allocate_arrays

  ! ==========================================================================
  ! SUBROUTINE: deallocate_arrays
  ! ==========================================================================
  ! Deallocates all dynamic arrays
  ! Should be called before program termination for clean memory management
  ! ==========================================================================
  subroutine deallocate_arrays()
    if (allocated(fdmon)) deallocate(fdmon)
    if (allocated(ebelam)) deallocate(ebelam)
    if (allocated(convp)) deallocate(convp)
    if (allocated(hvec)) deallocate(hvec)
    if (allocated(fem)) deallocate(fem)
    if (allocated(ehbclam)) deallocate(ehbclam)
    if (allocated(cdmonm)) deallocate(cdmonm)
    if (allocated(ae1)) deallocate(ae1)
    if (allocated(ae2)) deallocate(ae2)
    if (allocated(edu)) deallocate(edu)
  end subroutine deallocate_arrays

end module polymer_dft_data
