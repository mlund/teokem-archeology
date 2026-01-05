! MIT License
!
! Copyright (c) 2025 Jan Forsman and Mikael Lund
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
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

  ! Grid dimensions and discretization parameters
  ! Computed from input_params_t during initialization
  type :: grid_params_t
    ! Primary grid loop bounds
    integer(int32) :: istart    ! Starting index for z-grid (always 0)
    integer(int32) :: istp1     ! First interior z-index (always 1)
    integer(int32) :: islut     ! Last z-index (= nfack)
    integer(int32) :: imitt     ! Midpoint z-index (= nfack/2)
    integer(int32) :: nfack     ! Total z-grid points for full system

    ! Secondary grid bounds for solvent
    integer(int32) :: istp1s    ! Solvent start index (= istp1)
    integer(int32) :: isluts    ! Solvent end index (= islut)

    ! Maximum grid extents
    integer(int32) :: mxrho     ! Maximum radial grid points

    ! Discretized molecular dimensions (in grid units)
    integer(int32) :: ism       ! Hard sphere diameter in z-direction grid units
    integer(int32) :: ksm       ! Hard sphere diameter in rho-direction grid units
    integer(int32) :: ibl       ! Bond length in z-direction grid units
    integer(int32) :: kbl       ! Bond length in rho-direction grid units

    ! Angular discretization
    integer(int32) :: nphi      ! Number of angular grid points

    ! Polymer chain parameters
    integer(int32) :: nmon      ! Number of monomers (copy from input)
  end type grid_params_t

  ! Bulk thermodynamic properties from Carnahan-Starling equation of state
  type :: bulk_properties_t
    real(real64) :: bdpol      ! Polymer bulk density (bdm/nmon)
    real(real64) :: bdt        ! Monomer packing density (bdm*dhs^3)
    real(real64) :: Pb         ! Bulk pressure
    real(real64) :: bFex       ! Bulk excess free energy
    real(real64) :: chempp     ! Polymer chemical potential
    real(real64) :: emtrams    ! End-segment Boltzmann factor term
    real(real64) :: cmtrams    ! Middle-segment Boltzmann factor term
  end type bulk_properties_t

  ! Derived parameters computed from input_params_t and grid_params_t
  type :: computed_params_t
    ! Grid spacing reciprocals
    real(real64) :: rdz, rdrho, rdphi

    ! Molecular geometry parameters
    real(real64) :: dhs2, dhs3, rdhs3  ! dhs^2, dhs^3, 1/dhs^3
    real(real64) :: bl2                ! bl^2

    ! Grid integration helpers
    real(real64) :: dzrfp       ! dz/(4*pi)
    real(real64) :: twopidz     ! 2*pi*dz

    ! Polymer chain parameters
    real(real64) :: rnmon       ! Real-valued nmon
    real(real64) :: rrnmon      ! 1/nmon
    real(real64) :: Yfact       ! (nmon-2)*Y

    ! Thermodynamic scaling factors
    real(real64) :: scalem, emscale

    ! Boltzmann weight factors
    real(real64) :: bebelam, behbclam

    ! Geometry parameters
    real(real64) :: Rcoll2      ! Rcoll^2
    real(real64) :: zc2         ! zc1 + collsep

    ! Contact density normalization
    real(real64) :: cdnorm

    ! Solvent chemical potential
    real(real64) :: chemps
  end type computed_params_t

  ! ========================================================================
  ! TYPE: fields_t
  ! ========================================================================
  ! Encapsulates all density, potential, and work arrays for DFT calculation
  ! All arrays use 0-based indexing: (0:nrho, 0:nz) or (0:nrho, 0:nrho, 0:nz_hvec)
  ! ========================================================================
  type :: fields_t
    ! Density fields
    real(real64), allocatable :: fdmon(:, :)    ! Monomer density field ρ_m(r,z)
    real(real64), allocatable :: fem(:, :)      ! End-segment density field ρ_e(r,z)

    ! Boltzmann weight factors
    real(real64), allocatable :: ebelam(:, :)   ! End-segment Boltzmann factor exp(-βμ_e)
    real(real64), allocatable :: ehbclam(:, :)  ! Half-bond Boltzmann factor exp(-βμ_hb)

    ! Contact density and excess free energy
    real(real64), allocatable :: cdmonm(:, :)   ! Contact density at hard-sphere diameter
    real(real64), allocatable :: ae1(:, :)      ! Excess free energy (C1 term)
    real(real64), allocatable :: ae2(:, :)      ! Excess free energy (C2 term)
    real(real64), allocatable :: convp(:, :)    ! Convolution term for DFT functional

    ! External potentials
    real(real64), allocatable :: edu(:, :)      ! Lennard-Jones external potential exp(-βU_LJ)

    ! Lennard-Jones potential lookup table (3D)
    real(real64), allocatable :: hvec(:, :, :)  ! Precomputed LJ interaction integrals
  end type fields_t

  ! ========================================================================
  ! Module variables - structured input parameters
  ! ========================================================================
  type(input_params_t) :: input      ! User input parameters
  type(grid_params_t) :: grid        ! Grid dimensions and discretization
  type(computed_params_t) :: computed ! Derived parameters
  type(fields_t) :: fields           ! Density and potential fields

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
  subroutine allocate_arrays(flds, grd)
    type(fields_t), intent(inout) :: flds  ! Fields structure to allocate
    type(grid_params_t), intent(in) :: grd ! Grid parameters

    ! Local variables for dimensions
    integer(int32) :: nrho, nz, nz_hvec

    ! Calculate dimensions from grid parameters
    nrho = grd%mxrho + grd%kbl      ! Maximum rho grid points (including boundary)
    nz = grd%imitt + grd%ibl        ! Maximum z grid points (including boundary)
    nz_hvec = grd%nfack - 1         ! Maximum z for hvec (LJ potential table)

    ! Allocate 2D arrays with 0-based indexing
    allocate (flds%fdmon(0:nrho, 0:nz))
    allocate (flds%ebelam(0:nrho, 0:nz))
    allocate (flds%convp(0:nrho, 0:nz))
    allocate (flds%fem(0:nrho, 0:nz))
    allocate (flds%ehbclam(0:nrho, 0:nz))
    allocate (flds%cdmonm(0:nrho, 0:nz))
    allocate (flds%ae1(0:nrho, 0:nz))
    allocate (flds%ae2(0:nrho, 0:nz))
    allocate (flds%edu(0:nrho, 0:nz))

    ! Allocate 3D array (hvec has different z-dimension for LJ potential table)
    allocate (flds%hvec(0:nrho, 0:nrho, 0:nz_hvec))

    ! Initialize arrays to zero
    flds%fdmon = 0.d0
    flds%ebelam = 0.d0
    flds%convp = 0.d0
    flds%hvec = 0.d0
    flds%fem = 0.d0
    flds%ehbclam = 0.d0
    flds%cdmonm = 0.d0
    flds%ae1 = 0.d0
    flds%ae2 = 0.d0
    flds%edu = 0.d0

  end subroutine allocate_arrays

  ! ==========================================================================
  ! SUBROUTINE: deallocate_arrays
  ! ==========================================================================
  ! Deallocates all dynamic arrays
  ! Should be called before program termination for clean memory management
  ! ==========================================================================
  subroutine deallocate_arrays(flds)
    type(fields_t), intent(inout) :: flds  ! Fields structure to deallocate

    if (allocated(flds%fdmon)) deallocate (flds%fdmon)
    if (allocated(flds%ebelam)) deallocate (flds%ebelam)
    if (allocated(flds%convp)) deallocate (flds%convp)
    if (allocated(flds%hvec)) deallocate (flds%hvec)
    if (allocated(flds%fem)) deallocate (flds%fem)
    if (allocated(flds%ehbclam)) deallocate (flds%ehbclam)
    if (allocated(flds%cdmonm)) deallocate (flds%cdmonm)
    if (allocated(flds%ae1)) deallocate (flds%ae1)
    if (allocated(flds%ae2)) deallocate (flds%ae2)
    if (allocated(flds%edu)) deallocate (flds%edu)
  end subroutine deallocate_arrays

  !-----------------------------------------------------------------------------
  ! read_input_parameters - Read simulation parameters from input files
  !
  ! Reads all simulation parameters from input files and populates the
  ! input_params_t structure. Opens and closes the files internally.
  !
  ! Input files:
  !   - input.tsph: Main simulation parameters (bdm, nmon, grid spacing, etc.)
  !   - epfil: Lennard-Jones interaction parameter (epslj)
  !
  ! Outputs:
  !   inp - Populated input parameters structure
  !   bds - Solvent bulk density (computed as bdtot - bdm)
  !-----------------------------------------------------------------------------
  subroutine read_input_parameters(inp, bds)
    use iso_fortran_env, only: real64, int32
    implicit none

    type(input_params_t), intent(out) :: inp
    real(real64), intent(out) :: bds

    ! Local file units
    integer(int32) :: ins, iep

    ! Open input files
    open (newunit=ins, file='input.tsph', form='formatted', status='old')
    open (newunit=iep, file='epfil', form='formatted', status='old')

    ! Read simulation parameters from input file
    read (ins, *) inp%bdm         ! Monomer bulk density
    read (ins, *) inp%bdtot       ! Total bulk density
    bds = inp%bdtot - inp%bdm     ! Solvent bulk density
    read (ins, *) inp%nmon        ! Number of monomers per polymer chain
    read (ins, *) inp%dz          ! Grid spacing in z direction
    read (ins, *) inp%drho        ! Grid spacing in radial direction
    read (ins, *) inp%dphi        ! Angular grid spacing (input in units of pi)
    inp%dphi = PI*inp%dphi        ! Convert to radians
    read (ins, *) inp%Rcoll       ! Colloid radius
    read (ins, *) inp%zc1         ! Position of first colloid center
    read (ins, *) inp%collsep     ! Separation between colloid centers
    read (ins, *) inp%Rcyl        ! Cylinder radius (system boundary)
    read (ins, *) inp%ioimaxm     ! Maximum number of iterations
    read (ins, *) inp%dmm, inp%dms  ! Density mixing parameters (monomer, solvent)
    read (ins, *) inp%kread       ! Read initial guess from file (0=no, 1=yes)
    read (ins, *) inp%bl          ! Bond length
    read (ins, *) inp%dhs         ! Hard sphere diameter (monomer)
    read (ins, *) inp%dpphi       ! Angular grid spacing for potential calculation

    ! Read Lennard-Jones energy parameter
    read (iep, *) inp%epslj

    ! Close input files
    close (ins)
    close (iep)

    return
  end subroutine read_input_parameters

  ! ==========================================================================
  ! SUBROUTINE: initialize_grid_params
  ! ==========================================================================
  ! Computes grid parameters from input parameters
  ! Must be called after input parameters are read
  ! ==========================================================================
  subroutine initialize_grid_params(inp, grid)
    type(input_params_t), intent(in) :: inp
    type(grid_params_t), intent(out) :: grid

    ! Compute total z-grid extent
    grid%nfack = int(2.d0*(inp%zc1 + 0.5d0*inp%collsep)/inp%dz + 0.01d0)

    ! Set loop bounds
    grid%istart = 0
    grid%istp1 = 1
    grid%islut = grid%nfack
    grid%imitt = grid%nfack/2

    ! Solvent grid bounds (currently same as main grid)
    grid%istp1s = 1
    grid%isluts = grid%islut

    ! Discretized molecular dimensions
    grid%ism = int(inp%dhs/inp%dz + 0.01d0)
    grid%ksm = int(inp%dhs/inp%drho + 0.01d0)
    grid%ibl = int(inp%bl/inp%dz + 0.01d0)
    grid%kbl = int(inp%bl/inp%drho + 0.01d0)

    ! Angular discretization
    grid%nphi = int(PI/inp%dphi + 0.01d0)

    ! Copy nmon for convenience
    grid%nmon = inp%nmon

    ! mxrho computed later (depends on rdrho)
    grid%mxrho = 0
  end subroutine initialize_grid_params

  ! ==========================================================================
  ! SUBROUTINE: initialize_computed_params
  ! ==========================================================================
  ! Computes derived parameters from input and grid parameters
  ! Must be called after bulk thermodynamic calculations
  ! ==========================================================================
  subroutine initialize_computed_params(inp, grid, computed, bulk)
    type(input_params_t), intent(in) :: inp
    type(grid_params_t), intent(inout) :: grid
    type(computed_params_t), intent(out) :: computed
    type(bulk_properties_t), intent(in) :: bulk

    ! Grid spacing reciprocals
    computed%rdz = 1.d0/inp%dz
    computed%rdrho = 1.d0/inp%drho
    computed%rdphi = 1.d0/inp%dphi

    ! Now compute mxrho (depends on rdrho)
    grid%mxrho = int((inp%Rcyl - 1.d0)*computed%rdrho) + 1

    ! Molecular geometry
    computed%dhs2 = inp%dhs*inp%dhs
    computed%dhs3 = computed%dhs2*inp%dhs
    computed%rdhs3 = 1.d0/computed%dhs3
    computed%bl2 = inp%bl*inp%bl

    ! Grid integration helpers
    computed%dzrfp = inp%dz/(4.d0*PI)
    computed%twopidz = TWOPI*inp%dz

    ! Polymer chain parameters
    computed%rnmon = dble(inp%nmon)
    computed%rrnmon = 1.d0/computed%rnmon
    computed%Yfact = (computed%rnmon - 2.d0)*Y

    ! Thermodynamic scaling
    computed%scalem = bulk%chempp/(2.d0*computed%rnmon)
    computed%emscale = 2.d0*computed%scalem

    ! Boltzmann weights
    computed%bebelam = dexp(-bulk%emtrams + computed%emscale)
    computed%behbclam = dexp(-0.5d0*bulk%cmtrams + computed%scalem)

    ! Geometry
    computed%Rcoll2 = inp%Rcoll*inp%Rcoll
    computed%zc2 = inp%zc1 + inp%collsep

    ! Placeholder values
    computed%cdnorm = 0.d0  ! Set by CDFACT
    computed%chemps = 0.d0  ! Unused
  end subroutine initialize_computed_params

  ! ==========================================================================
  ! ==========================================================================
  ! SUBROUTINE: calculate_contact_density_normalization
  ! ==========================================================================
  ! Calculates the normalization constant (cdnorm) for the contact density
  ! functional. This is computed via a three-dimensional integral over the
  ! hard sphere volume using trapezoidal integration in cylindrical coordinates.
  ! ==========================================================================
  subroutine calculate_contact_density_normalization(inp, grd, comp, cos_phi_table)
    use iso_fortran_env, only: real64, int32
    implicit none

    ! Arguments
    type(input_params_t), intent(in) :: inp
    type(grid_params_t), intent(in) :: grd
    type(computed_params_t), intent(inout) :: comp
    real(real64), intent(in) :: cos_phi_table(:)

    ! Local variables
    integer(int32) :: iz, jz, iphi, irho, krhopmax, krhop
    real(real64) :: strho0, rho0, z, zpst, sume, zp, delz2, sumrhop
    real(real64) :: rhopmax, rho, rho02, rhop, rhomax2, fphi, phisum, rho2, fact, tcd
    strho0 = 0.5d0*inp%drho
    rho0 = strho0
    iz = 2*grd%ism
    z = 2.d0*inp%dhs - inp%dz
    zpst = z - inp%dhs - inp%dz
    sume = 0.d0
    zp = zpst
    ! Integrate over sphere of diameter dhs centered at test point
    ! Triple integration: z', rho', phi in cylindrical coordinates
    do jz = iz - grd%ism, iz + grd%ism
      zp = zp + inp%dz
      delz2 = (zp - z)**2
      sumrhop = 0.d0
      rhopmax = dsqrt(dabs(comp%dhs2 - delz2))
      krhopmax = nint(rhopmax*comp%rdrho)
      rho = rho0
      rho02 = rho0**2
      rhop = -0.5d0*inp%drho
      do krhop = 1, krhopmax
        rhop = rhop + inp%drho
        rhomax2 = rho02 + rhop*rhop
        fphi = 2.d0*rho0*rhop
        phisum = 0.d0
!$omp simd reduction(+:phisum)
        do iphi = 1, grd%nphi
!     Plus or minus sign doesn't matter for the value of the integral
          rho2 = rhomax2 - fphi*cos_phi_table(iphi)
          rho = dsqrt(rho2)
          irho = int(rho*comp%rdrho) + 1
          phisum = 1.d0 + phisum
        end do
!$omp end simd
        sumrhop = rhop*phisum*inp%dphi + sumrhop
      end do
      fact = 1.d0
      if (iabs(jz - iz) .eq. grd%ism) fact = 0.5d0
      sume = 2.d0*sumrhop*inp%drho*fact + sume
    end do
    tcd = 3.d0*sume*comp%dzrfp*comp%rdhs3
    comp%cdnorm = 1.d0/tcd
    return
  end subroutine calculate_contact_density_normalization

  ! ==========================================================================
  ! SUBROUTINE: calculate_contact_density
  ! ==========================================================================
  ! Calculates the contact density (cdmonm) at each grid point. The contact
  ! density is a weighted average of the monomer density around a sphere of
  ! diameter dhs, computed using 3D integration in cylindrical coordinates
  ! with angular averaging.
  ! ==========================================================================
  subroutine calculate_contact_density(inp, grd, comp, flds, cos_phi_table)
    use iso_fortran_env, only: real64, int32
    implicit none

    ! Arguments
    type(input_params_t), intent(in) :: inp
    type(grid_params_t), intent(in) :: grd
    type(computed_params_t), intent(in) :: comp
    type(fields_t), intent(inout) :: flds
    real(real64), intent(in) :: cos_phi_table(:)

    ! Local variables
    integer(int32) :: iz, jz, kz, iphi, irho, krhop, krhopmax
    real(real64) :: z, zpst, rho0, sume, zp, delz2, sumrhop, rhopmax, rho02
    real(real64) :: rhop, rhomax2, fphi, phisum, rho2, rho, fact

    ! Loop over all grid points to calculate contact density
!$omp parallel do private(z, zpst, rho0, kz, sume, zp, jz, delz2, sumrhop, rhopmax, krhopmax, rho02, rhop, rhomax2, fphi, phisum, iphi, rho2, rho, irho, fact) schedule(static)
    do iz = grd%istp1 + grd%ism, grd%imitt
      z = inp%dhs - 0.5d0*inp%dz + dble(iz - (grd%istp1 + grd%ism) + 1)*inp%dz
      zpst = z - inp%dhs - inp%dz
      rho0 = -0.5d0*inp%drho
      do kz = 1, grd%mxrho - grd%ksm
        rho0 = rho0 + inp%drho
        sume = 0.d0
        zp = zpst
        ! Integrate density over hard sphere volume
        do jz = iz - grd%ism, iz + grd%ism
          zp = zp + inp%dz
          delz2 = (zp - z)**2
          sumrhop = 0.d0
          rhopmax = dsqrt(dabs(comp%dhs2 - delz2))
          krhopmax = nint(rhopmax*comp%rdrho)
          rho02 = rho0**2
          rhop = -0.5d0*inp%drho
          do krhop = 1, krhopmax
            rhop = rhop + inp%drho
            rhomax2 = rho02 + rhop*rhop
            fphi = 2.d0*rho0*rhop
            phisum = 0.d0
!$omp simd reduction(+:phisum)
            do iphi = 1, grd%nphi
!     Plus or minus sign doesn't matter for the value of the integral
              rho2 = rhomax2 - fphi*cos_phi_table(iphi)
              rho = dsqrt(rho2)
              irho = int(rho*comp%rdrho) + 1
              phisum = flds%fdmon(irho, jz) + phisum
            end do
!$omp end simd
            sumrhop = rhop*phisum*inp%dphi + sumrhop
          end do
          fact = 1.d0
          if (iabs(jz - iz) .eq. grd%ism) fact = 0.5d0
          sume = 2.d0*sumrhop*inp%drho*fact + sume
        end do
        flds%cdmonm(kz, iz) = 3.d0*sume*comp%dzrfp*comp%cdnorm*comp%rdhs3
      end do
    end do
!$omp end parallel do
    return
  end subroutine calculate_contact_density

  ! ==========================================================================
  ! SUBROUTINE: calculate_excess_free_energy
  ! ==========================================================================
  ! Calculates the excess free energy arrays (ae1, ae2) and the convolution
  ! term (convp) from the contact density using the Carnahan-Starling equation
  ! of state. These quantities are used in the density functional expressions
  ! for the polymer chain propagators.
  ! ==========================================================================
  subroutine calculate_excess_free_energy(grd, comp, flds)
    use iso_fortran_env, only: real64, int32
    implicit none

    ! Arguments
    type(grid_params_t), intent(in) :: grd
    type(computed_params_t), intent(in) :: comp
    type(fields_t), intent(inout) :: flds

    ! Local variables
    integer(int32) :: iz, kz, jz
    real(real64) :: cdt, pcdt, xsi, rxsi, sqrxsi, flog, daex1, daex2

    ! Calculate excess free energy from contact density using Carnahan-Starling EOS
!$omp parallel do private(kz, cdt, pcdt, xsi, rxsi, sqrxsi, flog, daex1, daex2) schedule(static)
    do iz = grd%istp1 + 2*grd%ism, grd%imitt
      do kz = 1, grd%mxrho - grd%kbl
        cdt = flds%cdmonm(kz, iz)*comp%dhs3
        pcdt = PIS*cdt
        xsi = (1.d0 - pcdt)
        rxsi = 1.d0/xsi
        sqrxsi = rxsi*rxsi
        flog = dlog(xsi)
        flds%ae1(kz, iz) = -(c1 + 1.d0)*flog - 0.5d0*(AA1 + BB1*pcdt)*pcdt*sqrxsi
        flds%ae2(kz, iz) = -(c2 + 1.d0)*flog - 0.5d0*(AA2 + BB2*pcdt)*pcdt*sqrxsi
        daex1 = rxsi*(c1 + 1.d0 - 0.5d0*AA1*rxsi*(1.d0 + 2.d0*pcdt*rxsi) - &
                      BB1*pcdt*rxsi*(1.d0 + pcdt*rxsi))
        daex2 = rxsi*(c2 + 1.d0 - 0.5d0*AA2*rxsi*(1.d0 + 2.d0*pcdt*rxsi) - &
                      BB2*pcdt*rxsi*(1.d0 + pcdt*rxsi))
        flds%convp(kz, iz) = (Y*(flds%fdmon(kz, iz) - flds%fem(kz, iz))*(daex2 - daex1) + &
                              0.5d0*flds%fem(kz, iz)*daex2)*pis*comp%dhs3
      end do
    end do
!$omp end parallel do

!$omp parallel do private(jz, iz) schedule(static)
    do kz = 1, grd%mxrho - grd%kbl
      jz = grd%imitt + 1
      do iz = grd%imitt + 1, grd%imitt + grd%ibl
        jz = jz - 1
        flds%ae1(kz, iz) = flds%ae1(kz, jz)
        flds%ae2(kz, iz) = flds%ae2(kz, jz)
        flds%convp(kz, iz) = flds%convp(kz, jz)
      end do
    end do
!$omp end parallel do
    return
  end subroutine calculate_excess_free_energy

  ! ==========================================================================
  ! SUBROUTINE: calculate_boltzmann_factors
  ! ==========================================================================
  ! Calculates the Boltzmann weight factors for polymer end segments (ebelam)
  ! and middle/hinge segments (ehbclam) at each grid point. These factors
  ! include contributions from the excess free energy and the convolution
  ! integrals. Excludes regions inside colloids.
  ! ==========================================================================
  subroutine calculate_boltzmann_factors(inp, grd, comp, flds, cos_phi_table)
    use iso_fortran_env, only: real64, int32
    implicit none

    ! Arguments
    type(input_params_t), intent(in) :: inp
    type(grid_params_t), intent(in) :: grd
    type(computed_params_t), intent(in) :: comp
    type(fields_t), intent(inout) :: flds
    real(real64), intent(in) :: cos_phi_table(:)

    ! Local variables
    integer(int32) :: iz, kz, jstart, irho0min, krhop, krhopmax, jz, iphi, irho
    real(real64) :: z, zpst, diffz2, strho0, rho0, rho02, rt2, sume, zp
    real(real64) :: delz2, zpcsq, zpc2sq, sumrhop, rhopmax, rhop, rhomax2
    real(real64) :: fphi, phisum, rho2, rsq, rho, fact, trams, emtrams, cmtrams

    ! Set bulk values at boundaries
!$omp parallel do private(kz)
    do iz = 1, grd%ibl
      do kz = 1, grd%mxrho + grd%kbl
        flds%ebelam(kz, iz) = comp%bebelam
        flds%ehbclam(kz, iz) = comp%behbclam
      end do
    end do
!$omp end parallel do

    ! Calculate Boltzmann factors from convolution integrals
!$omp parallel do private(z, jstart, zpst, diffz2, irho0min, strho0, rho0, kz, rho02, rt2, sume, zp, jz, delz2, zpcsq, zpc2sq, sumrhop, rhopmax, krhopmax, rhop, rhomax2, fphi, phisum, iphi, rho2, rsq, rho, irho, fact, trams, emtrams, cmtrams) schedule(static)
    do iz = grd%ibl + 1, grd%imitt
      z = -0.5d0*inp%dz + dble(iz)*inp%dz
      jstart = iz - grd%ism
      zpst = z - inp%dhs - inp%dz
      diffz2 = (inp%zc1 - z)**2
      irho0min = 1
      strho0 = -0.5d0*inp%drho

      rho0 = strho0
      ! Loop over radial positions to calculate convolution integrals
      do kz = irho0min, grd%mxrho
        rho0 = rho0 + inp%drho
        rho02 = rho0**2

        ! Skip points inside first colloid
        rt2 = rho02 + (z - inp%zc1)**2
        if (rt2 .lt. comp%Rcoll2) then
          flds%ehbclam(kz, iz) = 0.d0
          flds%ebelam(kz, iz) = 0.d0
          cycle
        end if
        ! Skip points inside second colloid
        rt2 = rho02 + (z - comp%zc2)**2
        if (rt2 .lt. comp%Rcoll2) then
          flds%ehbclam(kz, iz) = 0.d0
          flds%ebelam(kz, iz) = 0.d0
          cycle
        end if

        ! Compute 3D convolution integral over hard sphere volume
        ! This gives the free energy contribution from local density variations
        sume = 0.d0
        zp = zpst
        ! Loop over z' within hard sphere diameter from current point
        do jz = jstart, iz + grd%ism
          zp = zp + inp%dz
          delz2 = (zp - z)**2
          zpcsq = (zp - inp%zc1)**2
          zpc2sq = (zp - comp%zc2)**2

          sumrhop = 0.d0
          ! Maximum rho' at this z' to stay within sphere of diameter dhs
          rhopmax = dsqrt(dabs(comp%dhs2 - delz2))
          krhopmax = nint(rhopmax*comp%rdrho)
          rho02 = rho0**2
          rhop = -0.5d0*inp%drho
          ! Loop over rho' (radial coordinate at integration point)
          do krhop = 1, krhopmax
            rhop = rhop + inp%drho
            rhomax2 = rho02 + rhop*rhop
            fphi = 2.d0*rho0*rhop
            phisum = 0.d0
            ! Loop over phi (angle between rho0 and rhop vectors)
            ! This completes the cylindrical coordinate integration
!$omp simd reduction(+:phisum)
            do iphi = 1, grd%nphi
              ! Plus or minus sign doesn't matter for the value of the integral
              ! Calculate rho at integration point using law of cosines
              rho2 = rhomax2 - fphi*cos_phi_table(iphi)
              ! Skip if integration point is inside first colloid
              rsq = rho2 + zpcsq
              if (rsq .lt. comp%Rcoll2) cycle
              ! Skip if integration point is inside second colloid
              rsq = rho2 + zpc2sq
              if (rsq .lt. comp%Rcoll2) cycle
              rho = dsqrt(rho2)
              irho = int(rho*comp%rdrho) + 1
              ! Sum convolution term (weighted density) over angular direction
              phisum = flds%convp(irho, jz) + phisum
            end do
!$omp end simd
            ! Integrate over rho': weight by rhop*dphi*drho (cylindrical volume element)
            sumrhop = rhop*phisum*inp%dphi + sumrhop
          end do
          ! Trapezoidal rule: half weight at boundaries
          fact = 1.d0
          if (iabs(jz - iz) .eq. grd%ism) fact = 0.5d0
          ! Integrate over z': weight by 2*drho*fact*dz
          sume = 2.d0*sumrhop*inp%drho*fact + sume
        end do
        ! Normalize convolution integral to get trams (local free energy contribution)
        trams = 3.d0*sume*comp%dzrfp*comp%cdnorm*comp%rdhs3

        ! Calculate Boltzmann factors from free energy contributions
        ! emtrams: end-segment contribution (half of middle segment)
        emtrams = trams + 0.5d0*flds%ae2(kz, iz)
        ! cmtrams: middle-segment contribution with Y factor for chain connectivity
        cmtrams = trams + Y*(flds%ae2(kz, iz) - flds%ae1(kz, iz))
        ! exp(-beta*F): Boltzmann weights with bulk chemical potential offset
        flds%ebelam(kz, iz) = dexp(-emtrams + comp%emscale)
        flds%ehbclam(kz, iz) = dexp(-0.5d0*(cmtrams) + comp%scalem)
      end do
    end do
!$omp end parallel do
    return
  end subroutine calculate_boltzmann_factors

  ! ==========================================================================
  ! SUBROUTINE: calculate_external_potential
  ! ==========================================================================
  ! Calculates the external potential contribution (edu) from Lennard-Jones
  ! interactions with the polymer solution. Computes exp(-U_LJ) where U_LJ is
  ! the interaction energy between a test particle at (rho,z) and the entire
  ! density field, using symmetry to include both colloids.
  ! ==========================================================================
  subroutine calculate_external_potential(inp, grd, comp, flds)
    use iso_fortran_env, only: real64, int32
    implicit none

    ! Arguments
    type(input_params_t), intent(in) :: inp
    type(grid_params_t), intent(in) :: grd
    type(computed_params_t), intent(in) :: comp
    type(fields_t), intent(inout) :: flds

    ! Local variables
    integer(int32) :: iz, kz, krho, kprho, ipz, itdz
    real(real64) :: z, sumpint, tz, tdz, sumrho

    ! Set boundary values to unity (no external potential at boundaries)
!$omp parallel do private(kz)
    do iz = 1, grd%ibl
      do kz = 1, grd%mxrho + grd%kbl
        flds%edu(kz, iz) = 1.d0
      end do
    end do
!$omp end parallel do

    ! Calculate Lennard-Jones potential energy at each grid point
!$omp parallel do private(z, krho, sumpint, tz, tdz, ipz, itdz, sumrho, kprho) schedule(static)
    do iz = grd%ibl + 1, grd%imitt
      z = (dble(iz) - 0.5d0)*inp%dz
      do krho = 1, grd%mxrho
        sumpint = 0.d0
        tz = inp%bl - 0.5d0*inp%dz
        ! Integrate LJ interaction with left half of density distribution
        do ipz = grd%ibl + 1, grd%imitt
          tz = tz + inp%dz
          tdz = z - tz
          itdz = nint(dabs(tdz*comp%rdz))
          sumrho = 0.d0
          do kprho = 1, grd%mxrho
            sumrho = sumrho + (flds%fdmon(kprho, ipz) - inp%bdm)*flds%hvec(kprho, krho, itdz)
          end do
          sumpint = 2.d0*sumrho*inp%drho + sumpint
        end do

        ! Integrate LJ interaction with right half (use symmetry)
        do ipz = grd%imitt + 1, grd%nfack - grd%ibl
          tz = tz + inp%dz
          tdz = z - tz
          itdz = nint(dabs(tdz*comp%rdz))
          sumrho = 0.d0
          do kprho = 1, grd%mxrho
            sumrho = sumrho + (flds%fdmon(kprho, grd%nfack + 1 - ipz) - inp%bdm)*flds%hvec(kprho, krho, itdz)
          end do
          sumpint = 2.d0*sumrho*inp%drho + sumpint
        end do
        sumpint = sumpint*inp%dz
        flds%edu(krho, iz) = dexp(-sumpint)
      end do
    end do
!$omp end parallel do
    return
  end subroutine calculate_external_potential

  ! ==========================================================================
  ! SUBROUTINE: calculate_lj_potential_table
  ! ==========================================================================
  ! Precomputes Lennard-Jones interaction potential on a grid for all
  ! distance combinations. This tabulation significantly speeds up the
  ! external potential calculation in EBDU.
  !
  ! The hvec array stores weighted LJ interaction integrals:
  !   hvec(rho', rho, dz) = rho' * integral[U_LJ(r) * dphi]
  ! where r is computed from (rho, rho', dz, phi) using cylindrical geometry.
  !
  ! Arguments:
  !   inp  - Input parameters (dz, drho, dpphi, epslj)
  !   grd  - Grid parameters (nfack, mxrho)
  !   comp - Computed parameters (dhs2, dhs3)
  !   flds - Fields structure (hvec will be filled)
  !
  ! Internal computations:
  !   - Calculates LJ coefficients alj and rlj from epslj and dhs3
  !   - Creates angular lookup table cos_pphi
  !   - Number of angular points (npphi) determined from dpphi
  ! ==========================================================================
  subroutine calculate_lj_potential_table(inp, grd, comp, flds)
    use iso_fortran_env, only: real64, int32
    implicit none

    ! Arguments
    type(input_params_t), intent(in) :: inp
    type(grid_params_t), intent(in) :: grd
    type(computed_params_t), intent(in) :: comp
    type(fields_t), intent(inout) :: flds

    ! Local variables
    integer(int32) :: itdz, iphi, krho, kprho, npphi
    real(real64) :: phi, tdz, tdzsq, rho, rhosq, use1, trho, trhosq
    real(real64) :: trmix, useful, pint, s2, dpphi_rad
    real(real64) :: alj, rlj
    real(real64), allocatable :: cos_pphi(:)

    ! Calculate Lennard-Jones coefficients from epslj and hard-sphere diameter
    alj = 4.d0*inp%epslj*comp%dhs3*comp%dhs3    ! Attractive (r^-6) coefficient
    rlj = 4.d0*inp%epslj*comp%dhs3**4           ! Repulsive (r^-12) coefficient

    ! Convert dpphi to radians and calculate number of angular points
    dpphi_rad = inp%dpphi*PI
    npphi = int(PI/dpphi_rad + 0.01d0)
    write (*, *) 'dpphi,npphi = ', inp%dpphi*PI, npphi

    ! Allocate and initialize cosine lookup table for angular integration
    allocate (cos_pphi(npphi))
    do iphi = 1, npphi
      phi = (dble(iphi) - 0.5d0)*dpphi_rad
      cos_pphi(iphi) = dcos(phi)
    end do

    ! Triple loop over rho, rho', z to compute pairwise LJ interaction integrals
    ! Loop order matches F77 for numerical consistency
!$omp parallel do private(tdz, tdzsq, rho, rhosq, use1, trho, trhosq, trmix, useful, pint, iphi, s2, krho, kprho) schedule(static)
    do itdz = 0, grd%nfack - 1
      tdz = -inp%dz + dble(itdz + 1)*inp%dz
      tdzsq = tdz*tdz
      rho = -0.5d0*inp%drho
      do krho = 1, grd%mxrho
        rho = rho + inp%drho
        rhosq = rho*rho
        use1 = tdzsq + rhosq
        trho = -0.5d0*inp%drho
        do kprho = 1, grd%mxrho
          trho = trho + inp%drho
          trhosq = trho*trho
          trmix = 2.d0*trho*rho
          useful = use1 + trhosq
          pint = 0.d0
          ! Angular integration for cylindrical geometry
!$omp simd reduction(+:pint)
          do iphi = 1, npphi
            s2 = useful - trmix*cos_pphi(iphi)
            if (s2 .gt. comp%dhs2) then
              ! Lennard-Jones potential: U(r) = 4*epsilon*[(sigma/r)^12 - (sigma/r)^6]
              pint = rlj/s2**6 - alj/s2**3 + pint
            end if
          end do
!$omp end simd
          flds%hvec(kprho, krho, itdz) = trho*pint*dpphi_rad
        end do
      end do
    end do
!$omp end parallel do

    ! Clean up temporary lookup table
    deallocate (cos_pphi)

    return
  end subroutine calculate_lj_potential_table

  ! ==========================================================================
  ! SUBROUTINE: initialize_density_fields
  ! ==========================================================================
  ! Initializes all density fields either from scratch (bulk with excluded
  ! volume) or by reading from a restart file. Also sets boundary conditions
  ! at all system boundaries.
  !
  ! For fresh start (kread=0):
  !   - Sets all interior points to bulk density
  !   - Zeros out colloid interiors (excluded volume)
  !
  ! For restart (kread=1):
  !   - Reads density fields from file unit ifc
  !
  ! Then applies boundary conditions:
  !   - z-boundaries: bulk values
  !   - Radial boundaries: bulk values
  !   - Midplane symmetry: mirror values from opposite side
  !
  ! Arguments:
  !   inp  - Input parameters (kread, dz, drho, bdm, zc1)
  !   grd  - Grid parameters (ibl, istp1, imitt, mxrho, kbl)
  !   comp - Computed parameters (zc2, Rcoll2, rrnmon, bebelam, behbclam)
  !   flds - Fields structure (output: fdmon, fem, ebelam, ehbclam, cdmonm)
  !   ifc  - File unit for reading restart data (if kread=1)
  ! ==========================================================================
  subroutine initialize_density_fields(inp, grd, comp, flds, ifc)
    use iso_fortran_env, only: real64, int32
    implicit none

    ! Arguments
    type(input_params_t), intent(in) :: inp
    type(grid_params_t), intent(in) :: grd
    type(computed_params_t), intent(in) :: comp
    type(fields_t), intent(inout) :: flds
    integer(int32), intent(in) :: ifc

    ! Local variables
    integer(int32) :: iz, jz, kz
    real(real64) :: z, z2, z22, rho, rt2, t1, t2

    ! Initialize density fields: either from scratch or read from file
    if (inp%kread .eq. 0) then
      ! Starting from bulk values with excluded volume for colloids
      z = -0.5d0*inp%dz
      ! Skip boundary region at z < 0
      do iz = 1, grd%ibl
        z = z + inp%dz
      end do
      ! Initialize all grid points to bulk values, then zero out colloid interiors
      do iz = grd%ibl + 1, grd%imitt
        z = z + inp%dz
        z2 = (z - inp%zc1)**2
        z22 = (z - comp%zc2)**2
        rho = -0.5d0*inp%drho
        ! Loop over radial positions
        do kz = 1, grd%mxrho
          rho = rho + inp%drho
          flds%fdmon(kz, iz) = inp%bdm
          flds%fem(kz, iz) = 2.d0*flds%fdmon(kz, iz)*comp%rrnmon
          flds%ebelam(kz, iz) = comp%bebelam
          flds%ehbclam(kz, iz) = comp%behbclam
          ! Check if point is inside first colloid
          rt2 = rho*rho + z2
          if (rt2 .lt. comp%Rcoll2) then
            flds%fdmon(kz, iz) = 0.d0
            flds%fem(kz, iz) = 0.d0
            flds%ebelam(kz, iz) = 0.d0
            flds%ehbclam(kz, iz) = 0.d0
          end if
          ! Check if point is inside second colloid
          rt2 = rho*rho + z22
          if (rt2 .lt. comp%Rcoll2) then
            flds%fdmon(kz, iz) = 0.d0
            flds%fem(kz, iz) = 0.d0
            flds%ebelam(kz, iz) = 0.d0
            flds%ehbclam(kz, iz) = 0.d0
          end if
        end do
      end do
    else
      ! Read initial guess from file (restart from previous calculation)
      rewind ifc
      do iz = grd%istp1, grd%imitt
        do kz = 1, grd%mxrho
          read (ifc, *) t1, t2, flds%fdmon(kz, iz), flds%fem(kz, iz)
        end do
      end do
    end if

    ! Set boundary conditions at z-boundaries (left edge)
    ! All densities set to bulk values
    do iz = grd%istp1, grd%ibl
      do kz = 1, grd%mxrho + grd%kbl
        flds%fdmon(kz, iz) = inp%bdm
        flds%fem(kz, iz) = 2.d0*flds%fdmon(kz, iz)*comp%rrnmon
        flds%ebelam(kz, iz) = comp%bebelam
        flds%ehbclam(kz, iz) = comp%behbclam
        flds%cdmonm(kz, iz) = inp%bdm
      end do
    end do

    ! Set boundary conditions at radial edge (outer cylinder boundary)
    ! All densities set to bulk values
    do iz = 1, grd%imitt
      do kz = grd%mxrho + 1, grd%mxrho + grd%kbl
        flds%fdmon(kz, iz) = inp%bdm
        flds%fem(kz, iz) = 2.d0*flds%fdmon(kz, iz)*comp%rrnmon
        flds%ebelam(kz, iz) = comp%bebelam
        flds%ehbclam(kz, iz) = comp%behbclam
        flds%cdmonm(kz, iz) = inp%bdm
      end do
    end do

    ! Apply symmetry boundary conditions at z = grd%imitt midplane
    jz = grd%imitt + 1
    do iz = grd%imitt + 1, grd%imitt + grd%ibl
      jz = jz - 1
      do kz = 1, grd%mxrho + grd%kbl
        flds%fdmon(kz, iz) = flds%fdmon(kz, jz)
        flds%fem(kz, iz) = flds%fem(kz, jz)
        flds%ebelam(kz, iz) = flds%ebelam(kz, jz)
        flds%ehbclam(kz, iz) = flds%ehbclam(kz, jz)
        flds%cdmonm(kz, iz) = inp%bdm
      end do
    end do

    return
  end subroutine initialize_density_fields

  ! ==========================================================================
  ! SUBROUTINE: initialize_boundary_excess_free_energy
  ! ==========================================================================
  ! Initializes excess free energy arrays (ae1, ae2, convp) at system
  ! boundaries where the density functional calculations need special
  ! treatment due to proximity to system edges.
  !
  ! The excess free energy is computed using the Carnahan-Starling equation
  ! of state for hard spheres, with corrections for polymer chain connectivity.
  !
  ! Two boundary regions are initialized:
  !   1. Near z-boundaries (left side): close to z = 0
  !   2. Near radial boundaries (outer edge): close to rho = Rcyl
  !
  ! Arguments:
  !   inp  - Input parameters (bdm, dhs)
  !   grd  - Grid parameters (istp1, ibl, mxrho, kbl, imitt)
  !   comp - Computed parameters (dhs3, rrnmon)
  !   flds - Fields structure (output: ae1, ae2, convp)
  ! ==========================================================================
  subroutine initialize_boundary_excess_free_energy(inp, grd, comp, flds)
    use iso_fortran_env, only: real64, int32
    implicit none

    ! Arguments
    type(input_params_t), intent(in) :: inp
    type(grid_params_t), intent(in) :: grd
    type(computed_params_t), intent(in) :: comp
    type(fields_t), intent(inout) :: flds

    ! Local variables
    integer(int32) :: iz, kz
    real(real64) :: cdt, pcdt, xsi, rxsi, sqrxsi, flog, daex1, daex2

    ! Initialize excess free energy arrays near z-boundaries (left side)
    ! These regions are near the system edge and require special treatment
    do iz = grd%istp1, grd%istp1 + 2*grd%ibl - 1
      do kz = 1, grd%mxrho + grd%kbl
        cdt = inp%bdm*comp%dhs3
        pcdt = PIS*cdt
        xsi = (1.d0 - pcdt)
        rxsi = 1.d0/xsi
        sqrxsi = rxsi*rxsi
        flog = dlog(xsi)
        flds%ae1(kz, iz) = -(C1 + 1.d0)*flog - 0.5d0*(AA1 + BB1*pcdt)*pcdt*sqrxsi
        flds%ae2(kz, iz) = -(C2 + 1.d0)*flog - 0.5d0*(AA2 + BB2*pcdt)*pcdt*sqrxsi
        daex1 = rxsi*(C1 + 1.d0 - 0.5d0*AA1*rxsi*(1.d0 + 2.d0*pcdt*rxsi) - &
                      BB1*pcdt*rxsi*(1.d0 + pcdt*rxsi))
        daex2 = rxsi*(C2 + 1.d0 - 0.5d0*AA2*rxsi*(1.d0 + 2.d0*pcdt*rxsi) - &
                      BB2*pcdt*rxsi*(1.d0 + pcdt*rxsi))
        flds%convp(kz, iz) = (Y*(inp%bdm - 2.d0*inp%bdm*comp%rrnmon)*(daex2 - daex1) + &
                              inp%bdm*comp%rrnmon*daex2)*PIS*comp%dhs3
      end do
    end do

    ! Initialize excess free energy arrays near radial boundaries (outer edge)
    do iz = grd%istp1 + 2*grd%ibl, grd%imitt + grd%ibl
      do kz = grd%mxrho - grd%kbl + 1, grd%mxrho + grd%kbl
        cdt = inp%bdm*comp%dhs3
        pcdt = PIS*cdt
        xsi = (1.d0 - pcdt)
        rxsi = 1.d0/xsi
        sqrxsi = rxsi*rxsi
        flog = dlog(xsi)
        flds%ae1(kz, iz) = -(C1 + 1.d0)*flog - 0.5d0*(AA1 + BB1*pcdt)*pcdt*sqrxsi
        flds%ae2(kz, iz) = -(C2 + 1.d0)*flog - 0.5d0*(AA2 + BB2*pcdt)*pcdt*sqrxsi
        daex1 = rxsi*(C1 + 1.d0 - 0.5d0*AA1*rxsi*(1.d0 + 2.d0*pcdt*rxsi) - &
                      BB1*pcdt*rxsi*(1.d0 + pcdt*rxsi))
        daex2 = rxsi*(C2 + 1.d0 - 0.5d0*AA2*rxsi*(1.d0 + 2.d0*pcdt*rxsi) - &
                      BB2*pcdt*rxsi*(1.d0 + pcdt*rxsi))
        flds%convp(kz, iz) = (Y*(inp%bdm - 2.d0*inp%bdm*comp%rrnmon)*(daex2 - daex1) + &
                              inp%bdm*comp%rrnmon*daex2)*PIS*comp%dhs3
      end do
    end do

    return
  end subroutine initialize_boundary_excess_free_energy

  ! ==========================================================================
  ! SUBROUTINE: calculate_bulk_properties
  ! ==========================================================================
  ! Calculates bulk thermodynamic properties using the Carnahan-Starling
  ! equation of state for hard spheres, with corrections for polymer chains.
  !
  ! Computes:
  !   - Packing fractions and related quantities
  !   - Excess free energy and its derivatives
  !   - Bulk pressure
  !   - Chemical potentials for polymer and solvent
  !   - Boltzmann weight factors (scalem, emscale)
  !   - Convolution terms for DFT functional
  !
  ! Also updates computed%scalem and computed%emscale which are needed
  ! throughout the calculation.
  !
  ! Arguments:
  !   inp   - Input parameters (bdm, nmon, dhs)
  !   comp  - Computed parameters (dhs3, rnmon, rrnmon, Yfact)
  !           Also updated: scalem, emscale (output)
  !   bulk  - Output: bulk thermodynamic properties
  ! ==========================================================================
  subroutine calculate_bulk_properties(inp, comp, bulk)
    use iso_fortran_env, only: real64, int32
    implicit none

    ! Arguments
    type(input_params_t), intent(in) :: inp
    type(computed_params_t), intent(inout) :: comp  ! scalem, emscale are outputs
    type(bulk_properties_t), intent(out) :: bulk

    ! Local variables
    real(real64) :: aeta, xsib, rxsib, rxsibsq
    real(real64) :: aex1, aex2, daex1, daex2, pdasum, bconvp, trams

    ! Bulk thermodynamic properties
    bulk%bdpol = inp%bdm/comp%rnmon  ! Polymer bulk density

    ! Hard sphere packing fraction and related quantities
    bulk%bdt = inp%bdm*comp%dhs3
    aeta = PIS*bulk%bdt
    xsib = 1.d0 - aeta
    rxsib = 1.d0/xsib
    rxsibsq = rxsib*rxsib

    ! Excess free energy terms (Carnahan-Starling)
    aex1 = -(C1 + 1.d0)*dlog(xsib) - &
           0.5d0*(AA1*PIS*bulk%bdt + BB1*(PIS*bulk%bdt)**2)*rxsibsq
    aex2 = -(C2 + 1.d0)*dlog(xsib) - &
           0.5d0*(AA2*PIS*bulk%bdt + BB2*(PIS*bulk%bdt)**2)*rxsibsq
    bulk%bFex = (inp%bdm - 2.d0*bulk%bdpol)*Y*(aex2 - aex1) + bulk%bdpol*aex2

    ! Derivatives of excess free energy
    daex1 = rxsib*(C1 + 1 - 0.5d0*(AA1 + 2.d0*BB1*aeta)*rxsib - &
                   aeta*(AA1 + BB1*aeta)*rxsibsq)
    daex2 = rxsib*(C2 + 1 - 0.5d0*(AA2 + 2.d0*BB2*aeta)*rxsib - &
                   aeta*(AA2 + BB2*aeta)*rxsibsq)
    pdasum = comp%Yfact*(daex2 - daex1) + daex2

    ! Bulk pressure and chemical potentials
    bulk%Pb = bulk%bdpol + bulk%bdpol*aeta*pdasum
    bulk%chempp = dlog(bulk%bdpol) + comp%Yfact*(aex2 - aex1) + aex2 + &
                  PIS*inp%bdm*pdasum*comp%dhs3
    comp%scalem = bulk%chempp/(2.d0*comp%rnmon)
    comp%emscale = 2.d0*comp%scalem

    ! Convolution terms for inhomogeneous density functional
    bconvp = (Y*(inp%bdm - 2.d0*inp%bdm*comp%rrnmon)*(daex2 - daex1) + &
              inp%bdm*comp%rrnmon*daex2)*PIS*comp%dhs3
    trams = bconvp
    bulk%emtrams = trams + 0.5d0*aex2
    bulk%cmtrams = trams + Y*(aex2 - aex1)

    return
  end subroutine calculate_bulk_properties

  !-----------------------------------------------------------------------------
  ! output_density_profiles - Write converged density profiles to output files
  !
  ! Writes density and propagator profiles to various output files:
  ! - fort.85: monomer and end-segment densities along z-axis at centerline (rho=0)
  ! - fort.89: propagator profiles along z-axis at centerline
  ! - fort.78: radially averaged density along z-axis (within radius 1.0)
  ! - fort.83: radial profile of densities at z = input%zc1
  ! - fort.87: radial profile of propagators at z = input%zc1
  !
  ! Opens files internally and closes them after writing.
  !
  ! Inputs:
  !   inp - Input parameters (dz, drho, zc1)
  !   grd - Grid parameters (istp1, imitt, mxrho)
  !   comp - Computed parameters (rdz)
  !   flds - Fields (fdmon, fem, ehbclam)
  !   c - Chain propagator array c(rho, z, segment)
  !-----------------------------------------------------------------------------
  subroutine output_density_profiles(inp, grd, comp, flds, c)
    use iso_fortran_env, only: real64, int32
    implicit none

    type(input_params_t), intent(in) :: inp
    type(grid_params_t), intent(in) :: grd
    type(computed_params_t), intent(in) :: comp
    type(fields_t), intent(in) :: flds
    real(real64), intent(in) :: c(0:, 0:, :)

    ! Local variables
    integer(int32) :: iz, i, kr, klm
    integer(int32) :: iout_zdens, iout_zprop, iout_zavg, iout_rdens, iout_rprop
    real(real64) :: z, rho, fsum

    ! Open output files for density profiles (use traditional fort.* names for compatibility)
    open (newunit=iout_zdens, file='fort.85', form='formatted', status='replace')
    open (newunit=iout_zprop, file='fort.89', form='formatted', status='replace')
    open (newunit=iout_zavg, file='fort.78', form='formatted', status='replace')
    open (newunit=iout_rdens, file='fort.83', form='formatted', status='replace')
    open (newunit=iout_rprop, file='fort.87', form='formatted', status='replace')

    z = -0.5d0*inp%dz

    ! Write density profiles along z-axis at rho=0 (centerline)
    do iz = grd%istp1, grd%imitt
      z = z + inp%dz
      ! iout_zdens: monomer and end-segment densities at centerline
      write (iout_zdens, *) z, flds%fdmon(1, iz), flds%fem(1, iz)
      ! iout_zprop: propagators for segments 1,3,5,9 and ehbclam at centerline
      write (iout_zprop, '(6f14.7)') z, c(1, iz, 1), c(1, iz, 3), c(1, iz, 5), &
        flds%ehbclam(1, iz), c(1, iz, 9)

      ! Radially integrate density within radius 1.0 to get average
      fsum = 0.d0
      klm = nint(1.d0/inp%drho)
      rho = -0.5d0*inp%drho
      do i = 1, klm
        rho = rho + inp%drho
        fsum = fsum + flds%fdmon(i, iz)*2.d0*PI*rho
      end do
      ! iout_zavg: z-position and radially averaged density
      write (iout_zavg, *) z, fsum*inp%drho/(PI*1.d0**2)
    end do

    ! Write radial profiles at z = inp%zc1 (first colloid center position)
    rewind iout_rdens
    rewind iout_rprop
    iz = int(inp%zc1*comp%rdz) + 1
    rho = -0.5d0*inp%drho
    do kr = 1, grd%mxrho
      rho = rho + inp%drho
      ! iout_rprop: radial profile of propagators for segments 1,3,5,9 and ehbclam
      write (iout_rprop, '(6f14.7)') rho, c(kr, iz, 1), c(kr, iz, 3), c(kr, iz, 5), &
        flds%ehbclam(kr, iz), c(kr, iz, 9)
      ! iout_rdens: radial profile of monomer and end-segment densities
      write (iout_rdens, *) rho, flds%fdmon(kr, iz), flds%fem(kr, iz)
    end do

    ! Close output files
    close (iout_zdens)
    close (iout_zprop)
    close (iout_zavg)
    close (iout_rdens)
    close (iout_rprop)

    return
  end subroutine output_density_profiles

  !-----------------------------------------------------------------------------
  ! calculate_grand_potential - Calculate grand potential (thermodynamic free energy)
  !
  ! Computes the grand potential Omega = F - mu*N, which measures the
  ! thermodynamic cost of the inhomogeneous density distribution relative
  ! to the bulk state. Integrates the grand potential density over the
  ! system volume, excluding colloid interiors.
  !
  ! Returns two formulations:
  ! - aW: Full grand potential including all terms
  ! - bW: Alternative formulation for comparison
  !
  ! Inputs:
  !   inp - Input parameters (bdm, dz, drho, zc1)
  !   grd - Grid parameters (istp1, imitt, mxrho)
  !   comp - Computed parameters (Rcoll2, emscale, scalem, rrnmon)
  !   flds - Fields (fdmon, fem, ebelam, ehbclam, ae1, ae2, edu)
  !   bulk - Bulk thermodynamic properties
  !
  ! Outputs:
  !   aW - Grand potential (primary formulation)
  !   bW - Grand potential (alternative formulation)
  !-----------------------------------------------------------------------------
  subroutine calculate_grand_potential(inp, grd, comp, flds, bulk, aW, bW)
    use iso_fortran_env, only: real64, int32
    implicit none

    type(input_params_t), intent(in) :: inp
    type(grid_params_t), intent(in) :: grd
    type(computed_params_t), intent(in) :: comp
    type(fields_t), intent(in) :: flds
    type(bulk_properties_t), intent(in) :: bulk
    real(real64), intent(out) :: aW, bW

    ! Local variables
    integer(int32) :: iz, kz
    real(real64) :: z, rho, rsq, diffz2
    real(real64) :: bfde, bfdc, asumW, bsumW
    real(real64) :: arsum, brsum
    real(real64) :: fdm, fde, fdc, belamb, bclamb, Fex, eexc

    bfde = 2.d0*bulk%bdpol  ! Bulk end-segment density
    bfdc = inp%bdm - bfde   ! Bulk internal-segment density
    asumW = 0.d0
    bsumW = 0.d0
    z = -0.5d0*inp%dz

    ! Integrate grand potential density over system volume
    do iz = grd%istp1, grd%imitt
      z = z + inp%dz
      arsum = 0.d0
      brsum = 0.d0
      diffz2 = (z - inp%zc1)**2
      rho = -0.5d0*inp%drho
      do kz = 1, grd%mxrho
        rho = rho + inp%drho
        rsq = rho*rho + diffz2
        fdm = flds%fdmon(kz, iz)

        ! Only integrate outside colloid volume
        if (rsq .ge. comp%Rcoll2) then
          ! Chemical potential contributions
          belamb = dlog(flds%ebelam(kz, iz)) - comp%emscale
          bclamb = 2.d0*(dlog(flds%ehbclam(kz, iz)) - comp%scalem)
          fde = flds%fem(kz, iz)
          fdc = fdm - fde
          ! Excess free energy from hard-sphere interactions
          Fex = fdc*Y*(flds%ae2(kz, iz) - flds%ae1(kz, iz)) + 0.5d0*fde*flds%ae2(kz, iz)

          ! Grand potential density omega(r) = f(r) - mu*rho(r)
          ! where f(r) is Helmholtz free energy density
          arsum = &
            rho*(fdc*bclamb + bfdc*bulk%cmtrams + fde*belamb + bfde*bulk%emtrams + &
                 bulk%bdpol - fdm*comp%rrnmon + Fex - bulk%bFex) + arsum
          brsum = &
            rho*(fdc*bclamb + fde*belamb - fdm*comp%rrnmon + Fex - bulk%bFex) + brsum
        end if

        ! Add Lennard-Jones contribution to grand potential
        eexc = -dlog(flds%edu(kz, iz))
        arsum = arsum - 0.5d0*rho*eexc*(fdm + inp%bdm)
        brsum = brsum + rho*(0.5d0*(fdm - inp%bdm)*eexc - fdm*eexc)
      end do
      ! Integrate radially: multiply by 2*pi*rho*inp%drho
      asumW = 2.d0*PI*arsum*inp%drho + asumW
      bsumW = 2.d0*PI*brsum*inp%drho + bsumW
    end do
    ! Integrate along z-axis: multiply by inp%dz
    asumW = asumW*inp%dz
    bsumW = bsumW*inp%dz
    ! Factor of 2 accounts for both halves of symmetric system
    aW = 2.d0*asumW
    bW = 2.d0*bsumW

    write (*, *)
    write (*, *) 'aW = ', aW
    write (*, *) 'bW = ', bW

    return
  end subroutine calculate_grand_potential

  !-----------------------------------------------------------------------------
  ! calculate_colloid_forces - Calculate forces on colloid from contact density
  !
  ! Integrates the contact density over the colloid surface to compute the
  ! net force on the colloid. Uses multiple integration methods:
  ! 1. Integration in z-direction (hemispheres)
  ! 2. Integration in rho-direction (alternative method)
  ! 3. Integration over cos(theta) coordinate
  !
  ! The contact density at the colloid surface is obtained by quadratic
  ! (Lagrange) interpolation from nearby grid points.
  !
  ! Inputs:
  !   inp - Input parameters (zc1, Rcoll, dz, drho)
  !   comp - Computed parameters (rdz, Rcoll2)
  !   flds - Fields (fdmon - monomer density field)
  !
  ! Outputs:
  !   rcliffF - Force from z-hemisphere integration
  !   ctF - Force from cos(theta) integration
  !   ch2 - Alternative integral for comparison
  !-----------------------------------------------------------------------------
  subroutine calculate_colloid_forces(inp, comp, flds, rcliffF, ctF, ch2)
    use iso_fortran_env, only: real64, int32
    implicit none

    type(input_params_t), intent(in) :: inp
    type(computed_params_t), intent(in) :: comp
    type(fields_t), intent(in) :: flds
    real(real64), intent(out) :: rcliffF, ctF, ch2

    ! Local variables
    integer(int32) :: iz, izmin, izmax, izc1, irho, irhomax, krho, kct, ict
    real(real64) :: z, zmin, zmax, rho, rhon, rhoc, rhomax, Rc
    real(real64) :: zsq, rhosq, deltazc
    real(real64) :: x, x1, x2, x3, y1, y2, y3
    real(real64) :: fdc, fdcn, fdcm1, fdcp1, fk
    real(real64) :: ctheta, ct, ctn, ctp
    real(real64) :: rhoFo, rcliffFo, rhoFi, rcliffFi
    real(real64) :: zFi, zFo, chi, cho, cliffFi, cliffFo
    real(real64) :: cckoll, ckoll, ckk, ccckoll, bordekoll, ccc
    real(real64) :: add, t, tn, th
    real(real64), dimension(0:1000) :: cdens, ctvec

    ! Determine integration limits for first colloid (centered at inp%zc1)
    izmin = nint((inp%zc1 + 0.5d0*inp%dz - inp%Rcoll)*comp%rdz + 0.5d0)
    zmin = (dfloat(izmin) - 0.5d0)*inp%dz
    izmax = nint((inp%zc1 - 0.5d0*inp%dz + inp%Rcoll)*comp%rdz + 0.5d0)
    zmax = (dfloat(izmax) - 0.5d0)*inp%dz
    izc1 = nint(inp%zc1*comp%rdz + 0.5d0)
    write (*, *) 'zmin,inp%zc1,zmax = ', zmin, inp%zc1, zmax
    write (*, *) 'izmin,izc1,izmax = ', izmin, izc1, izmax
    write (*, *) dfloat(izmin)*inp%dz - 0.5d0*inp%dz, dfloat(izmax)*inp%dz - 0.5d0*inp%dz
    write (*, *) dfloat(izc1)*inp%dz - 0.5d0*inp%dz
    ict = 0

    ! Calculate force on outer hemisphere (z < inp%zc1) of first colloid
    rhoFo = 0.d0
    rcliffFo = 0.d0
    z = zmin - inp%dz
    do iz = izmin, izc1 - 1
      z = z + inp%dz
      zsq = (z - inp%zc1)**2
      ! Only process z-slices that intersect the colloid
      if (zsq .le. comp%Rcoll2) then
        rho = -0.5d0*inp%drho
        irho = 0
        ! Find first grid point outside colloid at this z
        do
          rho = rho + inp%drho
          irho = irho + 1
          if ((rho*rho + zsq) .gt. comp%Rcoll2) exit
        end do
        Rc = dsqrt(rho*rho + zsq)
        rhoc = dsqrt(comp%Rcoll2 - zsq)

        ! Quadratic interpolation to get density at exact colloid surface
        ! Use 3 points near boundary (irho, irho+1, irho+2)
        if (dabs(flds%fdmon(irho, iz)) .gt. 0.00000001d0) then
          y3 = flds%fdmon(irho, iz)
          y2 = flds%fdmon(irho + 1, iz)
          y1 = flds%fdmon(irho + 2, iz)
          x3 = rho
          x2 = rho + inp%drho
          x1 = rho + 2.d0*inp%drho
        else
          y3 = flds%fdmon(irho + 1, iz)
          y2 = flds%fdmon(irho + 2, iz)
          y1 = flds%fdmon(irho + 3, iz)
          x3 = rho + inp%drho
          x2 = rho + 2.d0*inp%drho
          x1 = rho + 3.d0*inp%drho
          write (*, *) 'TJOHO!'
        end if

        ! Lagrange interpolation to get density at colloid surface
        x = rhoc
        fdc = y1*(x - x2)*(x - x3)/((x1 - x2)*(x1 - x3)) + &
              y2*(x - x1)*(x - x3)/((x2 - x1)*(x2 - x3)) + &
              y3*(x - x1)*(x - x2)/((x3 - x1)*(x3 - x2))
        ! cos(theta) = (z - inp%zc1)/inp%Rcoll for surface normal direction
        ctheta = (z - inp%zc1)/inp%Rcoll
        rhoFo = 2.d0*PI*rhoc*ctheta*fdc + rhoFo
        rcliffFo = 2.d0*PI*ctheta*fdc + rcliffFo
        ict = ict + 1
        ctvec(ict) = ctheta
        cdens(ict) = fdc
      end if
    end do
    write (*, *) 'rcliffFo = ', inp%Rcoll*rcliffFo*inp%dz
    write (*, *) 'z = ', z

    ! Calculate force on inner hemisphere (z > inp%zc1) of first colloid
    rhoFi = 0.d0
    rcliffFi = 0.d0
    z = inp%zc1 - 0.5d0*inp%dz
    do iz = izc1, izmax
      z = z + inp%dz
      zsq = (z - inp%zc1)**2
      ! Only process z-slices that intersect the colloid
      if (zsq .le. comp%Rcoll2) then
        rho = -0.5d0*inp%drho
        irho = 0
        ! Find first grid point outside colloid at this z
        do
          rho = rho + inp%drho
          irho = irho + 1
          if ((rho*rho + zsq) .gt. comp%Rcoll2) exit
        end do
        Rc = dsqrt(rho*rho + zsq)
        rhoc = dsqrt(comp%Rcoll2 - zsq)

        if (dabs(flds%fdmon(irho, iz)) .gt. 0.00000001d0) then
          y3 = flds%fdmon(irho, iz)
          y2 = flds%fdmon(irho + 1, iz)
          y1 = flds%fdmon(irho + 2, iz)
          x3 = rho
          x2 = rho + inp%drho
          x1 = rho + 2.d0*inp%drho
        else
          y3 = flds%fdmon(irho + 1, iz)
          y2 = flds%fdmon(irho + 2, iz)
          y1 = flds%fdmon(irho + 3, iz)
          x3 = rho + inp%drho
          x2 = rho + 2.d0*inp%drho
          x1 = rho + 3.d0*inp%drho
          write (*, *) 'TJOHO!!!!', flds%fdmon(irho, iz), rho
        end if

        x = rhoc
        fdc = y1*(x - x2)*(x - x3)/((x1 - x2)*(x1 - x3)) + &
              y2*(x - x1)*(x - x3)/((x2 - x1)*(x2 - x3)) + &
              y3*(x - x1)*(x - x2)/((x3 - x1)*(x3 - x2))
        ctheta = (z - inp%zc1)/inp%Rcoll
        rhoFi = 2.d0*PI*rhoc*ctheta*fdc + rhoFi
        rcliffFi = 2.d0*PI*ctheta*fdc + rcliffFi
        ict = ict + 1
        ctvec(ict) = ctheta
        cdens(ict) = fdc
      end if
    end do
    write (*, *) 'rcliffFi = ', inp%Rcoll*rcliffFi*inp%dz
    write (*, *)
    write (*, *) 'rcliffF = ', inp%Rcoll*(rcliffFi + rcliffFo)*inp%dz
    write (*, *)
    write (*, *) 'z = ', z

    ! Integrate force over colloid surface using cos(theta) as coordinate
    ! Extrapolate density to poles (theta = ±1) using linear interpolation
    ctF = 0.d0
    fdcm1 = &
      cdens(1) + (cdens(2) - cdens(1))*(-1.d0 - ctvec(1))/(ctvec(2) - ctvec(1))
    write (*, *) 'fdcm1 = ', fdcm1
    write (*, *) 'cdens(1),cdens(2) = ', cdens(1), cdens(2)
    cdens(0) = fdcm1
    ctvec(0) = -1.d0
    fdcp1 = &
      cdens(ict) + &
      (cdens(ict) - cdens(ict - 1))*(1.d0 - ctvec(ict))/(ctvec(ict) - &
                                                         ctvec(ict - 1))
    write (*, *) 'fdcp1 = ', fdcp1
    write (*, *) 'cdens(ict),cdens(ict-1) = ', cdens(ict), cdens(ict - 1)
    cdens(ict + 1) = fdcp1
    ctvec(ict + 1) = 1.d0

    ! Trapezoidal integration over cos(theta) from -1 to +1
    ch2 = 0.d0
    do kct = 0, ict
      ct = ctvec(kct)
      ctn = ctvec(kct + 1)
      fdc = cdens(kct)
      fdcn = cdens(kct + 1)
      fk = (fdcn - fdc)/(ctn - ct)
      ctF = 0.5d0*(fdc - fk*ct)*(ctn**2 - ct**2) + fk*(ctn**3 - ct**3)/3.d0 + ctF
      ch2 = 0.25d0*(fdc + fdcn)*(ctn**2 - ct**2) + ch2
    end do
    ctF = 2.d0*PI*comp%Rcoll2*ctF
    ch2 = 2.d0*PI*comp%Rcoll2*ch2
    write (*, *)
    write (*, *) 'ctF = ', ctF
    write (*, *)
    write (*, *) 'ch2 = ', ch2
    write (*, *)

    ! Alternative force calculation: integrate over rho slices at constant z

    ! Determine maximum radial index inside sphere
    irhomax = nint(inp%Rcoll*comp%rdrho + 1.d0)
    rhomax = (dfloat(irhomax) - 0.5d0)*inp%drho
    ! Stay inside the sphere to avoid boundary issues
    irhomax = irhomax - 1
    rhomax = rhomax - inp%drho
    write (*, *) 'rhomax,irhomax = ', rhomax, irhomax
    ict = 0
    zFo = 0.d0
    cho = 0.d0
    cliffFo = 0.d0
    rho = -0.5d0*inp%dz

    ! Loop over radial slices from center outward (outer hemisphere in z)
    do irho = 1, irhomax
      rho = rho + inp%drho
      rhosq = rho*rho
      ! Only process rho values that intersect the colloid
      if (rhosq .le. comp%Rcoll2) then
        z = inp%zc1 + 0.5d0*inp%dz
        iz = izc1
        ! Find first grid point outside colloid at this rho (moving down in z)
        do
          z = z - inp%dz
          iz = iz - 1
          zsq = (z - inp%zc1)**2
          if ((rhosq + zsq) .gt. comp%Rcoll2) exit
        end do
        Rc = dsqrt(rhosq + zsq)
        deltazc = dsqrt(comp%Rcoll2 - rhosq)

        ! Ensure we have enough points for 3-point interpolation
        ! Skip if iz is too close to lower bound
        if (iz - 3 .lt. 0) cycle

        ! Quadratic interpolation in z-direction to get surface density
        if (dabs(flds%fdmon(irho, iz)) .gt. 0.00000001d0) then
          y3 = flds%fdmon(irho, iz)
          y2 = flds%fdmon(irho, iz - 1)
          y1 = flds%fdmon(irho, iz - 2)
          x3 = dabs(z - inp%zc1)
          x2 = x3 + inp%dz
          x1 = x3 + 2.d0*inp%dz
        else
          y3 = flds%fdmon(irho, iz - 1)
          y2 = flds%fdmon(irho, iz - 2)
          y1 = flds%fdmon(irho, iz - 3)
          x3 = dabs(z - inp%zc1) + inp%dz
          x2 = x3 + inp%dz
          x1 = x3 + 2.d0*inp%dz
          write (*, *) 'TJOHO1!'
        end if

        x = deltazc
        fdc = y1*(x - x2)*(x - x3)/((x1 - x2)*(x1 - x3)) + &
              y2*(x - x1)*(x - x3)/((x2 - x1)*(x2 - x3)) + &
              y3*(x - x1)*(x - x2)/((x3 - x1)*(x3 - x2))
        ctheta = -deltazc/inp%Rcoll
        ict = ict + 1
        ctvec(ict) = ctheta
        cdens(ict) = fdc
      end if
    end do

    ! Loop over radial slices in reverse (inner hemisphere in z)
    rho = rho + inp%drho
    irho = irhomax + 1
    zFi = 0.d0
    chi = 0.d0
    cliffFi = 0.d0
    do krho = 1, irhomax
      irho = irho - 1
      rho = rho - inp%drho
      rhosq = rho*rho
      ! Only process rho values that intersect the colloid
      if (rhosq .le. comp%Rcoll2) then
        z = inp%zc1 - 0.5d0*inp%dz
        iz = izc1 - 1
        ! Find first grid point outside colloid at this rho (moving up in z)
        do
          z = z + inp%dz
          iz = iz + 1
          zsq = (z - inp%zc1)**2
          if ((rhosq + zsq) .gt. comp%Rcoll2) exit
        end do
        Rc = dsqrt(rhosq + zsq)
        deltazc = dsqrt(comp%Rcoll2 - rhosq)

        ! Ensure we have enough points for 3-point interpolation
        ! Skip if iz is too close to upper bound
        if (iz + 3 .gt. ubound(flds%fdmon, 2)) cycle

        if (dabs(flds%fdmon(irho, iz)) .gt. 0.00000001d0) then
          y3 = flds%fdmon(irho, iz)
          y2 = flds%fdmon(irho, iz + 1)
          y1 = flds%fdmon(irho, iz + 2)
          x3 = dabs(z - inp%zc1)
          x2 = x3 + inp%dz
          x1 = x3 + 2.d0*inp%dz
        else
          y3 = flds%fdmon(irho, iz + 1)
          y2 = flds%fdmon(irho, iz + 2)
          y1 = flds%fdmon(irho, iz + 3)
          x3 = dabs(z - inp%zc1) + inp%dz
          x2 = x3 + inp%dz
          x1 = x3 + 2.d0*inp%dz
          write (*, *) 'TJOHO2!'
        end if

        x = deltazc
        fdc = y1*(x - x2)*(x - x3)/((x1 - x2)*(x1 - x3)) + &
              y2*(x - x1)*(x - x3)/((x2 - x1)*(x2 - x3)) + &
              y3*(x - x1)*(x - x2)/((x3 - x1)*(x3 - x2))
        ctheta = deltazc/inp%Rcoll
        zFi = 2.d0*PI*rho*ctheta*fdc + zFi
        chi = 2.d0*PI*rho*ctheta + chi
        cliffFi = 2.d0*PI*ctheta*fdc + cliffFi
        ict = ict + 1
        ctvec(ict) = ctheta
        cdens(ict) = fdc
      end if
    end do

    ! Second integration over cos(theta) from alternative method
    ! Extrapolate to poles as before
    ctF = 0.d0
    th = 1.5d0
    fdcm1 = &
      cdens(1) + (cdens(2) - cdens(1))*(-1.d0 - ctvec(1))/(ctvec(2) - ctvec(1))
    write (*, *) 'fdcm1 = ', fdcm1
    write (*, *) 'cdens(1),cdens(2) = ', cdens(1), cdens(2)
    cdens(0) = fdcm1
    ctvec(0) = -1.d0
    fdcp1 = &
      cdens(ict) + &
      (cdens(ict) - cdens(ict - 1))*(1.d0 - ctvec(ict))/(ctvec(ict) - &
                                                         ctvec(ict - 1))
    write (*, *) 'fdcp1 = ', fdcp1
    write (*, *) 'cdens(ict),cdens(ict-1) = ', cdens(ict), cdens(ict - 1)
    ctvec(ict + 1) = 1.d0
    cdens(ict + 1) = fdcp1

    ! Initialize various force and integral accumulators
    ckoll = 0.d0
    cckoll = 0.d0
    ccckoll = 0.d0
    ch2 = 0.d0
    ccc = 0.d0
    bordekoll = 0.d0

    ! Loop over all theta intervals to compute multiple integrals
    do kct = 0, ict
      ct = ctvec(kct)
      ctn = ctvec(kct + 1)
      fdc = cdens(kct)
      fdcn = cdens(kct + 1)
      fk = (fdcn - fdc)/(ctn - ct)

      ! Trapezoidal integration of force in cos(theta) coordinate
      ctF = 0.5d0*(fdc - fk*ct)*(ctn**2 - ct**2) + fk*(ctn**3 - ct**3)/3.d0 + ctF
      ch2 = 0.25d0*(fdc + fdcn)*(ctn**2 - ct**2) + ch2

      ! Additional integral using sin(theta)^3 = (1 - cos^2(theta))^(3/2)
      tn = -(1.d0 - ctn*ctn)**1.5d0
      t = -(1.d0 - ct*ct)**1.5d0
      add = 0.5d0*(fdc + fdcn)*(tn - t)/3.d0
      bordekoll = add + bordekoll

      ! Various auxiliary integrals for force calculation checks
      if (kct .gt. 0) then
        cckoll = ct*fdc + cckoll
        if (ct .lt. 0.d0) then
          ccckoll = -fdc*ct*ct*(ctn - ct) + ccckoll
        else
          ccckoll = fdc*ct*ct*(ct - ctp) + ccckoll
        end if
        ckoll = ct*fdc*dsqrt(1.d0 - ct*ct) + ckoll
        ckk = dabs(ct)*dsqrt(1.d0 - ct*ct) + ckk
        ctp = ct
      end if

      ! Convert cos(theta) to rho coordinate: rho = R*sin(theta) = R*sqrt(1-cos^2)
      rho = inp%Rcoll*dsqrt(1.d0 - ct*ct)
      rhon = inp%Rcoll*dsqrt(1.d0 - ctn*ctn)
      if (kct .eq. 0) rho = inp%Rcoll
      if (kct .eq. ict) rhon = inp%Rcoll
      if (dabs(rhon - rho) .lt. 0.00000001d0) then
        fk = 0.d0
      else
        fk = (fdcn - fdc)/(rhon - rho)
      end if
      ccc = (fdc - fk*rho)*(rhon - rho) + 0.5d0*fk*(rhon*rhon - rho*rho) + ccc
    end do

    ctF = 2.d0*PI*comp%Rcoll2*ctF
    ch2 = 2.d0*PI*comp%Rcoll2*ch2
    cckoll = 2.d0*PI*cckoll*inp%drho
    ckoll = 2.d0*PI*ckoll*inp%Rcoll*inp%drho
    ckk = 2.d0*PI*ckk*inp%Rcoll*inp%drho
    ccc = 2.d0*PI*ccc
    write (*, *)
    write (*, *) 'ctF = ', ctF
    write (*, *)
    write (*, *) 'ch2 = ', ch2

    ! Return final values
    rcliffF = inp%Rcoll*(rcliffFi + rcliffFo)*inp%dz

    return
  end subroutine calculate_colloid_forces

  ! ==========================================================================
  ! SUBROUTINE: apply_boundary_conditions
  ! ==========================================================================
  ! Applies boundary conditions to field variables during SCF iteration
  !
  ! 1. Radial boundary (rho > Rcyl): Sets fields to bulk values
  ! 2. Z-symmetry at midplane: Mirrors fields across z = imitt
  !
  ! The z-symmetry is applied twice to ensure full consistency
  ! ==========================================================================
  subroutine apply_boundary_conditions(grd, comp, flds)
    use iso_fortran_env, only: real64, int32
    implicit none

    ! Arguments
    type(grid_params_t), intent(in) :: grd
    type(computed_params_t), intent(in) :: comp
    type(fields_t), intent(inout) :: flds

    ! Local variables
    integer(int32) :: iz, kz, jz

    ! Apply boundary conditions at outer radial edge (rho > Rcyl)
    ! Set to bulk values since density should approach bulk far from colloids
    do iz = grd%istp1, grd%imitt
    do kz = grd%mxrho + 1, grd%mxrho + grd%kbl
      flds%ebelam(kz, iz) = comp%bebelam
      flds%ehbclam(kz, iz) = comp%behbclam
      flds%edu(kz, iz) = 1.d0
    end do
    end do

    ! Apply symmetry boundary conditions at z = imitt (midplane between colloids)
    ! The system is symmetric about the midplane, so mirror the field values
    ! This loop copies from iz = imitt down to iz = 1 (reverse order)
    jz = grd%imitt + 1
    do iz = grd%imitt + 1, grd%imitt + grd%ibl
      jz = jz - 1
      do kz = 1, grd%mxrho + grd%kbl
        flds%ebelam(kz, iz) = flds%ebelam(kz, jz)
        flds%ehbclam(kz, iz) = flds%ehbclam(kz, jz)
        flds%edu(kz, iz) = flds%edu(kz, jz)
      end do
    end do

    ! Second symmetry application (appears redundant but ensures consistency)
    jz = grd%imitt + 1
    do iz = grd%imitt + 1, grd%imitt + grd%ibl
      jz = jz - 1
      do kz = 1, grd%mxrho + grd%kbl
        flds%ebelam(kz, iz) = flds%ebelam(kz, jz)
        flds%ehbclam(kz, iz) = flds%ehbclam(kz, jz)
        flds%edu(kz, iz) = flds%edu(kz, jz)
      end do
    end do

  end subroutine apply_boundary_conditions

  ! ==========================================================================
  ! SUBROUTINE: propagate_polymer_chain
  ! ==========================================================================
  ! Propagates polymer chains segment by segment using Green's function approach
  !
  ! Calculates chain propagators c(r,i) for each segment i at each position r
  ! by integrating over bond orientations. Uses OpenMP parallelization for
  ! performance.
  !
  ! The propagation proceeds backward from end segments (imon = nmon-1 down to 1)
  ! ==========================================================================
  subroutine propagate_polymer_chain(inp, grd, comp, flds, cos_phi_table, c, cA, cB)
    use iso_fortran_env, only: real64, int32
    implicit none

    ! Arguments
    type(input_params_t), intent(in) :: inp
    type(grid_params_t), intent(in) :: grd
    type(computed_params_t), intent(in) :: comp
    type(fields_t), intent(in) :: flds
    real(real64), intent(in) :: cos_phi_table(:)
    real(real64), intent(inout) :: c(0:, 0:, :)  ! Chain propagators
    real(real64), intent(inout) :: cA(0:, 0:)    ! Forward propagator (working array)
    real(real64), intent(inout) :: cB(0:, 0:)    ! Backward propagator (working array)

    ! Local variables
    integer(int32) :: iz, kz, jz, iphi, irho, imon, kmon, jstart, irho0min
    real(real64) :: z, rho0, rho02, rt2, sume, zp, delz2, zpcsq, zpc2sq
    real(real64) :: phisum, zfact, rhoz2, fphi, rho2, rsq1, rsq2, valid
    real(real64) :: rho, fact, efact, ffact, bebbe, zpst, strho0

    ! Calculate cA: the propagator for polymer end segments
    ! cA(r) = exp(-beta*mu_end)*exp(-U_LJ) where:
    !   ebelam = exp(-beta*mu_end) from hard-sphere and chain connectivity
    !   edu = exp(-U_LJ) from Lennard-Jones interactions
    do iz = grd%istp1, grd%imitt + grd%ibl
    do kz = 1, grd%mxrho + grd%kbl
      cA(kz, iz) = flds%ebelam(kz, iz)*flds%edu(kz, iz)
    end do
    end do

    ! Propagate polymer chains segment by segment from end to end
    ! This calculates c(r,i) = propagator for segment i at position r
    imon = inp%nmon
    do kmon = 1, inp%nmon - 1
      imon = imon - 1

      ! Merge all chain propagation loops into single parallel region
      ! to eliminate thread synchronization barriers
!$omp parallel private(iz, z, jstart, zpst, irho0min, strho0, rho0, kz, rho02, rt2, sume, zp, jz, delz2, zpcsq, zpc2sq, phisum, zfact, rhoz2, fphi, iphi, rho2, rsq1, rsq2, valid, rho, irho, fact, efact, ffact, bebbe)

      ! Loop over all spatial grid points
!$omp do schedule(static)
      do iz = grd%istp1 + grd%ibl, grd%imitt
        z = inp%bl - 0.5d0*inp%dz + dble(iz - (grd%istp1 + grd%ibl) + 1)*inp%dz
        jstart = iz - grd%ibl
        zpst = z - inp%bl - inp%dz
        irho0min = 1
        strho0 = -0.5d0*inp%drho
        rho0 = strho0
        do kz = irho0min, grd%mxrho - grd%ibl
          rho0 = rho0 + inp%drho

          rho02 = rho0**2
          rt2 = rho02 + (z - inp%zc1)**2
          if (rt2 .lt. comp%Rcoll2) then
            c(kz, iz, imon) = 0.d0
            if (iz .gt. grd%imitt - grd%ibl - 1) cB(kz, grd%islut + 1 - iz) = 0.d0
            cB(kz, iz) = 0.d0
            cycle
          end if
          rt2 = rho02 + (z - comp%zc2)**2
          if (rt2 .lt. comp%Rcoll2) then
            c(kz, iz, imon) = 0.d0
            cB(kz, grd%islut + 1 - iz) = 0.d0
            cB(kz, iz) = 0.d0
            cycle
          end if

          ! Integrate over bond orientations: sum contributions from all points
          ! within bond length inp%bl of current position (rho0, z)
          sume = 0.d0
          zp = zpst
          do jz = jstart, iz + grd%ibl
            zp = zp + inp%dz
            delz2 = (zp - z)**2
            zpcsq = (zp - inp%zc1)**2
            zpc2sq = (zp - comp%zc2)**2
            phisum = 0.d0
            zfact = dabs(comp%bl2 - delz2)
            rhoz2 = rho0**2 + zfact
            fphi = 2.d0*rho0*dsqrt(zfact)
!$omp simd reduction(+:phisum)
            do iphi = 1, grd%nphi
!     Plus or minus sign doesn't matter for the value of the integral
              rho2 = rhoz2 - fphi*cos_phi_table(iphi)
              rsq1 = rho2 + zpcsq
              rsq2 = rho2 + zpc2sq

              ! Mask-based approach: 1.0 if outside both colloids, 0.0 if inside either
              ! This eliminates conditional exits (cycle) for full SIMD vectorization
              valid = merge(1.0d0, 0.0d0, rsq1 >= comp%Rcoll2 .and. rsq2 >= comp%Rcoll2)

              rho = dsqrt(rho2)
              irho = int(rho*comp%rdrho) + 1
              phisum = phisum + valid*cA(irho, jz)
            end do
!$omp end simd
            fact = 1.d0
            if (iabs(jz - iz) .eq. grd%ibl) fact = 0.5d0
            sume = 2.d0*phisum*inp%dphi*fact + sume
          end do
          efact = dsqrt(flds%edu(kz, iz))
          ffact = sume*comp%dzrfp*flds%ehbclam(kz, iz)/inp%bl*efact
          c(kz, iz, imon) = ffact
          if (iz .gt. grd%imitt - grd%ibl - 1) cB(kz, grd%islut + 1 - iz) = ffact*flds%ehbclam(kz, iz)*efact
          cB(kz, iz) = ffact*flds%ehbclam(kz, iz)*efact
        end do
      end do
!$omp end do
      ! Implicit barrier: Loop 3 reads cB(kz,iz) which Loop 1 writes

      ! Handle boundary regions: propagators at z-boundaries
!$omp do schedule(static)
      do iz = grd%istp1, grd%ibl
      do kz = 1, grd%mxrho + grd%kbl
        bebbe = comp%behbclam*cA(kz, iz)
        c(kz, iz, imon) = bebbe
        cA(kz, iz) = comp%behbclam*bebbe
      end do
      end do
!$omp end do nowait

      ! Handle radial boundaries and update propagators
!$omp do schedule(static)
      do iz = grd%ibl + 1, grd%imitt
        ! Outer radial boundary: use bulk propagators
        do kz = grd%mxrho - grd%kbl, grd%mxrho + grd%kbl
          bebbe = comp%behbclam*cA(kz, iz)
          c(kz, iz, imon) = bebbe
          cA(kz, iz) = comp%behbclam*bebbe
        end do
        ! Interior region: use backward propagator for next iteration
        do kz = 1, grd%mxrho - grd%ibl - 1
          cA(kz, iz) = cB(kz, iz)
        end do
      end do
!$omp end do
      ! Implicit barrier: Loop 4 reads cA which Loops 2 & 3 write

      ! Apply symmetry to propagators at midplane
!$omp do schedule(static)
      do iz = grd%imitt + 1, grd%imitt + grd%ibl
        jz = grd%imitt + 1 - (iz - grd%imitt)
        do kz = 1, grd%mxrho + grd%kbl
          cA(kz, iz) = cA(kz, jz)
        end do
      end do
!$omp end do

!$omp end parallel

    end do

  end subroutine propagate_polymer_chain

  ! ==========================================================================
  ! SUBROUTINE: calculate_adaptive_mixing
  ! ==========================================================================
  ! Calculates adaptive mixing parameters for SCF convergence
  !
  ! Implements smart mixing that:
  ! - Detects restart vs fresh start scenarios
  ! - Adjusts mixing aggressiveness based on convergence state
  ! - Detects and dampens oscillations
  !
  ! The mixing parameters control how much of the new vs old density is used
  ! ==========================================================================
  subroutine calculate_adaptive_mixing(niter, ddmax, inp, use_adaptive, &
                                       osc_count, ddmax_prev, dmm_out, dms_out)
    use iso_fortran_env, only: real64, int32
    implicit none

    ! Arguments
    integer(int32), intent(in) :: niter
    real(real64), intent(in) :: ddmax
    type(input_params_t), intent(in) :: inp
    logical, intent(inout) :: use_adaptive
    integer(int32), intent(inout) :: osc_count
    real(real64), intent(inout) :: ddmax_prev
    real(real64), intent(out) :: dmm_out, dms_out

    ! Smart adaptive mixing: detect restart scenarios and adjust accordingly
    ! Fresh start: ddmax ~100 at iteration 2 → use adaptive mixing
    ! Restart: ddmax ~0.005 at iteration 2 → use conservative mixing
    if (niter == 2 .and. ddmax < 0.1d0) then
      ! Restart scenario detected: ddmax already small at iteration 2
      use_adaptive = .false.
      write (*, *) 'Restart detected (ddmax < 0.1 at iter 2): using conservative mixing'
    end if

    if (.not. use_adaptive) then
      ! Conservative mixing for restart scenarios
      dmm_out = inp%dmm
      dms_out = inp%dms
    else
      ! Detect oscillations: ddmax increased after decreasing
      if (niter .gt. 10 .and. ddmax .gt. ddmax_prev .and. ddmax .lt. 0.3d0) then
        osc_count = osc_count + 1
      else if (ddmax .lt. ddmax_prev) then
        ! Reset oscillation counter when making progress
        osc_count = max(0, osc_count - 1)
      end if

      ! Adaptive mixing for fresh starts: adjust based on convergence state
      ! Add conservative bias when oscillations detected
      if (ddmax .gt. 1.0d0) then
        dmm_out = 0.90d0  ! Standard mixing when far from solution
        dms_out = 0.50d0
      else if (ddmax .gt. 0.1d0) then
        dmm_out = 0.85d0  ! Slightly more aggressive in mid-range
        dms_out = 0.45d0
      else if (ddmax .gt. 0.01d0) then
        ! When approaching convergence, use conservative mixing if oscillating
        if (osc_count .gt. 3) then
          dmm_out = 0.92d0  ! Very conservative to dampen oscillations
          dms_out = 0.52d0
        else
          dmm_out = 0.75d0  ! More aggressive when not oscillating
          dms_out = 0.35d0
        end if
      else if (ddmax .gt. 0.001d0) then
        if (osc_count .gt. 3) then
          dmm_out = 0.90d0  ! Very conservative near solution if oscillating
          dms_out = 0.50d0
        else
          dmm_out = 0.60d0  ! Aggressive near solution
          dms_out = 0.25d0
        end if
      else
        if (osc_count .gt. 3) then
          dmm_out = 0.88d0  ! Very conservative very close if oscillating
          dms_out = 0.48d0
        else
          dmm_out = 0.40d0  ! Extremely aggressive very close to solution
          dms_out = 0.15d0
        end if
      end if

      if (niter .gt. 60) then
        if (osc_count .gt. 3) then
          write (*, '(A,E12.5,A,F5.3,A,I3)') 'Adaptive mixing (osc): ddmax=', ddmax, ', dmm=', &
            dmm_out, ', osc_count=', osc_count
        else
          write (*, '(A,E12.5,A,F5.3)') 'Adaptive mixing: ddmax=', ddmax, ', dmm=', dmm_out
        end if
      end if

      ! Update previous ddmax for next iteration
      ddmax_prev = ddmax
    end if

  end subroutine calculate_adaptive_mixing

  ! ==========================================================================
  ! SUBROUTINE: update_densities_and_check_convergence
  ! ==========================================================================
  ! Updates density fields using mixing scheme and calculates convergence metric
  !
  ! Calculates new densities from chain propagators, mixes with old values,
  ! and computes the maximum relative change (ddmax) for convergence checking.
  ! Also applies symmetry boundary conditions at the midplane.
  ! ==========================================================================
  subroutine update_densities_and_check_convergence(inp, grd, comp, flds, c, &
                                                    dmm_mix, ddmax_out)
    use iso_fortran_env, only: real64, int32
    implicit none

    ! Arguments
    type(input_params_t), intent(in) :: inp
    type(grid_params_t), intent(in) :: grd
    type(computed_params_t), intent(in) :: comp
    type(fields_t), intent(inout) :: flds
    real(real64), intent(in) :: c(0:, 0:, :)
    real(real64), intent(in) :: dmm_mix
    real(real64), intent(out) :: ddmax_out

    ! Local variables
    integer(int32) :: i, j, k, iz, kz, jz
    real(real64) :: z, rho, rsq, diffz2, dumsum, tfem, tfdm, ddiff, tdmm

    tdmm = 1.d0 - dmm_mix

    ! Update densities using mixing scheme and check convergence
    ! Calculate new densities from propagators and mix with old values
    ddmax_out = 0.d0
    z = -0.5d0*inp%dz
    do i = grd%istp1, grd%imitt
      z = z + inp%dz
      diffz2 = (z - inp%zc1)**2
      rho = -0.5d0*inp%drho
      do j = 1, grd%mxrho
        rho = rho + inp%drho
        rsq = rho*rho + diffz2
        if (rsq .lt. comp%Rcoll2) then
          flds%fem(j, i) = 0.d0
          flds%fdmon(j, i) = 0.d0
        else
          ! Calculate total monomer density from chain propagators
          ! Sum over all internal segments (convolution of forward and backward propagators)
          dumsum = 0.d0
          do k = 2, inp%nmon - 1
            dumsum = c(j, i, k)*c(j, i, inp%nmon + 1 - k) + dumsum
          end do
          tfem = 2.d0*c(j, i, 1)*flds%ebelam(j, i)*dsqrt(flds%edu(j, i))/flds%ehbclam(j, i)
          tfdm = dumsum + tfem
          if (dabs(tfdm) .gt. 1.0d-14) then
            ddiff = abs(tfdm - flds%fdmon(j, i))/tfdm
            if (ddiff .gt. ddmax_out) ddmax_out = ddiff
          end if
          flds%fem(j, i) = flds%fem(j, i)*dmm_mix + tdmm*tfem
          flds%fdmon(j, i) = flds%fdmon(j, i)*dmm_mix + tdmm*tfdm
        end if
      end do
    end do

    ! Apply symmetry to updated densities at midplane
    ! Mirror fdmon and fem across z = imitt to maintain symmetry
    jz = grd%imitt + 1
    do iz = grd%imitt + 1, grd%imitt + grd%ibl
      jz = jz - 1
      do kz = 1, grd%mxrho + grd%kbl
        flds%fdmon(kz, iz) = flds%fdmon(kz, jz)
        flds%fem(kz, iz) = flds%fem(kz, jz)
      end do
    end do

  end subroutine update_densities_and_check_convergence

end module polymer_dft_data
