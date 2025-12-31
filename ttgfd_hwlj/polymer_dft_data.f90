! ============================================================================
! MODULE: polymer_dft_data
! ============================================================================
! Contains all shared data for the polymer DFT calculation.
! This module replaces the old COMMON blocks and include file approach.
!
! Migrated from: t2.inc.f90
! Date: 2025-12-31
! ============================================================================
module polymer_dft_data
  use iso_fortran_env, only: real64, int32
  implicit none

  ! ========================================================================
  ! Array dimensions - fixed size (can be made allocatable in future)
  ! ========================================================================
  integer, parameter :: maxel = 1001
  integer, parameter :: maxrho = 321
  integer, parameter :: maxphi = 5000  ! Maximum phi grid points for lookup tables

  ! ========================================================================
  ! Large arrays (formerly in COMMON/VECT/)
  ! ========================================================================
  ! Note: Using explicit double precision for COMMON block compatibility
  ! These will be converted to real(real64) when COMMON blocks are removed
  ! Note: SAVE attribute implicit in COMMON, so not specified
  double precision :: fdmon(0:maxrho, 0:maxel)
  double precision :: ebelam(0:maxrho, 0:maxel)
  double precision :: convp(0:maxrho, 0:maxel)
  double precision :: hvec(0:maxrho, 0:maxrho, 0:maxel)
  double precision :: fem(0:maxrho, 0:maxel)
  double precision :: ehbclam(0:maxrho, 0:maxel)
  double precision :: cdmonm(0:maxrho, 0:maxel)
  double precision :: ae1(0:maxrho, 0:maxel)
  double precision :: ae2(0:maxrho, 0:maxel)
  double precision :: edu(0:maxrho, 0:maxel)

  ! Cosine lookup tables for performance optimization
  double precision :: cos_phi(maxphi)
  double precision :: cos_pphi(maxphi)

  ! COMMON block for compatibility (will be removed in future)
  common /VECT/ fdmon, ebelam, convp, hvec, fem, ehbclam, cdmonm, ae1, ae2, edu, cos_phi, cos_pphi

  ! ========================================================================
  ! Grid and physical parameters (formerly in COMMON/VAR/)
  ! ========================================================================
  double precision :: dz, scalem
  double precision :: emscale
  double precision :: rdz
  double precision :: Yfact, rnmon, rrnmon
  double precision :: bdtot, chemps
  double precision :: drho, zc1, Rcoll2
  double precision :: Rcoll, behbclam, bebelam
  double precision :: zc2
  double precision :: dphi, dzrfp, rdrho, cdnorm
  double precision :: bdm, bl, dhs, dhs2, rdhs3, dhs3

  common /VAR/ dz, scalem, emscale, rdz, Yfact, rnmon, rrnmon, &
    bdtot, chemps, drho, zc1, Rcoll2, Rcoll, behbclam, bebelam, zc2, &
    dphi, dzrfp, rdrho, cdnorm, bdm, bl, dhs, dhs2, rdhs3, dhs3

  ! ========================================================================
  ! Integer grid dimensions (formerly in COMMON/HELTAL/)
  ! ========================================================================
  integer :: istart, istp1, islut, ism, nfack, imitt, nmon
  integer :: istp1s, isluts, mxrho, ksm, nphi
  integer :: ibl, kbl

  common /HELTAL/ istart, istp1, islut, ism, nfack, imitt, nmon, &
    istp1s, isluts, mxrho, ksm, nphi, ibl, kbl

  ! ========================================================================
  ! Physical and mathematical constants (compile-time parameters)
  ! ========================================================================

  ! Mathematical constants
  ! Note: Precision limited to double precision (15-17 significant digits)
  double precision, parameter :: PI = 3.14159265358979d0
  double precision, parameter :: TWOPI = 2.d0*PI
  double precision, parameter :: FOURPI = 4.d0*PI

  ! Carnahan-Starling equation of state parameters for hard spheres
  double precision, parameter :: A1_CS = 1.d0
  double precision, parameter :: A2_CS = 2.45696d0
  double precision, parameter :: B1_CS = 1.d0
  double precision, parameter :: B2_CS = 4.10386d0

  ! Convergence tolerance
  double precision, parameter :: CONV_TOL = 0.00001d0

  ! Derived constants (computed from fundamental constants)
  double precision, parameter :: PIS = PI/6.d0
  double precision, parameter :: C1 = -1.d0
  double precision, parameter :: C2 = -3.75503d0
  double precision, parameter :: AA1 = 2.d0*C1 - 2.d0*A1_CS - 4.d0
  double precision, parameter :: AA2 = 2.d0*C2 - 2.d0*A2_CS - 4.d0
  double precision, parameter :: BB1 = 3.d0 - B1_CS + A1_CS - 3.d0*C1
  double precision, parameter :: BB2 = 3.d0 - B2_CS + A2_CS - 3.d0*C2
  double precision, parameter :: Y = (9.82605d0 - 9.d0*PI*0.25d0)/(9.d0*PI*0.25d0 - 4.d0*PI/3.d0)

end module polymer_dft_data
