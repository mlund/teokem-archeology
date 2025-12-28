! Array dimensions - fixed size for COMMON block compatibility
INTEGER, PARAMETER :: maxel = 1001
INTEGER, PARAMETER :: maxrho = 321
INTEGER, PARAMETER :: maxphi = 5000  ! Maximum phi grid points for lookup tables

! Large arrays in COMMON block (subroutines need access via COMMON)
DOUBLE PRECISION :: fdmon(0:maxrho, 0:maxel), ebelam(0:maxrho, 0:maxel)
DOUBLE PRECISION :: convp(0:maxrho, 0:maxel), hvec(0:maxel, 0:maxrho, 0:maxrho)
DOUBLE PRECISION :: fem(0:maxrho, 0:maxel), ehbclam(0:maxrho, 0:maxel)
DOUBLE PRECISION :: cdmonm(0:maxrho, 0:maxel)
DOUBLE PRECISION :: ae1(0:maxrho, 0:maxel), ae2(0:maxrho, 0:maxel)
DOUBLE PRECISION :: edu(0:maxrho, 0:maxel)

! Cosine lookup tables for performance optimization (precomputed trigonometric values)
DOUBLE PRECISION :: cos_phi(maxphi), cos_pphi(maxphi)

COMMON/VECT/fdmon, ebelam, convp, hvec, fem, ehbclam, cdmonm, ae1, ae2, edu, cos_phi, cos_pphi

! Type declarations for COMMON block variables (required for implicit none)
! COMMON/VAR/ - double precision variables
DOUBLE PRECISION :: dz, scalem
DOUBLE PRECISION :: emscale
DOUBLE PRECISION :: rdz
DOUBLE PRECISION :: Yfact, rnmon, rrnmon
DOUBLE PRECISION :: bdtot, chemps
DOUBLE PRECISION :: drho, zc1, Rcoll2
DOUBLE PRECISION :: Rcoll, behbclam, bebelam
DOUBLE PRECISION :: zc2
DOUBLE PRECISION :: dphi, dzrfp, rdrho, cdnorm
DOUBLE PRECISION :: bdm, bl, dhs, dhs2, rdhs3, dhs3

! COMMON/HELTAL/ - integer variables
INTEGER :: istart, istp1, islut, ism, nfack, imitt, nmon
INTEGER :: istp1s, isluts, mxrho, ksm, nphi
INTEGER :: ibl, kbl

COMMON/VAR/dz, scalem, emscale, rdz, Yfact, rnmon, rrnmon, &
  bdtot, chemps, drho, zc1, Rcoll2, Rcoll, behbclam, bebelam, zc2, &
  dphi, dzrfp, rdrho, cdnorm, bdm, bl, dhs, dhs2, rdhs3, dhs3
COMMON/HELTAL/istart, istp1, islut, ism, nfack, imitt, nmon, &
  istp1s, isluts, mxrho, ksm, nphi, ibl, kbl

! ========================================================================
! Physical and mathematical constants (compile-time parameters)
! ========================================================================

! Mathematical constants
DOUBLE PRECISION, PARAMETER :: PI = 3.141592653589793238462643383279502884197d0
DOUBLE PRECISION, PARAMETER :: TWOPI = 2.d0*PI
DOUBLE PRECISION, PARAMETER :: FOURPI = 4.d0*PI
DOUBLE PRECISION, PARAMETER :: VOLFACT = FOURPI/3.d0
DOUBLE PRECISION, PARAMETER :: RVOLFACT = 1.d0/VOLFACT

! Physical constants
DOUBLE PRECISION, PARAMETER :: BK = 1.38066D-23        ! Boltzmann constant (J/K)
DOUBLE PRECISION, PARAMETER :: AVNO = 6.02214D23       ! Avogadro's number (1/mol)
DOUBLE PRECISION, PARAMETER :: ELCH = 1.602D-19        ! Elementary charge (C)
DOUBLE PRECISION, PARAMETER :: VACUUM_PERMITTIVITY = 8.85418782D-12  ! F/m
DOUBLE PRECISION, PARAMETER :: DIELC_WATER = 78.3d0    ! Dielectric constant (water)

! Carnahan-Starling equation of state parameters for hard spheres
DOUBLE PRECISION, PARAMETER :: A1_CS = 1.d0
DOUBLE PRECISION, PARAMETER :: A2_CS = 2.45696d0
DOUBLE PRECISION, PARAMETER :: B1_CS = 1.d0
DOUBLE PRECISION, PARAMETER :: B2_CS = 4.10386d0

! Convergence tolerance
DOUBLE PRECISION, PARAMETER :: CONV_TOL = 0.00001d0

! Derived constants (computed from fundamental constants)
DOUBLE PRECISION, PARAMETER :: PIS = PI/6.d0
DOUBLE PRECISION, PARAMETER :: PIT = PI/3.d0
DOUBLE PRECISION, PARAMETER :: PIF = PI/4.d0
DOUBLE PRECISION, PARAMETER :: C1 = -1.d0
DOUBLE PRECISION, PARAMETER :: C2 = -3.75503d0
DOUBLE PRECISION, PARAMETER :: AA1 = 2.d0*C1 - 2.d0*A1_CS - 4.d0
DOUBLE PRECISION, PARAMETER :: AA2 = 2.d0*C2 - 2.d0*A2_CS - 4.d0
DOUBLE PRECISION, PARAMETER :: BB1 = 3.d0 - B1_CS + A1_CS - 3.d0*C1
DOUBLE PRECISION, PARAMETER :: BB2 = 3.d0 - B2_CS + A2_CS - 3.d0*C2
DOUBLE PRECISION, PARAMETER :: Y = (9.82605d0 - 9.d0*PI*0.25d0)/(9.d0*PI*0.25d0 - 4.d0*PI/3.d0)
