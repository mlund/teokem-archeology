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
DOUBLE PRECISION :: dz, closew, PIS, PIF, PIT, vk, rrT, vkrrT, hvk, scalem, scales
DOUBLE PRECISION :: dzpie, AA1, AA2, BB1, BB2, C1, C2, Y, emscale, eblemb, ehblcmb, eblsmb, bcdt
DOUBLE PRECISION :: rrjkdiff, threqz, rtwelve, pie, rthree, rdz, btrams
DOUBLE PRECISION :: sclosew, q1, q2, q3, p1, p2, p3, r2, r1, r0, s2, s1, s0, b2, b1, b0, r2sq, r1sq
DOUBLE PRECISION :: r0sq, Yfact, veq, rnmon, rrnmon, rrcmon, rq3, cdmbulk, cdsbulk
DOUBLE PRECISION :: cdmlbulk, cdslbulk, elblemb, elhblcmb, elblsmb, distp1, dmitt, eplag
DOUBLE PRECISION :: seplag, bdtot, rrnarm, rnarm, tscalem, donn, csurf, chemps
DOUBLE PRECISION :: drho, zc1, Rcoll2
DOUBLE PRECISION :: dr, rdr, Rcoll, dct, behbclam, bebelam
DOUBLE PRECISION :: Rcollref, zc2, Rtr, utr, Rtr2, amp, Rpot
DOUBLE PRECISION :: dphi, dzrfp, rdrho, cdnorm
DOUBLE PRECISION :: eRtr, eutr, eRtr2, eamp
DOUBLE PRECISION :: rbl, Ascmm, tkappa, bdm, Ascpm, bl, admhs, dhs, dhs2, rdhs3, dhs3

! COMMON/HELTAL/ - integer variables
INTEGER :: istart, istp1, islut, ism, inw, nfack, imitt, nmon
INTEGER :: ist, ifin, istp1s, isluts, isms, inws, kst, kfin, mxrho, ksm, nct, nphi
INTEGER :: ibl, kbl

COMMON/VAR/dz, closew, PIS, PIF, PIT, vk, rrT, vkrrT, hvk, scalem, scales, &
  dzpie, AA1, AA2, BB1, BB2, C1, C2, Y, emscale, eblemb, ehblcmb, eblsmb, bcdt, &
  rrjkdiff, threqz, rtwelve, pie, rthree, rdz, btrams, &
  sclosew, q1, q2, q3, p1, p2, p3, r2, r1, r0, s2, s1, s0, b2, b1, b0, r2sq, r1sq, &
  r0sq, Yfact, veq, rnmon, rrnmon, rrcmon, rq3, cdmbulk, cdsbulk, &
  cdmlbulk, cdslbulk, elblemb, elhblcmb, elblsmb, distp1, dmitt, eplag, &
  seplag, bdtot, rrnarm, rnarm, tscalem, donn, csurf, chemps, &
  drho, zc1, Rcoll2, &
  dr, rdr, Rcoll, dct, behbclam, bebelam, &
  Rcollref, zc2, Rtr, utr, Rtr2, amp, Rpot, &
  dphi, dzrfp, rdrho, cdnorm, &
  eRtr, eutr, eRtr2, eamp, &
  rbl, Ascmm, tkappa, bdm, Ascpm, bl, admhs, dhs, dhs2, rdhs3, dhs3
COMMON/HELTAL/istart, istp1, islut, ism, inw, nfack, imitt, nmon, &
  ist, ifin, istp1s, isluts, isms, inws, kst, kfin, mxrho, ksm, nct, nphi, &
  ibl, kbl

