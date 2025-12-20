! =============================================================================
! Common blocks for Monte Carlo program
! All variables are explicitly typed (no reliance on implicit typing rules)
! =============================================================================

! Parameters
integer, parameter :: mion   = 600
integer, parameter :: mxspec = 10

! Character variable
character*16 :: title

! Array variables
double precision :: esa(mion, mion)
double precision :: ei(mion)
double precision :: caver(mxspec)
double precision :: hmik(mxspec, mxspec)
double precision :: cwi(25, mxspec, mxspec)
double precision :: dum(25)
double precision :: x6(mion)
double precision :: y6(mion)
double precision :: z6(mion)
double precision :: chv(mion)
double precision :: dp(mion)
double precision :: rw2(mion)
double precision :: rwi(mion)
double precision :: hc2v(mion, mxspec)
double precision :: hion(mxspec, 5)

! Integer arrays
integer :: ispc(mion)

! Scalar double precision variables
double precision :: pi
double precision :: avno
double precision :: bk
double precision :: eps
double precision :: epsx
double precision :: ecf
double precision :: ech
double precision :: abeta
double precision :: dtemp
double precision :: tx6
double precision :: ty6
double precision :: tz6
double precision :: xww1
double precision :: pwiwal
double precision :: pwimid
double precision :: pwiwalv
double precision :: pwimidv
double precision :: pcollav
double precision :: pcollv
double precision :: box
double precision :: box2
double precision :: box2i
double precision :: qfww1

! Integer scalar variables
integer :: ny1
integer :: ny2
integer :: ny3
integer :: my3
integer :: nspec
integer :: ispec
integer :: npart
integer :: nwins
integer :: nwint
integer :: nfix
integer :: il
integer :: isos
integer :: islu
integer :: ink
integer :: lll
integer :: mmm
integer :: jjj
integer :: iii
integer :: kkk

! Common block declaration
common /arr1/ esa, ei, &
              caver, hmik, &
              cwi, dum, &
              x6, y6, z6, chv, dp, &
              rw2, rwi, hc2v, hion, ispc, &
              pi, avno, bk, eps, epsx, ecf, ech, abeta, dtemp, &
              tx6, ty6, tz6, xww1, &
              pwiwal, pwimid, pwiwalv, pwimidv, pcollav, pcollv, &
              box, box2, box2i, &
              ny1, ny2, ny3, my3, nspec, ispec, npart, nwins, nwint, &
              nfix, il, isos, islu, ink, &
              lll, mmm, jjj, iii, kkk, &
              qfww1
