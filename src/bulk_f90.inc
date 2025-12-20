! =============================================================================
! Common blocks for Monte Carlo program
! =============================================================================

parameter (mion = 600, mxspec = 10)

character*16 title

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

double precision :: esa(mion, mion), ei(mion)
double precision :: caver(mxspec), hmik(mxspec, mxspec)
double precision :: cwi(25, mxspec, mxspec), dum(25)
double precision :: x6(mion), y6(mion), z6(mion), chv(mion), dp(mion)
double precision :: rw2(mion), rwi(mion), hc2v(mion, mxspec), hion(mxspec, 5)
integer          :: ispc(mion)
double precision :: pi, avno, bk, eps, epsx, ecf, ech, abeta, dtemp
double precision :: tx6, ty6, tz6, xww1
double precision :: pwiwal, pwimid, pwiwalv, pwimidv, pcollav, pcollv
double precision :: box, box2, box2i
integer          :: ny1, ny2, ny3, my3, nspec, ispec, npart, nwins, nwint
integer          :: nfix, il, isos, islu, ink
integer          :: lll, mmm, jjj, iii, kkk
double precision :: qfww1
