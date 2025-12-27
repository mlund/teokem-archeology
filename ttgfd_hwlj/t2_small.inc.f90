PARAMETER (maxel=201,maxrho=81)
COMMON/VECT/fdmon(0:maxrho,0:maxel),ebelam(0:maxrho,0:maxel), &
  convp(0:maxrho,0:maxel),hvec(0:maxel,0:maxrho,0:maxrho), &
  Vexm(0:maxrho,0:maxel),Virm(0:maxrho,0:maxel), &
  fem(0:maxrho,0:maxel), &
  ehbclam(0:maxrho,0:maxel),cdmonm(0:maxrho,0:maxel), &
  ae1(0:maxrho,0:maxel),ae2(0:maxrho,0:maxel), &
  Vexs(0:maxrho,0:maxel),Virs(0:maxrho,0:maxel), &
  econvp(0:maxrho,0:maxel),edu(0:maxrho,0:maxel)
COMMON/VAR/dz,closew,pis,pif,pit,vk,rrT,vkrrT,hvk,scalem,scales, &
  dzpie,AA1,AA2,BB1,BB2,c1,c2,Y,emscale,eblemb,ehblcmb,eblsmb,bcdt, &
  rrjkdiff,threqz,rtwelve,pie,rthree,rdz,btrams, &
  sclosew,q1,q2,q3,p1,p2,p3,r2,r1,r0,s2,s1,s0,b2,b1,b0,r2sq,r1sq, &
  r0sq,Yfact,veq,rnmon,rrnmon,rrcmon,rq3,cdmbulk,cdsbulk, &
  cdmlbulk,cdslbulk,elblemb,elhblcmb,elblsmb,distp1,dmitt,eplag, &
  seplag,bdtot,rrnarm,rnarm,tscalem,donn,csurf,chemps, &
  drho,zc1,Rcoll2, &
  dr,rdr,Rcoll,dct,behbclam,bebelam, &
  Rcollref,zc2,Rtr,utr,Rtr2,amp,Rpot, &
  dphi,dzrfp,rdrho,cdnorm, &
  eRtr,eutr,eRtr2,eamp, &
  rbl,Ascmm,tkappa,bdm,Ascpm,bl,admhs,dhs,dhs2,rdhs3,dhs3
COMMON/HELTAL/istart,istp1,islut,ism,inw,nfack,imitt,nmon, &
  ist,ifin,istp1s,isluts,isms,inws,kst,kfin,mxrho,ksm,nct,nphi, &
  ibl,kbl



