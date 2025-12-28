program platem
  implicit double precision(a - h, o - z)
  parameter(maxmon=2401)
  include 't2.inc.f90'
  dimension c(0:maxrho, 0:maxel, maxmon), &
    cA(0:maxrho, 0:maxel), cB(0:maxrho, 0:maxel), &
    cdens(0:1000), ctvec(0:1000)
  ifc = 38
  ins = 49
  iep = 50
  pi = acos(-1.d0)
  bk = 1.38066D-23
  avno = 6.02214D23
  elch = 1.602D-19
  faraday = 8.85418782D-12
  dielc = 78.3d0
  fourpi = 4.d0*pi
  twopi = 2.d0*pi
  fourpi = 4.d0*pi
  rtwelve = 1.d0/12.d0
  rthree = 1.d0/3.d0
  volfact = fourpi/3.d0
  rnine = 1.d0/9.d0
  rthree = 1.d0/3.d0
  rphi = fourpi*0.2d0
  aphi = fourpi*0.5d0
  rvolfact = 1.d0/volfact
  a1 = 1.d0
  a2 = 2.45696d0
  b1 = 1.d0
  b2 = 4.10386d0
  c1 = -1.d0
  c2 = -3.75503d0
  AA1 = 2.d0*c1 - 2.d0*a1 - 4.d0
  AA2 = 2.d0*c2 - 2.d0*a2 - 4.d0
  BB1 = 3.d0 - b1 + a1 - 3.d0*c1
  BB2 = 3.d0 - b2 + a2 - 3.d0*c2
  Y = (9.82605d0 - 9.d0*pi*0.25d0)/(9.d0*pi*0.25d0 - 4.d0*pi/3.d0)
  pis = pi/6.d0
  pit = pi/3.d0
  pif = pi/4.d0
  ddtol = 0.00001d0
  open (ifc, file='fcdfil', form='formatted')
  open (ins, file='input.tsph', form='formatted')
  open (iep, file='epfil', form='formatted')
  rewind ifc
  rewind ins
  rewind iep
  read (ins, *) bdm
  read (ins, *) bdtot
  bds = bdtot - bdm
  read (ins, *) nmon
  read (ins, *) dz
  read (ins, *) drho
  read (ins, *) dphi
  dphi = pi*dphi
  read (ins, *) Rcoll
  read (ins, *) zc1
  read (ins, *) collsep
  read (ins, *) Rcyl
  read (ins, *) ioimaxm
!      ioimaxm = 10
  read (ins, *) dmm, dms
  read (ins, *) kread
!      kread = 0
  read (ins, *) bl
  read (ins, *) dhs
  read (ins, *) dpphi
  rbl = 1.d0/bl
  bl2 = bl*bl
  dhs2 = dhs*dhs
  dhs3 = dhs2*dhs
  rdhs3 = 1.d0/dhs3

  read (iep, *) epslj
  alj = 4.d0*epslj*dhs3*dhs3
  rlj = 4.d0*epslj*dhs3**4

  Rcyl2 = Rcyl*Rcyl
  tdmm = 1.d0 - dmm
  tdms = 1.d0 - dms
  rdz = 1.d0/dz
  rdrho = 1.d0/drho
  rdphi = 1.d0/dphi
  twopidz = twopi*dz
  dzrfp = dz/(4.d0*pi)
!      dzrfp = dz/(4.d0*pi*dhs**3)
  nphi = int(pi/dphi + 0.01d0)
  irdz = int(rdz + 0.001d0)
  nfack = int(2.d0*(zc1 + 0.5d0*collsep)/dz + 0.01d0)
  istart = 0
  istp1 = 1
  islut = nfack
  istp1s = 1
  isluts = islut
!      ism = int(1.d0/dz+0.01d0)
!      ksm = int(1.d0/drho+0.01d0)
  ism = int(dhs/dz + 0.01d0)
  ksm = int(dhs/drho + 0.01d0)
  ibl = int(bl/dz + 0.01d0)
  kbl = int(bl/drho + 0.01d0)
  pie = pi/8.d0
  dzpie = pie*dz
  rnmon = dble(nmon)
  rrnmon = 1.d0/rnmon
  rrcmon = 1.d0/(rnmon - 2.d0)

  Yfact = (rnmon - 2.d0)*Y
  bdpol = bdm/rnmon

  bdt = bdm*dhs3
  aeta = pis*bdt
  xsib = 1.d0 - aeta
  rxsib = 1.d0/xsib
  rxsibsq = rxsib*rxsib

  aex1 = -(c1 + 1.d0)*dlog(xsib) - &
         0.5d0*(AA1*pis*bdt + BB1*(pis*bdt)**2)*rxsibsq
  aex2 = -(c2 + 1.d0)*dlog(xsib) - &
         0.5d0*(AA2*pis*bdt + BB2*(pis*bdt)**2)*rxsibsq
  bFex = (bdm - 2.d0*bdpol)*Y*(aex2 - aex1) + bdpol*aex2

  daex1 = rxsib*(c1 + 1 - 0.5d0*(AA1 + 2.d0*BB1*aeta)*rxsib - &
                 aeta*(AA1 + BB1*aeta)*rxsibsq)
  daex2 = rxsib*(c2 + 1 - 0.5d0*(AA2 + 2.d0*BB2*aeta)*rxsib - &
                 aeta*(AA2 + BB2*aeta)*rxsibsq)
  pdasum = Yfact*(daex2 - daex1) + daex2
  baex1 = aex1
  baex2 = aex2
  bdaex1 = daex1
  bdaex2 = daex2

  Pb = bdpol + bdpol*aeta*pdasum
  chempp = dlog(bdpol) + Yfact*(aex2 - aex1) + aex2 + pis*bdm*pdasum*dhs3
  scalem = chempp/(2.d0*rnmon)
  emscale = 2.d0*scalem
  bconvp = (Y*(bdm - 2.d0*bdm*rrnmon)*(daex2 - daex1) + &
            bdm*rrnmon*daex2)*pis*dhs3
  trams = bconvp
  emtrams = trams + 0.5d0*aex2
  cmtrams = trams + Y*(aex2 - aex1)
  bemtrams = emtrams
  bcmtrams = cmtrams
  bebelam = dexp(-emtrams + emscale)
  behbclam = dexp(-0.5d0*cmtrams + scalem)

  imitt = nfack/2
  write (*, *) 'GFD POLYMER SOLUTION MODEL!'
  write (*, *) 'bdm,bdpol =', bdm, bdpol
  write (*, *) 'monomer density  = ', bdm
  write (*, *) 'bdt = ', bdt
  write (*, *) 'collsep,dz = ', collsep, dz
  write (*, *) 'Rcoll = ', Rcoll
  write (*, *) 'bond length (bl): ', bl
  write (*, *) 'monomer hs diameter (bl): ', dhs
  write (*, *) 'no. of monomers/polymer = ', nmon
  write (*, *) 'max no. of iterations = ', ioimaxm
  write (*, *) 'polymer chemical pot. (betamu) = ', chempp
  write (*, *) 'solvent chemical pot. (betamu) = ', chemps
  write (*, *) 'total bulk pressure = ', Pb
  zc2 = zc1 + collsep
  write (*, *) 'zc1,zc2 = ', zc1, zc2
  write (*, *) 'nfack,imitt = ', nfack, imitt
  write (*, *) 'istp1,islut = ', istp1, islut
  write (*, *) 'istp1s,isluts = ', istp1s, isluts
  write (*, *) 'ism,ibl = ', ism, ibl
  write (*, *) 'ksm,kbl = ', ksm, kbl
  write (*, *) 'dmm,dms (density mixing param. mon.,solv.) = ', dmm, dms
  write (*, *) 'Rcyl  = ', Rcyl
  write (*, *) 'bFex = ', bFex
  write (*, *) 'bebelam,behbclam = ', bebelam, behbclam
  Rcoll2 = Rcoll*Rcoll
  mxrho = int((Rcyl - 1.d0)*rdrho) + 1

  CALL CDFACT
  write (*, *) 'cdnorm = ', cdnorm
  do iz = istp1, istp1 + 2*ibl - 1
  do kz = 1, mxrho + kbl
    cdt = bdm*dhs3
    pcdt = pis*cdt
    xsi = (1.d0 - pcdt)
    rxsi = 1.d0/xsi
    sqrxsi = rxsi*rxsi
    flog = dlog(xsi)
    ae1(kz, iz) = -(c1 + 1.d0)*flog - 0.5d0*(AA1 + BB1*pcdt)*pcdt*sqrxsi
    ae2(kz, iz) = -(c2 + 1.d0)*flog - 0.5d0*(AA2 + BB2*pcdt)*pcdt*sqrxsi
    daex1 = rxsi*(c1 + 1.d0 - 0.5d0*AA1*rxsi*(1.d0 + 2.d0*pcdt*rxsi) - &
                  BB1*pcdt*rxsi*(1.d0 + pcdt*rxsi))
    daex2 = rxsi*(c2 + 1.d0 - 0.5d0*AA2*rxsi*(1.d0 + 2.d0*pcdt*rxsi) - &
                  BB2*pcdt*rxsi*(1.d0 + pcdt*rxsi))
    convp(kz, iz) = (Y*(bdm - 2.d0*bdm*rrnmon)*(daex2 - daex1) + &
                     bdm*rrnmon*daex2)*pis*dhs3
!     *bdm*rrnmon*daex2)*pis
  end do
  end do
  do iz = istp1 + 2*ibl, imitt + ibl
  do kz = mxrho - kbl + 1, mxrho + kbl
    cdt = bdm*dhs3
    pcdt = pis*cdt
    xsi = (1.d0 - pcdt)
    rxsi = 1.d0/xsi
    sqrxsi = rxsi*rxsi
    flog = dlog(xsi)
    ae1(kz, iz) = -(c1 + 1.d0)*flog - 0.5d0*(AA1 + BB1*pcdt)*pcdt*sqrxsi
    ae2(kz, iz) = -(c2 + 1.d0)*flog - 0.5d0*(AA2 + BB2*pcdt)*pcdt*sqrxsi
    daex1 = rxsi*(c1 + 1.d0 - 0.5d0*AA1*rxsi*(1.d0 + 2.d0*pcdt*rxsi) - &
                  BB1*pcdt*rxsi*(1.d0 + pcdt*rxsi))
    daex2 = rxsi*(c2 + 1.d0 - 0.5d0*AA2*rxsi*(1.d0 + 2.d0*pcdt*rxsi) - &
                  BB2*pcdt*rxsi*(1.d0 + pcdt*rxsi))
    convp(kz, iz) = (Y*(bdm - 2.d0*bdm*rrnmon)*(daex2 - daex1) + &
                     bdm*rrnmon*daex2)*pis*dhs3
!     *bdm*rrnmon*daex2)*pis
  end do
  end do

  if (kread .eq. 0) then
    z = -0.5d0*dz
    do iz = 1, ibl
      z = z + dz
    end do
    do iz = ibl + 1, imitt
      z = z + dz
      z2 = (z - zc1)**2
      z22 = (z - zc2)**2
      rho = -0.5d0*drho
      do kz = 1, mxrho
        rho = rho + drho
        fdmon(kz, iz) = bdm
        fem(kz, iz) = 2.d0*fdmon(kz, iz)*rrnmon
        ebelam(kz, iz) = bebelam
        ehbclam(kz, iz) = behbclam
        rt2 = rho*rho + z2
        if (rt2 .lt. Rcoll2) then
          fdmon(kz, iz) = 0.d0
          fem(kz, iz) = 0.d0
          ebelam(kz, iz) = 0.d0
          ehbclam(kz, iz) = 0.d0
        end if
        rt2 = rho*rho + z22
        if (rt2 .lt. Rcoll2) then
          fdmon(kz, iz) = 0.d0
          fem(kz, iz) = 0.d0
          ebelam(kz, iz) = 0.d0
          ehbclam(kz, iz) = 0.d0
        end if
      end do
    end do

  else
    rewind ifc
    do iz = istp1, imitt
    do kz = 1, mxrho
      read (ifc, *) t1, t2, fdmon(kz, iz), trams
      fem(kz, iz) = 2.d0*fdmon(kz, iz)*rrnmon
    end do
    end do
  end if

  do iz = istp1, ibl
  do kz = 1, mxrho + kbl
    fdmon(kz, iz) = bdm
    fem(kz, iz) = 2.d0*fdmon(kz, iz)*rrnmon
    ebelam(kz, iz) = bebelam
    ehbclam(kz, iz) = behbclam
    cdmonm(kz, iz) = bdm
  end do
  end do
  do iz = 1, imitt
  do kz = mxrho + 1, mxrho + kbl
    fdmon(kz, iz) = bdm
    fem(kz, iz) = 2.d0*fdmon(kz, iz)*rrnmon
    ebelam(kz, iz) = bebelam
    ehbclam(kz, iz) = behbclam
    cdmonm(kz, iz) = bdm
  end do
  end do
  jz = imitt + 1
  do iz = imitt + 1, imitt + ibl
    jz = jz - 1
    do kz = 1, mxrho + kbl
      fdmon(kz, iz) = fdmon(kz, jz)
      fem(kz, iz) = fem(kz, jz)
      ebelam(kz, iz) = ebelam(kz, jz)
      ehbclam(kz, iz) = ehbclam(kz, jz)
      cdmonm(kz, iz) = bdm
    end do
  end do
    write (*, *) 'fdmon(1,1) = ', fdmon(1, 1)
    write (*, *) 'fdmon(1,11) = ', fdmon(1, 11)

    write (*, *) 'dpphi = ', dpphi
    dpphi = dpphi*pi
    npphi = int(pi/dpphi + 0.01d0)
    write (*, *) 'dpphi,npphi = ', dpphi, npphi

    kcm = nfack
    tdz = -dz
    do itdz = 0, kcm - 1
!     tdz loop
      tdz = tdz + dz
      tdzsq = tdz*tdz
      rho = -0.5d0*drho
      do krho = 1, mxrho
!     rho loop
        rho = rho + drho
        rhosq = rho*rho
        use1 = tdzsq + rhosq
        trho = -0.5d0*drho
        do kprho = 1, mxrho
!     rhop loop
          trho = trho + drho
          trhosq = trho*trho
          trmix = 2.d0*trho*rho
          useful = use1 + trhosq
          pint = 0.d0
          phi = -0.5d0*dpphi
          do iphi = 1, npphi
            phi = phi + dpphi
!      s = dsqrt(useful-trmix*dcos(phi))
            s2 = useful - trmix*dcos(phi)
            if (s2 .gt. dhs2) then
!      s = dsqrt(s2)
!      pint = dexp(-tkappa*s)/s+pint
              pint = rlj/s2**6 - alj/s2**3 + pint
!      uuu = rlj/s2**6-alj/s2**3
!      if (uuu.lt.50.d0) pint = uuu+pint
!      if (uuu.lt.16.d0) pint = 0.01d0+pint
!      pint = 0.01d0+pint
!      pint = uuu+pint
            end if
          end do
          hvec(itdz, krho, kprho) = trho*pint*dpphi
        end do
      end do
    end do
    write (*, *) 'hvec fixad'

    ddmax = 10000.d0
    niter = 0

    ! Main iteration loop
    do while (.true.)
      niter = niter + 1
      write (*, *) 'ddmax,niter = ', ddmax, niter
      write (*, *)
      if (niter .gt. ioimaxm) then
        write (*, *) 'NITER.GT.IOIMAXM !', niter
        exit
      end if

      CALL CDCALC
    CALL AVEC
    CALL EBLMNEW
    CALL EBDU

    do iz = istp1, imitt
    do kz = mxrho + 1, mxrho + kbl
      ebelam(kz, iz) = bebelam
      ehbclam(kz, iz) = behbclam
      edu(kz, iz) = 1.d0
    end do
    end do
    jz = imitt + 1
    do iz = imitt + 1, imitt + ibl
      jz = jz - 1
      do kz = 1, mxrho + kbl
        ebelam(kz, iz) = ebelam(kz, jz)
        ehbclam(kz, iz) = ehbclam(kz, jz)
        edu(kz, iz) = edu(kz, jz)
      end do
    end do
    jz = imitt + 1
    do iz = imitt + 1, imitt + ibl
      jz = jz - 1
      do kz = 1, mxrho + kbl
        ebelam(kz, iz) = ebelam(kz, jz)
        ehbclam(kz, iz) = ehbclam(kz, jz)
        edu(kz, iz) = edu(kz, jz)
      end do
    end do
    do iz = istp1, imitt + ibl
    do kz = 1, mxrho + kbl
      cA(kz, iz) = ebelam(kz, iz)*edu(kz, iz)
    end do
    end do

    imon = nmon
    do kmon = 1, nmon - 1
      imon = imon - 1

      z = bl - 0.5d0*dz
      do iz = istp1 + ibl, imitt
        z = z + dz
        jstart = iz - ibl
        zpst = z - bl - dz
        irho0min = 1
        strho0 = -0.5d0*drho
        rho0 = strho0
        do kz = irho0min, mxrho - ibl
          rho0 = rho0 + drho

          rho02 = rho0**2
          rt2 = rho02 + (z - zc1)**2
          if (rt2 .lt. Rcoll2) then
            c(kz, iz, imon) = 0.d0
            if (iz .gt. imitt - ibl - 1) cB(kz, islut + 1 - iz) = 0.d0
            cB(kz, iz) = 0.d0
            cycle
          end if
          rt2 = rho02 + (z - zc2)**2
          if (rt2 .lt. Rcoll2) then
            c(kz, iz, imon) = 0.d0
            cB(kz, islut + 1 - iz) = 0.d0
            cB(kz, iz) = 0.d0
            cycle
          end if

          sume = 0.d0
          zp = zpst
          do jz = jstart, iz + ibl
            zp = zp + dz
            delz2 = (zp - z)**2
            zpcsq = (zp - zc1)**2
            zpc2sq = (zp - zc2)**2
            phisum = 0.d0
            phi = -0.5d0*dphi
            zfact = dabs(bl2 - delz2)
            rhoz2 = rho0**2 + zfact
            fphi = 2.d0*rho0*dsqrt(zfact)
            do iphi = 1, nphi
              phi = phi + dphi
!     plus eller minus spelar ingen roll foer integralens vaerde
              rho2 = rhoz2 - fphi*dcos(phi)
              rsq = rho2 + zpcsq
              if (rsq .lt. Rcoll2) cycle
              rsq = rho2 + zpc2sq
              if (rsq .lt. Rcoll2) cycle
              rho = dsqrt(rho2)
              irho = int(rho*rdrho) + 1
              phisum = cA(irho, jz) + phisum
            end do
            fact = 1.d0
            if (iabs(jz - iz) .eq. ibl) fact = 0.5d0
            sume = 2.d0*phisum*dphi*fact + sume
          end do
          efact = dsqrt(edu(kz, iz))
          ffact = sume*dzrfp*ehbclam(kz, iz)*rbl*efact
          c(kz, iz, imon) = ffact
          if (iz .gt. imitt - ibl - 1) cB(kz, islut + 1 - iz) = ffact*ehbclam(kz, iz)*efact
          cB(kz, iz) = ffact*ehbclam(kz, iz)*efact
!     rho0 svepet faerdigt!
        end do
!     dags foer nytt z0 vaerde!
      end do
!     dags foer ny monomer laengs kedjan
        do iz = istp1, ibl
        do kz = 1, mxrho + kbl
          bebbe = behbclam*cA(kz, iz)
          c(kz, iz, imon) = bebbe
          cA(kz, iz) = behbclam*bebbe
        end do
        end do
        do iz = ibl + 1, imitt
        do kz = mxrho - kbl, mxrho + kbl
          bebbe = behbclam*cA(kz, iz)
          c(kz, iz, imon) = bebbe
          cA(kz, iz) = behbclam*bebbe
        end do
        do kz = 1, mxrho - ibl - 1
          cA(kz, iz) = cB(kz, iz)
        end do
        end do

        jz = imitt + 1
        do iz = imitt + 1, imitt + ibl
          jz = jz - 1
          do kz = 1, mxrho + kbl
            cA(kz, iz) = cA(kz, jz)
          end do
        end do

      end do

      if (ddmax .lt. ddtol) exit  ! Converged

      ddmax = 0.d0
      z = -0.5d0*dz
      do i = istp1, imitt
        z = z + dz
        diffz2 = (z - zc1)**2
        rho = -0.5d0*drho
        do j = 1, mxrho
          rho = rho + drho
          rsq = rho*rho + diffz2
          if (rsq .lt. Rcoll2) then
            fem(j, i) = 0.d0
            fdmon(j, i) = 0.d0
          else
            dumsum = 0.d0
            do k = 2, nmon - 1
              dumsum = c(j, i, k)*c(j, i, nmon + 1 - k) + dumsum
            end do
!      tfem = 2.d0*c(j,i,1)*ehbclam(j,i)
            tfem = 2.d0*c(j, i, 1)*ebelam(j, i)*dsqrt(edu(j, i))/ehbclam(j, i)
            tfdm = dumsum + tfem
            if (dabs(tfdm) .gt. 1.0d-14) then
              ddiff = abs(tfdm - fdmon(j, i))/tfdm
              if (ddiff .gt. ddmax) ddmax = ddiff
            end if
!      if (tfdm.gt.(10.*fdm)) then
!      tfdm = 10.d0*fdm
!      tfem = 10.d0*fem(j,i)
!      endif
            fem(j, i) = fem(j, i)*dmm + tdmm*tfem
            fdmon(j, i) = fdmon(j, i)*dmm + tdmm*tfdm
          end if
        end do
      end do

      jz = imitt + 1
      do iz = imitt + 1, imitt + ibl
        jz = jz - 1
        do kz = 1, mxrho + kbl
          fdmon(kz, iz) = fdmon(kz, jz)
          fem(kz, iz) = fem(kz, jz)
        end do
      end do
    end do  ! End of main iteration loop

    ! Check if maximum iterations exceeded
    if (niter .gt. ioimaxm) then
      stop
    end if

    ! Output results (formerly label 200)
    rewind 89
          rewind 78
          rewind 85
          sumW = 0.d0
          z = -0.5d0*dz
          do iz = istp1, imitt
            z = z + dz
            write (85, *) z, fdmon(1, iz), fem(1, iz)
            write (89, '(6f14.7)') z, c(1, iz, 1), c(1, iz, 3), c(1, iz, 5), &
              ehbclam(1, iz), c(1, iz, 9)
            fsum = 0.d0
            klm = nint(1.d0/drho)
            rho = -0.5d0*drho
            do i = 1, klm
              rho = rho + drho
              fsum = fsum + fdmon(i, iz)*2.d0*pi*rho
            end do
            write (78, *) z, fsum*drho/(pi*1.d0**2)
          end do

          rewind 83
          rewind 87
          iz = int(zc1*rdz) + 1
          rho = -0.5d0*drho
          do kr = 1, mxrho
            rho = rho + drho
            write (87, '(6f14.7)') rho, c(kr, iz, 1), c(kr, iz, 3), c(kr, iz, 5), &
              ehbclam(kr, iz), c(kr, iz, 9)
            write (83, *) rho, fdmon(kr, iz), fem(kr, iz)
          end do

          bfde = 2.d0*bdpol
          bfdc = bdm - bfde
          asumW = 0.d0
          sumsp = 0.d0
          sumsn = 0.d0
          bsumW = 0.d0
          chvol = 0.d0
          z = -0.5d0*dz
          do iz = istp1, imitt
            z = z + dz
            arsum = 0.d0
            brsum = 0.d0
            cv = 0.d0
            diffz2 = (z - zc1)**2
            rho = -0.5d0*drho
            do kz = 1, mxrho
              rho = rho + drho
              rsq = rho*rho + diffz2
              fdm = fdmon(kz, iz)
              if (rsq .ge. Rcoll2) then
                belamb = dlog(ebelam(kz, iz)) - emscale
                bclamb = 2.d0*(dlog(ehbclam(kz, iz)) - scalem)
                fde = fem(kz, iz)
                fdc = fdm - fde
                Fex = fdc*Y*(ae2(kz, iz) - ae1(kz, iz)) + 0.5d0*fde*ae2(kz, iz)
                arsum = &
                  rho*(fdc*bclamb + bfdc*bcmtrams + fde*belamb + bfde*bemtrams + &
                       bdpol - fdm*rrnmon + Fex - bFex) + arsum
                brsum = &
                  rho*(fdc*bclamb + fde*belamb - fdm*rrnmon + Fex - bFex) + brsum
              end if
              eexc = -dlog(edu(kz, iz))
              arsum = arsum - 0.5d0*rho*eexc*(fdm + bdm)
              brsum = brsum + rho*(0.5d0*(fdm - bdm)*eexc - fdm*eexc)
            end do
            asumW = 2.d0*pi*arsum*drho + asumW
            bsumW = 2.d0*pi*brsum*drho + bsumW
          end do
          asumW = asumW*dz
          bsumW = bsumW*dz
          aW = 2.d0*asumW
          bW = 2.d0*bsumW
          write (*, *)
          write (*, *) 'aW = ', aW
          write (*, *) 'bW = ', bW

!     net force from rho iterations:

          izmin = nint((zc1 + 0.5d0*dz - Rcoll)*rdz + 0.5d0)
          zmin = (dfloat(izmin) - 0.5d0)*dz
          izmax = nint((zc1 - 0.5d0*dz + Rcoll)*rdz + 0.5d0)
          zmax = (dfloat(izmax) - 0.5d0)*dz
          izc1 = nint(zc1*rdz + 0.5d0)
          write (*, *) 'zmin,zc1,zmax = ', zmin, zc1, zmax
          write (*, *) 'izmin,izc1,izmax = ', izmin, izc1, izmax
          write (*, *) dfloat(izmin)*dz - 0.5d0*dz, dfloat(izmax)*dz - 0.5d0*dz
          write (*, *) dfloat(izc1)*dz - 0.5d0*dz
          ict = 0

          rhoFo = 0.d0
          rcliffFo = 0.d0
          z = zmin - dz
          do iz = izmin, izc1 - 1
            z = z + dz
            zsq = (z - zc1)**2
            if (zsq .le. Rcoll2) then
              rho = -0.5d0*drho
              irho = 0
              ! Find boundary point
              do while ((rho*rho + zsq) .le. Rcoll2)
                rho = rho + drho
                irho = irho + 1
              end do
              Rc = dsqrt(rho*rho + zsq)
              rhoc = dsqrt(Rcoll2 - zsq)

            if (dabs(fdmon(irho, iz)) .gt. 0.00000001d0) then
              y3 = fdmon(irho, iz)
              y2 = fdmon(irho + 1, iz)
              y1 = fdmon(irho + 2, iz)
              x3 = rho
              x2 = rho + drho
              x1 = rho + 2.d0*drho
            else
              y3 = fdmon(irho + 1, iz)
              y2 = fdmon(irho + 2, iz)
              y1 = fdmon(irho + 3, iz)
              x3 = rho + drho
              x2 = rho + 2.d0*drho
              x1 = rho + 3.d0*drho
              write (*, *) 'TJOHO!'
            end if

            x = rhoc
            fdc = y1*(x - x2)*(x - x3)/((x1 - x2)*(x1 - x3)) + &
                  y2*(x - x1)*(x - x3)/((x2 - x1)*(x2 - x3)) + &
                  y3*(x - x1)*(x - x2)/((x3 - x1)*(x3 - x2))
              ctheta = (z - zc1)/Rcoll
              rhoFo = 2.d0*pi*rhoc*ctheta*fdc + rhoFo
              rcliffFo = 2.d0*pi*ctheta*fdc + rcliffFo
              ict = ict + 1
              ctvec(ict) = ctheta
              cdens(ict) = fdc
            end if
          end do
! 4554 write(*,*) 'rhoFo = ',rhoFo*dz
          write (*, *) 'rcliffFo = ', Rcoll*rcliffFo*dz
          write (*, *) 'z = ', z

          rhoFi = 0.d0
          rcliffFi = 0.d0
          z = zc1 - 0.5d0*dz
          do iz = izc1, izmax
            z = z + dz
            zsq = (z - zc1)**2
            if (zsq .le. Rcoll2) then
              rho = -0.5d0*drho
              irho = 0
              ! Find boundary point
              do while ((rho*rho + zsq) .le. Rcoll2)
                rho = rho + drho
                irho = irho + 1
              end do
              Rc = dsqrt(rho*rho + zsq)
              rhoc = dsqrt(Rcoll2 - zsq)

            if (dabs(fdmon(irho, iz)) .gt. 0.00000001d0) then
              y3 = fdmon(irho, iz)
              y2 = fdmon(irho + 1, iz)
              y1 = fdmon(irho + 2, iz)
              x3 = rho
              x2 = rho + drho
              x1 = rho + 2.d0*drho
            else
              y3 = fdmon(irho + 1, iz)
              y2 = fdmon(irho + 2, iz)
              y1 = fdmon(irho + 3, iz)
              x3 = rho + drho
              x2 = rho + 2.d0*drho
              x1 = rho + 3.d0*drho
              write (*, *) 'TJOHO!!!!', fdmon(irho, iz), rho
            end if

            x = rhoc
            fdc = y1*(x - x2)*(x - x3)/((x1 - x2)*(x1 - x3)) + &
                  y2*(x - x1)*(x - x3)/((x2 - x1)*(x2 - x3)) + &
                  y3*(x - x1)*(x - x2)/((x3 - x1)*(x3 - x2))
              ctheta = (z - zc1)/Rcoll
              rhoFi = 2.d0*pi*rhoc*ctheta*fdc + rhoFi
              rcliffFi = 2.d0*pi*ctheta*fdc + rcliffFi
              ict = ict + 1
              ctvec(ict) = ctheta
              cdens(ict) = fdc
            end if
          end do
! 5554 write(*,*) 'rhoFi = ',rhoFi*dz
          write (*, *) 'rcliffFi = ', Rcoll*rcliffFi*dz
          rhoF = (rhoFi + rhoFo)*dz
!      write(*,*) 'rhoF = ',rhoF
          write (*, *)
          write (*, *) 'rcliffF = ', Rcoll*(rcliffFi + rcliffFo)*dz
          write (*, *)
          write (*, *) 'z = ', z

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
          ctF = 2.d0*pi*Rcoll2*ctF
          ch2 = 2.d0*pi*Rcoll2*ch2
          write (*, *)
          write (*, *) 'ctF = ', ctF
          write (*, *)
          write (*, *) 'ch2 = ', ch2
          write (*, *)

!     net force from z iterations:

          irhomax = nint(Rcoll*rdrho + 1.d0)
          rhomax = (dfloat(irhomax) - 0.5d0)*drho
!     stay inside the sphere
          irhomax = irhomax - 1
          rhomax = rhomax - drho
          write (*, *) 'rhomax,irhomax = ', rhomax, irhomax
          ict = 0
          zFo = 0.d0
          cho = 0.d0
          cliffFo = 0.d0
          rho = -0.5d0*dz
          do irho = 1, irhomax
            rho = rho + drho
            rhosq = rho*rho
            if (rhosq .le. Rcoll2) then
              z = zc1 + 0.5d0*dz
              iz = izc1
              ! Find boundary point
              do while ((rhosq + (z - zc1)**2) .le. Rcoll2)
                z = z - dz
                iz = iz - 1
              end do
              zsq = (z - zc1)**2
              Rc = dsqrt(rhosq + zsq)
              deltazc = dsqrt(Rcoll2 - rhosq)

            if (dabs(fdmon(irho, iz)) .gt. 0.00000001d0) then
              y3 = fdmon(irho, iz)
              y2 = fdmon(irho, iz - 1)
              y1 = fdmon(irho, iz - 2)
              x3 = dabs(z - zc1)
              x2 = x3 + dz
              x1 = x3 + 2.d0*dz
            else
              y3 = fdmon(irho, iz - 1)
              y2 = fdmon(irho, iz - 2)
              y1 = fdmon(irho, iz - 3)
              x3 = dabs(z - zc1) + dz
              x2 = x3 + dz
              x1 = x3 + 2.d0*dz
              write (*, *) 'TJOHO1!'
            end if

            x = deltazc
            fdc = y1*(x - x2)*(x - x3)/((x1 - x2)*(x1 - x3)) + &
                  y2*(x - x1)*(x - x3)/((x2 - x1)*(x2 - x3)) + &
                  y3*(x - x1)*(x - x2)/((x3 - x1)*(x3 - x2))
              ctheta = -deltazc/Rcoll
              ict = ict + 1
              ctvec(ict) = ctheta
              cdens(ict) = fdc
            end if
          end do

          rho = rho + drho
          irho = irhomax + 1
!      ict = 0
          zFi = 0.d0
          chi = 0.d0
          cliffFi = 0.d0
!      rho = -0.5d0*dz
          do krho = 1, irhomax
!      rho = rho+drho
            irho = irho - 1
            rho = rho - drho
            rhosq = rho*rho
            if (rhosq .le. Rcoll2) then
              z = zc1 - 0.5d0*dz
              iz = izc1 - 1
              ! Find boundary point
              do while ((rhosq + (z - zc1)**2) .le. Rcoll2)
                z = z + dz
                iz = iz + 1
              end do
              zsq = (z - zc1)**2
              Rc = dsqrt(rhosq + zsq)
              deltazc = dsqrt(Rcoll2 - rhosq)

            if (dabs(fdmon(irho, iz)) .gt. 0.00000001d0) then
              y3 = fdmon(irho, iz)
              y2 = fdmon(irho, iz + 1)
              y1 = fdmon(irho, iz + 2)
              x3 = dabs(z - zc1)
              x2 = x3 + dz
              x1 = x3 + 2.d0*dz
            else
              y3 = fdmon(irho, iz + 1)
              y2 = fdmon(irho, iz + 2)
              y1 = fdmon(irho, iz + 3)
              x3 = dabs(z - zc1) + dz
              x2 = x3 + dz
              x1 = x3 + 2.d0*dz
              write (*, *) 'TJOHO2!'
            end if

            x = deltazc
            fdc = y1*(x - x2)*(x - x3)/((x1 - x2)*(x1 - x3)) + &
                  y2*(x - x1)*(x - x3)/((x2 - x1)*(x2 - x3)) + &
                  y3*(x - x1)*(x - x2)/((x3 - x1)*(x3 - x2))
            ctheta = deltazc/Rcoll
              zFi = 2.d0*pi*rho*ctheta*fdc + zFi
              chi = 2.d0*pi*rho*ctheta + chi
              cliffFi = 2.d0*pi*ctheta*fdc + cliffFi
              ict = ict + 1
              ctvec(ict) = ctheta
              cdens(ict) = fdc
            end if
          end do

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
          ckoll = 0.d0
          cckoll = 0.d0
          ccckoll = 0.d0
          ch2 = 0.d0
          ccc = 0.d0
          bordekoll = 0.d0
          do kct = 0, ict
            ct = ctvec(kct)
            ctn = ctvec(kct + 1)
            fdc = cdens(kct)
            fdcn = cdens(kct + 1)
            fk = (fdcn - fdc)/(ctn - ct)
            ctF = 0.5d0*(fdc - fk*ct)*(ctn**2 - ct**2) + fk*(ctn**3 - ct**3)/3.d0 + ctF
            ch2 = 0.25d0*(fdc + fdcn)*(ctn**2 - ct**2) + ch2

            tn = -(1.d0 - ctn*ctn)**1.5d0
            t = -(1.d0 - ct*ct)**1.5d0
            add = 0.5d0*(fdc + fdcn)*(tn - t)/3.d0
            bordekoll = add + bordekoll

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
            rho = Rcoll*dsqrt(1.d0 - ct*ct)
            rhon = Rcoll*dsqrt(1.d0 - ctn*ctn)
            if (kct .eq. 0) rho = Rcoll
            if (kct .eq. ict) rhon = Rcoll
            if (dabs(rhon - rho) .lt. 0.00000001d0) then
              fk = 0.d0
            else
              fk = (fdcn - fdc)/(rhon - rho)
            end if
            ccc = (fdc - fk*rho)*(rhon - rho) + 0.5d0*fk*(rhon*rhon - rho*rho) + ccc
          end do

          ctF = 2.d0*pi*Rcoll2*ctF
          ch2 = 2.d0*pi*Rcoll2*ch2
          cckoll = 2.d0*pi*cckoll*drho
          ckoll = 2.d0*pi*ckoll*Rcoll*drho
          ckk = 2.d0*pi*ckk*Rcoll*drho
          ccc = 2.d0*pi*ccc
          write (*, *)
          write (*, *) 'ctF = ', ctF
          write (*, *)
          write (*, *) 'ch2 = ', ch2

          rewind ifc
          z = -0.5d0*dz
          do iz = istp1, imitt
            z = z + dz
            rho = -0.5d0*drho
            do kz = 1, mxrho
              rho = rho + drho
              write (ifc, '(2f12.5,2f21.12)') &
                z, rho, fdmon(kz, iz), fem(kz, iz)
            end do
          end do

          STOP
        END

        subroutine CDFACT
          implicit double precision(a - h, o - z)
          include 't2.inc.f90'
          strho0 = 0.5d0*drho
          rho0 = strho0
          iz = 2*ism
!      z = 2.d0-dz
!      zpst = z-1.d0-dz
          z = 2.d0*dhs - dz
          zpst = z - dhs - dz
          sume = 0.d0
          zp = zpst
          do jz = iz - ism, iz + ism
            zp = zp + dz
            delz2 = (zp - z)**2
            sumrhop = 0.d0
!      rhopmax = dsqrt(dabs(1.d0-delz2))
            rhopmax = dsqrt(dabs(dhs2 - delz2))
            krhopmax = nint(rhopmax*rdrho)
            rho = rho0
            rho02 = rho0**2
            rhop = -0.5d0*drho
            do krhop = 1, krhopmax
              rhop = rhop + drho
              rhomax2 = rho02 + rhop*rhop
              fphi = 2.d0*rho0*rhop
              phisum = 0.d0
              phi = -0.5d0*dphi
              do iphi = 1, nphi
                phi = phi + dphi
!     plus eller minus spelar ingen roll foer integralens vaerde
                rho2 = rhomax2 - fphi*dcos(phi)
                rho = dsqrt(rho2)
                irho = int(rho*rdrho) + 1
                phisum = 1.d0 + phisum
              end do
              sumrhop = rhop*phisum*dphi + sumrhop
            end do
            write (*, *) rho, rhopmax, dabs(zp - z)
            fact = 1.d0
            if (iabs(jz - iz) .eq. ism) fact = 0.5d0
            sume = 2.d0*sumrhop*drho*fact + sume
          end do
          tcd = 3.d0*sume*dzrfp*rdhs3
          write (*, *) 'tcd = ', tcd, nphi
          cdnorm = 1.d0/tcd
          return
        end

        subroutine CDCALC
          implicit double precision(a - h, o - z)
          include 't2.inc.f90'
!      z = 1.d0-0.5d0*dz
          z = dhs - 0.5d0*dz
          do iz = istp1 + ism, imitt
            z = z + dz
!      zpst = z-1.d0-dz
            zpst = z - dhs - dz
            rho0 = -0.5d0*drho
            do kz = 1, mxrho - ksm
              rho0 = rho0 + drho
              sume = 0.d0
              zp = zpst
              do jz = iz - ism, iz + ism
                zp = zp + dz
                delz2 = (zp - z)**2
                sumrhop = 0.d0
!      rhopmax = dsqrt(dabs(1.d0-delz2))
                rhopmax = dsqrt(dabs(dhs2 - delz2))
                krhopmax = nint(rhopmax*rdrho)
                rho02 = rho0**2
                rhop = -0.5d0*drho
                do krhop = 1, krhopmax
                  rhop = rhop + drho
                  rhomax2 = rho02 + rhop*rhop
                  fphi = 2.d0*rho0*rhop
                  phisum = 0.d0
                  phi = -0.5d0*dphi
                  do iphi = 1, nphi
                    phi = phi + dphi
!     plus eller minus spelar ingen roll foer integralens vaerde
                    rho2 = rhomax2 - fphi*dcos(phi)
                    rho = dsqrt(rho2)
                    irho = int(rho*rdrho) + 1
                    phisum = fdmon(irho, jz) + phisum
                  end do
                  sumrhop = rhop*phisum*dphi + sumrhop
                end do
                fact = 1.d0
                if (iabs(jz - iz) .eq. ism) fact = 0.5d0
                sume = 2.d0*sumrhop*drho*fact + sume
              end do
              cdmonm(kz, iz) = 3.d0*sume*dzrfp*cdnorm*rdhs3
!     rho0 svepet faerdigt!
            end do
!     dags foer nytt z0 vaerde!
          end do
          return
        end

          subroutine AVEC
            implicit double precision(a - h, o - z)
            include 't2.inc.f90'
            do iz = istp1 + 2*ism, imitt
!      do kz = 1,mxrho-ksm
              do kz = 1, mxrho - kbl
                cdt = cdmonm(kz, iz)*dhs3
                pcdt = pis*cdt
                xsi = (1.d0 - pcdt)
                rxsi = 1.d0/xsi
                sqrxsi = rxsi*rxsi
                flog = dlog(xsi)
                ae1(kz, iz) = -(c1 + 1.d0)*flog - 0.5d0*(AA1 + BB1*pcdt)*pcdt*sqrxsi
                ae2(kz, iz) = -(c2 + 1.d0)*flog - 0.5d0*(AA2 + BB2*pcdt)*pcdt*sqrxsi
                daex1 = rxsi*(c1 + 1.d0 - 0.5d0*AA1*rxsi*(1.d0 + 2.d0*pcdt*rxsi) - &
                              BB1*pcdt*rxsi*(1.d0 + pcdt*rxsi))
                daex2 = rxsi*(c2 + 1.d0 - 0.5d0*AA2*rxsi*(1.d0 + 2.d0*pcdt*rxsi) - &
                              BB2*pcdt*rxsi*(1.d0 + pcdt*rxsi))
                convp(kz, iz) = (Y*(fdmon(kz, iz) - fem(kz, iz))*(daex2 - daex1) + &
                                 0.5d0*fem(kz, iz)*daex2)*pis*dhs3
!     *0.5d0*fem(kz,iz)*daex2)*pis
              end do
            end do

!      do kz = 1,mxrho-ksm
            do kz = 1, mxrho - kbl
              jz = imitt + 1
!      do iz = imitt+1,imitt+ism
              do iz = imitt + 1, imitt + ibl
                jz = jz - 1
                ae1(kz, iz) = ae1(kz, jz)
                ae2(kz, iz) = ae2(kz, jz)
                convp(kz, iz) = convp(kz, jz)
              end do
            end do
            return
          end

          subroutine EBLMNEW
            implicit double precision(a - h, o - z)
            include 't2.inc.f90'
            z = -0.5d0*dz
            do iz = 1, ibl
              z = z + dz
              do kz = 1, mxrho + kbl
                ebelam(kz, iz) = bebelam
                ehbclam(kz, iz) = behbclam
              end do
            end do

            do iz = ibl + 1, imitt
              z = z + dz
              jstart = iz - ism
!      zpst = z-1.d0-dz
              zpst = z - dhs - dz
              diffz2 = (zc1 - z)**2
              irho0min = 1
              strho0 = -0.5d0*drho

              rho0 = strho0
              do kz = irho0min, mxrho
                rho0 = rho0 + drho
                rho02 = rho0**2

                rt2 = rho02 + (z - zc1)**2
                if (rt2 .lt. Rcoll2) then
                  ehbclam(kz, iz) = 0.d0
                  ebelam(kz, iz) = 0.d0
                  cycle
                end if
                rt2 = rho02 + (z - zc2)**2
                if (rt2 .lt. Rcoll2) then
                  ehbclam(kz, iz) = 0.d0
                  ebelam(kz, iz) = 0.d0
                  cycle
                end if

                sume = 0.d0
                zp = zpst
                do jz = jstart, iz + ism
                  zp = zp + dz
                  delz2 = (zp - z)**2
                  zpcsq = (zp - zc1)**2
                  zpc2sq = (zp - zc2)**2

                  sumrhop = 0.d0
!      rhopmax = dsqrt(dabs(1.d0-delz2))
                  rhopmax = dsqrt(dabs(dhs2 - delz2))
                  krhopmax = nint(rhopmax*rdrho)
                  rho02 = rho0**2
                  rhop = -0.5d0*drho
                  do krhop = 1, krhopmax
                    rhop = rhop + drho
                    rhomax2 = rho02 + rhop*rhop
                    fphi = 2.d0*rho0*rhop
                    phisum = 0.d0
                    phi = -0.5d0*dphi
                    do iphi = 1, nphi
                      phi = phi + dphi
!     plus eller minus spelar ingen roll foer integralens vaerde
                      rho2 = rhomax2 - fphi*dcos(phi)
                      rsq = rho2 + zpcsq
                      if (rsq .lt. Rcoll2) cycle
                      rsq = rho2 + zpc2sq
                      if (rsq .lt. Rcoll2) cycle
                      rho = dsqrt(rho2)
                      irho = int(rho*rdrho) + 1
                      phisum = convp(irho, jz) + phisum
                    end do
                    sumrhop = rhop*phisum*dphi + sumrhop
                  end do
                  fact = 1.d0
                  if (iabs(jz - iz) .eq. ism) fact = 0.5d0
                  sume = 2.d0*sumrhop*drho*fact + sume
                end do
                trams = 3.d0*sume*dzrfp*cdnorm*rdhs3

                emtrams = trams + 0.5d0*ae2(kz, iz)
                cmtrams = trams + Y*(ae2(kz, iz) - ae1(kz, iz))
                ebelam(kz, iz) = dexp(-emtrams + emscale)
                ehbclam(kz, iz) = dexp(-0.5d0*(cmtrams) + scalem)
!     rho0 svepet faerdigt!
              end do
!     dags foer nytt z0 vaerde!
            end do
              return
            end

            subroutine EBDU
              implicit double precision(a - h, o - z)
              include 't2.inc.f90'
              z = -0.5d0*dz
              do iz = 1, ibl
                z = z + dz
                do kz = 1, mxrho + kbl
                  edu(kz, iz) = 1.d0
                  edu(kz, iz) = 1.d0
                end do
              end do

              do iz = ibl + 1, imitt
!     z loop
                z = z + dz
                rho = -0.5d0*drho
                do krho = 1, mxrho
!     rho loop
                  rho = rho + drho
                  sumpint = 0.d0
                  tz = bl - 0.5d0*dz
                  do ipz = ibl + 1, imitt
!     zp loop, LHS
                    tz = tz + dz
                    tdz = z - tz
                    itdz = nint(dabs(tdz*rdz))
                    sumrho = 0.d0
                    do kprho = 1, mxrho
                      sumrho = (fdmon(kprho, ipz) - bdm)*hvec(itdz, krho, kprho) + sumrho
                    end do
                    sumpint = 2.d0*sumrho*drho + sumpint
!     zp loop, LHS, finished
                  end do

                  do ipz = imitt + 1, nfack - ibl
!     zp loop, RHS
                    tz = tz + dz
                    tdz = z - tz
                    itdz = nint(dabs(tdz*rdz))
                    sumrho = 0.d0
                    do kprho = 1, mxrho
                      sumrho = &
                        (fdmon(kprho, nfack + 1 - ipz) - bdm)*hvec(itdz, krho, kprho) + sumrho
                    end do
                    sumpint = 2.d0*sumrho*drho + sumpint
!     zp loop finished
                  end do
                  sumpint = sumpint*dz
                  edu(krho, iz) = dexp(-sumpint)
!     rho loop finished
                end do
!     z loop finished
              end do
                return
              end

