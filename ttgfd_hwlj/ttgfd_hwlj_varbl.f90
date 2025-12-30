! ============================================================================
! TTGFD_HWLJ - Generalized Flory-Dimer Theory for Polymer Solutions
! ============================================================================
! Calculates polymer density profiles near hard-wall surfaces with
! Lennard-Jones interactions using density functional theory.
! Implements the Generalized Flory-Dimer (GFD) approximation for
! inhomogeneous polymer solutions in cylindrical geometry.
! ============================================================================
program platem
  implicit none
  include 't2.inc.f90'

  ! Constants
  INTEGER, PARAMETER :: MAXMON = 151

  ! Additional arrays for main program - also dynamically allocated
  double precision, allocatable :: c(:, :, :), cA(:, :), cB(:, :)
  double precision :: cdens(0:1000), ctvec(0:1000)

  ! Integer variables
  integer :: i, j, k, iz, jz, kz, irho, krho, kprho, iphi, imon, kmon
  integer :: itdz, ict, kct, irdz, izc1, izmin, izmax, irho0min, irhomax, jstart
  integer :: ifc, ins, iep, ioimaxm, kread, niter, kcm, klm, kr, npphi

  ! Double precision variables
  double precision :: add, aeta, aex1, aex2, alj, arsum, asumw, aw
  double precision :: baex1, baex2, bclamb, bcmtrams, bconvp, bdaex1, bdaex2
  double precision :: bdpol, bds, bdt, bebbe, belamb, bemtrams, bfdc, bfde, bfex
  double precision :: bl2, bordekoll, brsum, bsumw, bw
  double precision :: ccc, ccckoll, cckoll, cdt, ch2, chempp, chi, cho, chvol
  double precision :: ckk, ckoll, clifffi, clifffo, cmtrams, collsep, ct, ctf
  double precision :: ctheta, ctn, ctp, cv
  double precision :: daex1, daex2, ddiff, ddmax, deltazc, delz2, diffz2, distance
  double precision :: dmm, dms, dmm_adaptive, dms_adaptive, dpphi, dumsum
  double precision :: interface_width, smooth_factor
  double precision :: eexc, efact, emtrams, epslj

  ! Logical variables
  logical :: use_adaptive_mixing
  integer :: oscillation_count
  double precision :: ddmax_prev
  double precision :: fact, fdc, fdcm1, fdcn, fdcp1, fde, fdm, fex, ffact, fk, flog, fphi, fsum
  double precision :: pb, pcdt, pdasum, phi, phisum, pint
  double precision :: rc, rclifffi, rclifffo, rcyl, rcyl2, rdphi
  double precision :: rho, rho0, rho02, rho2, rhoc, rhof, rhofi, rhofo, rhomax, rhon, rhosq, rhoz2, rlj
  double precision :: rsq, rt2, rxsi, rxsib, rxsibsq, s2, sqrxsi, strho0
  double precision :: sume, sumsn, sumsp, sumw
  double precision :: t, t1, t2, tdmm, tdms, tdz, tdzsq, tfdm, tfem, th, tn, trams
  double precision :: trho, trhosq, trmix, twopidz, use1, useful
  double precision :: x, x1, x2, x3, xsi, xsib, y1, y2, y3
  double precision :: z, z2, z22, zfact, zfi, zfo, zmax, zmin, zp, zpc2sq, zpcsq, zpst, zsq

  ! ========================================================================
  ! File unit numbers
  ifc = 38  ! Output: concentration profiles
  ins = 49  ! Input: simulation parameters
  iep = 50  ! Input: Lennard-Jones epsilon parameter

  ! Open input and output files
  open (ifc, file='fcdfil', form='formatted')
  open (ins, file='input.tsph', form='formatted')
  open (iep, file='epfil', form='formatted')
  rewind ifc
  rewind ins
  rewind iep

  ! Read simulation parameters from input file
  read (ins, *) bdm         ! Monomer bulk density
  read (ins, *) bdtot       ! Total bulk density
  bds = bdtot - bdm         ! Solvent bulk density
  read (ins, *) nmon        ! Number of monomers per polymer chain
  read (ins, *) dz          ! Grid spacing in z direction
  read (ins, *) drho        ! Grid spacing in radial direction
  read (ins, *) dphi        ! Angular grid spacing (input in units of pi)
  dphi = PI*dphi
  read (ins, *) Rcoll       ! Colloid radius
  read (ins, *) zc1         ! Position of first colloid center
  read (ins, *) collsep     ! Separation between colloid centers
  read (ins, *) Rcyl        ! Cylinder radius (system boundary)
  read (ins, *) ioimaxm     ! Maximum number of iterations
  read (ins, *) dmm, dms    ! Density mixing parameters (monomer, solvent)
  read (ins, *) kread       ! Read initial guess from file (0=no, 1=yes)
  read (ins, *) bl          ! Bond length
  read (ins, *) dhs         ! Hard sphere diameter (monomer)
  read (ins, *) dpphi       ! Angular grid spacing for potential calculation
  ! Compute derived geometric parameters
  bl2 = bl*bl
  dhs2 = dhs*dhs
  dhs3 = dhs2*dhs
  rdhs3 = 1.d0/dhs3

  ! Read Lennard-Jones energy parameter and compute LJ coefficients
  read (iep, *) epslj
  alj = 4.d0*epslj*dhs3*dhs3    ! Attractive (r^-6) coefficient
  rlj = 4.d0*epslj*dhs3**4      ! Repulsive (r^-12) coefficient

  ! Compute additional derived parameters
  Rcyl2 = Rcyl*Rcyl
  tdmm = 1.d0 - dmm
  tdms = 1.d0 - dms
  rdz = 1.d0/dz
  rdrho = 1.d0/drho
  rdphi = 1.d0/dphi
  twopidz = TWOPI*dz
  dzrfp = dz/(4.d0*PI)

  ! Discretization parameters
  nphi = int(PI/dphi + 0.01d0)

  ! Initialize cosine lookup table for dphi
  do iphi = 1, nphi
    phi = (dble(iphi) - 0.5d0)*dphi
    cos_phi(iphi) = dcos(phi)
  end do

  irdz = int(rdz + 0.001d0)
  nfack = int(2.d0*(zc1 + 0.5d0*collsep)/dz + 0.01d0)
  istart = 0
  istp1 = 1
  islut = nfack
  istp1s = 1
  isluts = islut
  ism = int(dhs/dz + 0.01d0)      ! Hard sphere diameter in grid units (z)
  ksm = int(dhs/drho + 0.01d0)    ! Hard sphere diameter in grid units (rho)
  ibl = int(bl/dz + 0.01d0)       ! Bond length in grid units (z)
  kbl = int(bl/drho + 0.01d0)     ! Bond length in grid units (rho)
  rnmon = dble(nmon)
  rrnmon = 1.d0/rnmon

  ! Bulk thermodynamic properties
  Yfact = (rnmon - 2.d0)*Y
  bdpol = bdm/rnmon  ! Polymer bulk density

  ! Hard sphere packing fraction and related quantities
  bdt = bdm*dhs3
  aeta = PIS*bdt
  xsib = 1.d0 - aeta
  rxsib = 1.d0/xsib
  rxsibsq = rxsib*rxsib

  ! Excess free energy terms (Carnahan-Starling)
  aex1 = -(C1 + 1.d0)*dlog(xsib) - &
         0.5d0*(AA1*PIS*bdt + BB1*(PIS*bdt)**2)*rxsibsq
  aex2 = -(C2 + 1.d0)*dlog(xsib) - &
         0.5d0*(AA2*PIS*bdt + BB2*(PIS*bdt)**2)*rxsibsq
  bFex = (bdm - 2.d0*bdpol)*Y*(aex2 - aex1) + bdpol*aex2

  ! Derivatives of excess free energy
  daex1 = rxsib*(C1 + 1 - 0.5d0*(AA1 + 2.d0*BB1*aeta)*rxsib - &
                 aeta*(AA1 + BB1*aeta)*rxsibsq)
  daex2 = rxsib*(C2 + 1 - 0.5d0*(AA2 + 2.d0*BB2*aeta)*rxsib - &
                 aeta*(AA2 + BB2*aeta)*rxsibsq)
  pdasum = Yfact*(daex2 - daex1) + daex2
  baex1 = aex1
  baex2 = aex2
  bdaex1 = daex1
  bdaex2 = daex2

  ! Bulk pressure and chemical potentials
  Pb = bdpol + bdpol*aeta*pdasum
  chempp = dlog(bdpol) + Yfact*(aex2 - aex1) + aex2 + PIS*bdm*pdasum*dhs3
  scalem = chempp/(2.d0*rnmon)
  emscale = 2.d0*scalem

  ! Convolution terms for inhomogeneous density functional
  bconvp = (Y*(bdm - 2.d0*bdm*rrnmon)*(daex2 - daex1) + &
            bdm*rrnmon*daex2)*PIS*dhs3
  trams = bconvp
  emtrams = trams + 0.5d0*aex2
  cmtrams = trams + Y*(aex2 - aex1)
  bemtrams = emtrams
  bcmtrams = cmtrams
  bebelam = dexp(-emtrams + emscale)
  behbclam = dexp(-0.5d0*cmtrams + scalem)

  ! Print simulation parameters
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

  ! Allocate main program arrays (COMMON block arrays are already static)
  allocate (c(0:maxrho, 0:maxel, MAXMON))
  allocate (cA(0:maxrho, 0:maxel))
  allocate (cB(0:maxrho, 0:maxel))

  ! Calculate normalization constant for contact density
  CALL CDFACT
  write (*, *) 'cdnorm = ', cdnorm

  ! Initialize excess free energy arrays near z-boundaries (left side)
  ! These regions are near the system edge and require special treatment
  do iz = istp1, istp1 + 2*ibl - 1
  do kz = 1, mxrho + kbl
    cdt = bdm*dhs3
    pcdt = PIS*cdt
    xsi = (1.d0 - pcdt)
    rxsi = 1.d0/xsi
    sqrxsi = rxsi*rxsi
    flog = dlog(xsi)
    ae1(kz, iz) = -(C1 + 1.d0)*flog - 0.5d0*(AA1 + BB1*pcdt)*pcdt*sqrxsi
    ae2(kz, iz) = -(C2 + 1.d0)*flog - 0.5d0*(AA2 + BB2*pcdt)*pcdt*sqrxsi
    daex1 = rxsi*(C1 + 1.d0 - 0.5d0*AA1*rxsi*(1.d0 + 2.d0*pcdt*rxsi) - &
                  BB1*pcdt*rxsi*(1.d0 + pcdt*rxsi))
    daex2 = rxsi*(C2 + 1.d0 - 0.5d0*AA2*rxsi*(1.d0 + 2.d0*pcdt*rxsi) - &
                  BB2*pcdt*rxsi*(1.d0 + pcdt*rxsi))
    convp(kz, iz) = (Y*(bdm - 2.d0*bdm*rrnmon)*(daex2 - daex1) + &
                     bdm*rrnmon*daex2)*PIS*dhs3
  end do
  end do
  ! Initialize excess free energy arrays near radial boundaries (outer edge)
  do iz = istp1 + 2*ibl, imitt + ibl
  do kz = mxrho - kbl + 1, mxrho + kbl
    cdt = bdm*dhs3
    pcdt = PIS*cdt
    xsi = (1.d0 - pcdt)
    rxsi = 1.d0/xsi
    sqrxsi = rxsi*rxsi
    flog = dlog(xsi)
    ae1(kz, iz) = -(C1 + 1.d0)*flog - 0.5d0*(AA1 + BB1*pcdt)*pcdt*sqrxsi
    ae2(kz, iz) = -(C2 + 1.d0)*flog - 0.5d0*(AA2 + BB2*pcdt)*pcdt*sqrxsi
    daex1 = rxsi*(C1 + 1.d0 - 0.5d0*AA1*rxsi*(1.d0 + 2.d0*pcdt*rxsi) - &
                  BB1*pcdt*rxsi*(1.d0 + pcdt*rxsi))
    daex2 = rxsi*(C2 + 1.d0 - 0.5d0*AA2*rxsi*(1.d0 + 2.d0*pcdt*rxsi) - &
                  BB2*pcdt*rxsi*(1.d0 + pcdt*rxsi))
    convp(kz, iz) = (Y*(bdm - 2.d0*bdm*rrnmon)*(daex2 - daex1) + &
                     bdm*rrnmon*daex2)*PIS*dhs3
  end do
  end do

  ! Initialize density fields: either from scratch or read from file
  if (kread .eq. 0) then
    ! Starting from bulk values with excluded volume for colloids
    z = -0.5d0*dz
    ! Skip boundary region at z < 0
    do iz = 1, ibl
      z = z + dz
    end do
    ! Initialize all grid points to bulk values, then zero out colloid interiors
    do iz = ibl + 1, imitt
      z = z + dz
      z2 = (z - zc1)**2
      z22 = (z - zc2)**2
      rho = -0.5d0*drho
      ! Loop over radial positions
      do kz = 1, mxrho
        rho = rho + drho
        fdmon(kz, iz) = bdm
        fem(kz, iz) = 2.d0*fdmon(kz, iz)*rrnmon
        ebelam(kz, iz) = bebelam
        ehbclam(kz, iz) = behbclam
        ! Check if point is inside first colloid
        rt2 = rho*rho + z2
        if (rt2 .lt. Rcoll2) then
          fdmon(kz, iz) = 0.d0
          fem(kz, iz) = 0.d0
          ebelam(kz, iz) = 0.d0
          ehbclam(kz, iz) = 0.d0
        end if
        ! Check if point is inside second colloid
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
    ! Read initial guess from file (restart from previous calculation)
    rewind ifc
    do iz = istp1, imitt
    do kz = 1, mxrho
      read (ifc, *) t1, t2, fdmon(kz, iz), fem(kz, iz)
    end do
    end do
  end if

  ! Set boundary conditions at z-boundaries (left edge)
  ! All densities set to bulk values
  do iz = istp1, ibl
  do kz = 1, mxrho + kbl
    fdmon(kz, iz) = bdm
    fem(kz, iz) = 2.d0*fdmon(kz, iz)*rrnmon
    ebelam(kz, iz) = bebelam
    ehbclam(kz, iz) = behbclam
    cdmonm(kz, iz) = bdm
  end do
  end do
  ! Set boundary conditions at radial edge (outer cylinder boundary)
  ! All densities set to bulk values
  do iz = 1, imitt
  do kz = mxrho + 1, mxrho + kbl
    fdmon(kz, iz) = bdm
    fem(kz, iz) = 2.d0*fdmon(kz, iz)*rrnmon
    ebelam(kz, iz) = bebelam
    ehbclam(kz, iz) = behbclam
    cdmonm(kz, iz) = bdm
  end do
  end do

  ! Apply symmetry boundary conditions at z = imitt midplane
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

  ! Precompute Lennard-Jones interaction potential on grid (hvec array)
  ! This tabulates U_LJ for all distance combinations to speed up later calculations
  write (*, *) 'dpphi = ', dpphi
  dpphi = dpphi*PI
  npphi = int(PI/dpphi + 0.01d0)
  write (*, *) 'dpphi,npphi = ', dpphi, npphi

  ! Initialize cosine lookup table for dpphi
  do iphi = 1, npphi
    phi = (dble(iphi) - 0.5d0)*dpphi
    cos_pphi(iphi) = dcos(phi)
  end do

  kcm = nfack
  ! Triple loop over rho, rho', z to compute pairwise LJ interaction integrals
  ! Loop order matches F77 for numerical consistency
!$omp parallel do private(tdz, tdzsq, rho, rhosq, use1, trho, trhosq, trmix, useful, pint, iphi, s2, krho, kprho) schedule(static)
  do itdz = 0, kcm - 1
    tdz = -dz + dble(itdz + 1)*dz
    tdzsq = tdz*tdz
    rho = -0.5d0*drho
    do krho = 1, mxrho
      rho = rho + drho
      rhosq = rho*rho
      use1 = tdzsq + rhosq
      trho = -0.5d0*drho
      do kprho = 1, mxrho
        trho = trho + drho
        trhosq = trho*trho
        trmix = 2.d0*trho*rho
        useful = use1 + trhosq
        pint = 0.d0
        ! Angular integration for cylindrical geometry
!$omp simd reduction(+:pint)
        do iphi = 1, npphi
          s2 = useful - trmix*cos_pphi(iphi)
          if (s2 .gt. dhs2) then
            ! Lennard-Jones potential: U(r) = 4*epsilon*[(sigma/r)^12 - (sigma/r)^6]
            pint = rlj/s2**6 - alj/s2**3 + pint
          end if
        end do
!$omp end simd
        hvec(kprho, krho, itdz) = trho*pint*dpphi
      end do
    end do
  end do
!$omp end parallel do
  write (*, *) 'hvec fixad'

  ! Initialize iteration
  ddmax = 10000.d0
  ddmax_prev = 10000.d0
  niter = 0
  oscillation_count = 0
  use_adaptive_mixing = .true.

  ! Main self-consistent field iteration loop
  do while (.true.)
    niter = niter + 1
    write (*, *) 'ddmax,niter = ', ddmax, niter
    write (*, *)
    if (niter .gt. ioimaxm) then
      write (*, *) 'NITER.GT.IOIMAXM !', niter
      exit
    end if

    ! Update fields in density functional theory calculation
    CALL CDCALC     ! Calculate contact density
    CALL AVEC       ! Calculate excess free energy
    CALL EBLMNEW    ! Calculate end-segment Boltzmann factors
    CALL EBDU       ! Calculate external potential contribution

    ! Apply boundary conditions at outer radial edge (rho > Rcyl)
    ! Set to bulk values since density should approach bulk far from colloids
    do iz = istp1, imitt
    do kz = mxrho + 1, mxrho + kbl
      ebelam(kz, iz) = bebelam
      ehbclam(kz, iz) = behbclam
      edu(kz, iz) = 1.d0
    end do
    end do

    ! Apply symmetry boundary conditions at z = imitt (midplane between colloids)
    ! The system is symmetric about the midplane, so mirror the field values
    ! This loop copies from iz = imitt down to iz = 1 (reverse order)
    jz = imitt + 1
    do iz = imitt + 1, imitt + ibl
      jz = jz - 1
      do kz = 1, mxrho + kbl
        ebelam(kz, iz) = ebelam(kz, jz)
        ehbclam(kz, iz) = ehbclam(kz, jz)
        edu(kz, iz) = edu(kz, jz)
      end do
    end do

    ! Second symmetry application (appears redundant but ensures consistency)
    jz = imitt + 1
    do iz = imitt + 1, imitt + ibl
      jz = jz - 1
      do kz = 1, mxrho + kbl
        ebelam(kz, iz) = ebelam(kz, jz)
        ehbclam(kz, iz) = ehbclam(kz, jz)
        edu(kz, iz) = edu(kz, jz)
      end do
    end do

    ! Calculate cA: the propagator for polymer end segments
    ! cA(r) = exp(-beta*mu_end)*exp(-U_LJ) where:
    !   ebelam = exp(-beta*mu_end) from hard-sphere and chain connectivity
    !   edu = exp(-U_LJ) from Lennard-Jones interactions
    do iz = istp1, imitt + ibl
    do kz = 1, mxrho + kbl
      cA(kz, iz) = ebelam(kz, iz)*edu(kz, iz)
    end do
    end do

    ! Propagate polymer chains segment by segment from end to end
    ! This calculates c(r,i) = propagator for segment i at position r
    imon = nmon
    do kmon = 1, nmon - 1
      imon = imon - 1

      ! Loop over all spatial grid points
!$omp parallel do private(z, jstart, zpst, irho0min, strho0, rho0, kz, rho02, rt2, sume, zp, jz, delz2, zpcsq, zpc2sq, phisum, zfact, rhoz2, fphi, iphi, rho2, rsq, rho, irho, fact, efact, ffact) schedule(static)
      do iz = istp1 + ibl, imitt
        z = bl - 0.5d0*dz + dble(iz - (istp1 + ibl) + 1)*dz
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

          ! Integrate over bond orientations: sum contributions from all points
          ! within bond length bl of current position (rho0, z)
          sume = 0.d0
          zp = zpst
          do jz = jstart, iz + ibl
            zp = zp + dz
            delz2 = (zp - z)**2
            zpcsq = (zp - zc1)**2
            zpc2sq = (zp - zc2)**2
            phisum = 0.d0
            zfact = dabs(bl2 - delz2)
            rhoz2 = rho0**2 + zfact
            fphi = 2.d0*rho0*dsqrt(zfact)
!$omp simd reduction(+:phisum)
            do iphi = 1, nphi
!     Plus or minus sign doesn't matter for the value of the integral
              rho2 = rhoz2 - fphi*cos_phi(iphi)
              rsq = rho2 + zpcsq
              if (rsq .lt. Rcoll2) cycle  ! Inside first colloid
              rsq = rho2 + zpc2sq
              if (rsq .lt. Rcoll2) cycle  ! Inside second colloid
              rho = dsqrt(rho2)
              irho = int(rho*rdrho) + 1
              phisum = cA(irho, jz) + phisum
            end do
!$omp end simd
            fact = 1.d0
            if (iabs(jz - iz) .eq. ibl) fact = 0.5d0
            sume = 2.d0*phisum*dphi*fact + sume
          end do
          efact = dsqrt(edu(kz, iz))
          ffact = sume*dzrfp*ehbclam(kz, iz)/bl*efact
          c(kz, iz, imon) = ffact
          if (iz .gt. imitt - ibl - 1) cB(kz, islut + 1 - iz) = ffact*ehbclam(kz, iz)*efact
          cB(kz, iz) = ffact*ehbclam(kz, iz)*efact
        end do
      end do
!$omp end parallel do

      ! Handle boundary regions: propagators at z-boundaries
!$omp parallel do private(kz, bebbe)
      do iz = istp1, ibl
      do kz = 1, mxrho + kbl
        bebbe = behbclam*cA(kz, iz)
        c(kz, iz, imon) = bebbe
        cA(kz, iz) = behbclam*bebbe
      end do
      end do
!$omp end parallel do

      ! Handle radial boundaries and update propagators
!$omp parallel do private(kz, bebbe)
      do iz = ibl + 1, imitt
        ! Outer radial boundary: use bulk propagators
        do kz = mxrho - kbl, mxrho + kbl
          bebbe = behbclam*cA(kz, iz)
          c(kz, iz, imon) = bebbe
          cA(kz, iz) = behbclam*bebbe
        end do
        ! Interior region: use backward propagator for next iteration
        do kz = 1, mxrho - ibl - 1
          cA(kz, iz) = cB(kz, iz)
        end do
      end do
!$omp end parallel do

      ! Apply symmetry to propagators at midplane
!$omp parallel do private(jz, kz)
      do iz = imitt + 1, imitt + ibl
        jz = imitt + 1 - (iz - imitt)
        do kz = 1, mxrho + kbl
          cA(kz, iz) = cA(kz, jz)
        end do
      end do
!$omp end parallel do

    end do

    if (ddmax .lt. CONV_TOL) exit  ! Converged

    ! Smart adaptive mixing: detect restart scenarios and adjust accordingly
    ! Fresh start: ddmax ~100 at iteration 2 → use adaptive mixing
    ! Restart: ddmax ~0.005 at iteration 2 → use conservative mixing
    if (niter == 2 .and. ddmax < 0.1d0) then
      ! Restart scenario detected: ddmax already small at iteration 2
      use_adaptive_mixing = .false.
      write (*, *) 'Restart detected (ddmax < 0.1 at iter 2): using conservative mixing'
    end if

    if (.not. use_adaptive_mixing) then
      ! Conservative mixing for restart scenarios
      dmm_adaptive = dmm
      dms_adaptive = dms
    else
      ! Detect oscillations: ddmax increased after decreasing
      if (niter .gt. 10 .and. ddmax .gt. ddmax_prev .and. ddmax .lt. 0.3d0) then
        oscillation_count = oscillation_count + 1
      else if (ddmax .lt. ddmax_prev) then
        ! Reset oscillation counter when making progress
        oscillation_count = max(0, oscillation_count - 1)
      end if

      ! Adaptive mixing for fresh starts: adjust based on convergence state
      ! Add conservative bias when oscillations detected
      if (ddmax .gt. 1.0d0) then
        dmm_adaptive = 0.90d0  ! Standard mixing when far from solution
        dms_adaptive = 0.50d0
      else if (ddmax .gt. 0.1d0) then
        dmm_adaptive = 0.85d0  ! Slightly more aggressive in mid-range
        dms_adaptive = 0.45d0
      else if (ddmax .gt. 0.01d0) then
        ! When approaching convergence, use conservative mixing if oscillating
        if (oscillation_count .gt. 3) then
          dmm_adaptive = 0.92d0  ! Very conservative to dampen oscillations
          dms_adaptive = 0.52d0
        else
          dmm_adaptive = 0.75d0  ! More aggressive when not oscillating
          dms_adaptive = 0.35d0
        end if
      else if (ddmax .gt. 0.001d0) then
        if (oscillation_count .gt. 3) then
          dmm_adaptive = 0.90d0  ! Very conservative near solution if oscillating
          dms_adaptive = 0.50d0
        else
          dmm_adaptive = 0.60d0  ! Aggressive near solution
          dms_adaptive = 0.25d0
        end if
      else
        if (oscillation_count .gt. 3) then
          dmm_adaptive = 0.88d0  ! Very conservative very close if oscillating
          dms_adaptive = 0.48d0
        else
          dmm_adaptive = 0.40d0  ! Extremely aggressive very close to solution
          dms_adaptive = 0.15d0
        end if
      end if

      if (niter .gt. 60) then
        if (oscillation_count .gt. 3) then
          write (*, '(A,E12.5,A,F5.3,A,I3)') 'Adaptive mixing (osc): ddmax=', ddmax, ', dmm=', &
                dmm_adaptive, ', osc_count=', oscillation_count
        else
          write (*, '(A,E12.5,A,F5.3)') 'Adaptive mixing: ddmax=', ddmax, ', dmm=', dmm_adaptive
        end if
      end if

      ! Update previous ddmax for next iteration
      ddmax_prev = ddmax
    end if

    tdmm = 1.d0 - dmm_adaptive
    tdms = 1.d0 - dms_adaptive

    ! Update densities using mixing scheme and check convergence
    ! Calculate new densities from propagators and mix with old values
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
          ! Calculate total monomer density from chain propagators
          ! Sum over all internal segments (convolution of forward and backward propagators)
          dumsum = 0.d0
          do k = 2, nmon - 1
            dumsum = c(j, i, k)*c(j, i, nmon + 1 - k) + dumsum
          end do
          tfem = 2.d0*c(j, i, 1)*ebelam(j, i)*dsqrt(edu(j, i))/ehbclam(j, i)
          tfdm = dumsum + tfem
          if (dabs(tfdm) .gt. 1.0d-14) then
            ddiff = abs(tfdm - fdmon(j, i))/tfdm
            if (ddiff .gt. ddmax) ddmax = ddiff
          end if
          fem(j, i) = fem(j, i)*dmm_adaptive + tdmm*tfem
          fdmon(j, i) = fdmon(j, i)*dmm_adaptive + tdmm*tfdm
        end if
      end do
    end do

    ! Apply symmetry to updated densities at midplane
    ! Mirror fdmon and fem across z = imitt to maintain symmetry
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

  ! ===== Output converged results =====
  ! Write density profiles and calculate thermodynamic properties
  rewind 89
  rewind 78
  rewind 85
  sumW = 0.d0
  z = -0.5d0*dz

  ! Write density profiles along z-axis at rho=0 (centerline)
  do iz = istp1, imitt
    z = z + dz
    ! File 85: monomer and end-segment densities at centerline
    write (85, *) z, fdmon(1, iz), fem(1, iz)
    ! File 89: propagators for segments 1,3,5,9 and ehbclam at centerline
    write (89, '(6f14.7)') z, c(1, iz, 1), c(1, iz, 3), c(1, iz, 5), &
      ehbclam(1, iz), c(1, iz, 9)

    ! Radially integrate density within radius 1.0 to get average
    fsum = 0.d0
    klm = nint(1.d0/drho)
    rho = -0.5d0*drho
    do i = 1, klm
      rho = rho + drho
      fsum = fsum + fdmon(i, iz)*2.d0*PI*rho
    end do
    ! File 78: z-position and radially averaged density
    write (78, *) z, fsum*drho/(PI*1.d0**2)
  end do

  ! Write radial profiles at z = zc1 (first colloid center position)
  rewind 83
  rewind 87
  iz = int(zc1*rdz) + 1
  rho = -0.5d0*drho
  do kr = 1, mxrho
    rho = rho + drho
    ! File 87: radial profile of propagators for segments 1,3,5,9 and ehbclam
    write (87, '(6f14.7)') rho, c(kr, iz, 1), c(kr, iz, 3), c(kr, iz, 5), &
      ehbclam(kr, iz), c(kr, iz, 9)
    ! File 83: radial profile of monomer and end-segment densities
    write (83, *) rho, fdmon(kr, iz), fem(kr, iz)
  end do

  ! ===== Calculate grand potential (thermodynamic potential) =====
  ! The grand potential Omega = F - mu*N measures the thermodynamic cost
  ! of the inhomogeneous density distribution relative to bulk
  bfde = 2.d0*bdpol  ! Bulk end-segment density
  bfdc = bdm - bfde  ! Bulk internal-segment density
  asumW = 0.d0
  sumsp = 0.d0
  sumsn = 0.d0
  bsumW = 0.d0
  chvol = 0.d0
  z = -0.5d0*dz

  ! Integrate grand potential density over system volume
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

      ! Only integrate outside colloid volume
      if (rsq .ge. Rcoll2) then
        ! Chemical potential contributions
        belamb = dlog(ebelam(kz, iz)) - emscale
        bclamb = 2.d0*(dlog(ehbclam(kz, iz)) - scalem)
        fde = fem(kz, iz)
        fdc = fdm - fde
        ! Excess free energy from hard-sphere interactions
        Fex = fdc*Y*(ae2(kz, iz) - ae1(kz, iz)) + 0.5d0*fde*ae2(kz, iz)

        ! Grand potential density omega(r) = f(r) - mu*rho(r)
        ! where f(r) is Helmholtz free energy density
        arsum = &
          rho*(fdc*bclamb + bfdc*bcmtrams + fde*belamb + bfde*bemtrams + &
               bdpol - fdm*rrnmon + Fex - bFex) + arsum
        brsum = &
          rho*(fdc*bclamb + fde*belamb - fdm*rrnmon + Fex - bFex) + brsum
      end if

      ! Add Lennard-Jones contribution to grand potential
      eexc = -dlog(edu(kz, iz))
      arsum = arsum - 0.5d0*rho*eexc*(fdm + bdm)
      brsum = brsum + rho*(0.5d0*(fdm - bdm)*eexc - fdm*eexc)
    end do
    ! Integrate radially: multiply by 2*pi*rho*drho
    asumW = 2.d0*PI*arsum*drho + asumW
    bsumW = 2.d0*PI*brsum*drho + bsumW
  end do
  ! Integrate along z-axis: multiply by dz
  asumW = asumW*dz
  bsumW = bsumW*dz
  ! Factor of 2 accounts for both halves of symmetric system
  aW = 2.d0*asumW
  bW = 2.d0*bsumW
  write (*, *)
  write (*, *) 'aW = ', aW
  write (*, *) 'bW = ', bW

  ! ===== Calculate forces on colloid from contact density =====
  ! Integrate contact density over colloid surface to get net force

  ! Determine integration limits for first colloid (centered at zc1)
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

  ! Calculate force on outer hemisphere (z < zc1) of first colloid
  rhoFo = 0.d0
  rcliffFo = 0.d0
  z = zmin - dz
  do iz = izmin, izc1 - 1
    z = z + dz
    zsq = (z - zc1)**2
    ! Only process z-slices that intersect the colloid
    if (zsq .le. Rcoll2) then
      rho = -0.5d0*drho
      irho = 0
      ! Find first grid point outside colloid at this z
      do
        rho = rho + drho
        irho = irho + 1
        if ((rho*rho + zsq) .gt. Rcoll2) exit
      end do
      Rc = dsqrt(rho*rho + zsq)
      rhoc = dsqrt(Rcoll2 - zsq)

      ! Quadratic interpolation to get density at exact colloid surface
      ! Use 3 points near boundary (irho, irho+1, irho+2)
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

      ! Lagrange interpolation to get density at colloid surface
      x = rhoc
      fdc = y1*(x - x2)*(x - x3)/((x1 - x2)*(x1 - x3)) + &
            y2*(x - x1)*(x - x3)/((x2 - x1)*(x2 - x3)) + &
            y3*(x - x1)*(x - x2)/((x3 - x1)*(x3 - x2))
      ! cos(theta) = (z - zc1)/Rcoll for surface normal direction
      ctheta = (z - zc1)/Rcoll
      rhoFo = 2.d0*PI*rhoc*ctheta*fdc + rhoFo
      rcliffFo = 2.d0*PI*ctheta*fdc + rcliffFo
      ict = ict + 1
      ctvec(ict) = ctheta
      cdens(ict) = fdc
    end if
  end do
  write (*, *) 'rcliffFo = ', Rcoll*rcliffFo*dz
  write (*, *) 'z = ', z

  ! Calculate force on inner hemisphere (z > zc1) of first colloid
  rhoFi = 0.d0
  rcliffFi = 0.d0
  z = zc1 - 0.5d0*dz
  do iz = izc1, izmax
    z = z + dz
    zsq = (z - zc1)**2
    ! Only process z-slices that intersect the colloid
    if (zsq .le. Rcoll2) then
      rho = -0.5d0*drho
      irho = 0
      ! Find first grid point outside colloid at this z
      do
        rho = rho + drho
        irho = irho + 1
        if ((rho*rho + zsq) .gt. Rcoll2) exit
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
      rhoFi = 2.d0*PI*rhoc*ctheta*fdc + rhoFi
      rcliffFi = 2.d0*PI*ctheta*fdc + rcliffFi
      ict = ict + 1
      ctvec(ict) = ctheta
      cdens(ict) = fdc
    end if
  end do
  write (*, *) 'rcliffFi = ', Rcoll*rcliffFi*dz
  rhoF = (rhoFi + rhoFo)*dz
  write (*, *)
  write (*, *) 'rcliffF = ', Rcoll*(rcliffFi + rcliffFo)*dz
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
  ctF = 2.d0*PI*Rcoll2*ctF
  ch2 = 2.d0*PI*Rcoll2*ch2
  write (*, *)
  write (*, *) 'ctF = ', ctF
  write (*, *)
  write (*, *) 'ch2 = ', ch2
  write (*, *)

  ! Alternative force calculation: integrate over rho slices at constant z

  ! Determine maximum radial index inside sphere
  irhomax = nint(Rcoll*rdrho + 1.d0)
  rhomax = (dfloat(irhomax) - 0.5d0)*drho
  ! Stay inside the sphere to avoid boundary issues
  irhomax = irhomax - 1
  rhomax = rhomax - drho
  write (*, *) 'rhomax,irhomax = ', rhomax, irhomax
  ict = 0
  zFo = 0.d0
  cho = 0.d0
  cliffFo = 0.d0
  rho = -0.5d0*dz

  ! Loop over radial slices from center outward (outer hemisphere in z)
  do irho = 1, irhomax
    rho = rho + drho
    rhosq = rho*rho
    ! Only process rho values that intersect the colloid
    if (rhosq .le. Rcoll2) then
      z = zc1 + 0.5d0*dz
      iz = izc1
      ! Find first grid point outside colloid at this rho (moving down in z)
      do
        z = z - dz
        iz = iz - 1
        zsq = (z - zc1)**2
        if ((rhosq + zsq) .gt. Rcoll2) exit
      end do
      Rc = dsqrt(rhosq + zsq)
      deltazc = dsqrt(Rcoll2 - rhosq)

      ! Quadratic interpolation in z-direction to get surface density
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

  ! Loop over radial slices in reverse (inner hemisphere in z)
  rho = rho + drho
  irho = irhomax + 1
  zFi = 0.d0
  chi = 0.d0
  cliffFi = 0.d0
  do krho = 1, irhomax
    irho = irho - 1
    rho = rho - drho
    rhosq = rho*rho
    ! Only process rho values that intersect the colloid
    if (rhosq .le. Rcoll2) then
      z = zc1 - 0.5d0*dz
      iz = izc1 - 1
      ! Find first grid point outside colloid at this rho (moving up in z)
      do
        z = z + dz
        iz = iz + 1
        zsq = (z - zc1)**2
        if ((rhosq + zsq) .gt. Rcoll2) exit
      end do
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

  ctF = 2.d0*PI*Rcoll2*ctF
  ch2 = 2.d0*PI*Rcoll2*ch2
  cckoll = 2.d0*PI*cckoll*drho
  ckoll = 2.d0*PI*ckoll*Rcoll*drho
  ckk = 2.d0*PI*ckk*Rcoll*drho
  ccc = 2.d0*PI*ccc
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

  ! Close files to ensure buffers are flushed
  close (ifc)
  close (ins)
  close (iep)

  ! Deallocate arrays before exit (COMMON block arrays are static)
  deallocate (c, cA, cB)

  STOP
END

! ============================================================================
! SUBROUTINE: CDFACT
! ============================================================================
! Calculates the normalization constant (cdnorm) for the contact density
! functional. This is computed via a three-dimensional integral over the
! hard sphere volume using trapezoidal integration in cylindrical coordinates.
! ============================================================================
subroutine CDFACT
  implicit none
  include 't2.inc.f90'

  ! Local variables
  integer :: iz, jz, iphi, irho, krhopmax, krhop
  double precision :: strho0, rho0, z, zpst, sume, zp, delz2, sumrhop
  double precision :: rhopmax, rho, rho02, rhop, rhomax2, fphi, phisum, phi, rho2, fact, tcd
  strho0 = 0.5d0*drho
  rho0 = strho0
  iz = 2*ism
  z = 2.d0*dhs - dz
  zpst = z - dhs - dz
  sume = 0.d0
  zp = zpst
  ! Integrate over sphere of diameter dhs centered at test point
  ! Triple integration: z', rho', phi in cylindrical coordinates
  do jz = iz - ism, iz + ism
    zp = zp + dz
    delz2 = (zp - z)**2
    sumrhop = 0.d0
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
!$omp simd reduction(+:phisum)
      do iphi = 1, nphi
!     Plus or minus sign doesn't matter for the value of the integral
        rho2 = rhomax2 - fphi*cos_phi(iphi)
        rho = dsqrt(rho2)
        irho = int(rho*rdrho) + 1
        phisum = 1.d0 + phisum
      end do
!$omp end simd
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

! ============================================================================
! SUBROUTINE: CDCALC
! ============================================================================
! Calculates the contact density (cdmonm) at each grid point. The contact
! density is a weighted average of the monomer density around a sphere of
! diameter dhs, computed using 3D integration in cylindrical coordinates
! with angular averaging.
! ============================================================================
subroutine CDCALC
  implicit none
  include 't2.inc.f90'

  ! Local variables
  integer :: iz, jz, kz, iphi, irho, krhop, krhopmax
  double precision :: z, zpst, rho0, sume, zp, delz2, sumrhop, rhopmax, rho02
  double precision :: rhop, rhomax2, fphi, phisum, phi, rho2, rho, fact

  ! Loop over all grid points to calculate contact density
!$omp parallel do private(z, zpst, rho0, kz, sume, zp, jz, delz2, sumrhop, rhopmax, krhopmax, rho02, rhop, rhomax2, fphi, phisum, iphi, rho2, rho, irho, fact) schedule(static)
  do iz = istp1 + ism, imitt
    z = dhs - 0.5d0*dz + dble(iz - (istp1 + ism) + 1)*dz
    zpst = z - dhs - dz
    rho0 = -0.5d0*drho
    do kz = 1, mxrho - ksm
      rho0 = rho0 + drho
      sume = 0.d0
      zp = zpst
      ! Integrate density over hard sphere volume
      do jz = iz - ism, iz + ism
        zp = zp + dz
        delz2 = (zp - z)**2
        sumrhop = 0.d0
        rhopmax = dsqrt(dabs(dhs2 - delz2))
        krhopmax = nint(rhopmax*rdrho)
        rho02 = rho0**2
        rhop = -0.5d0*drho
        do krhop = 1, krhopmax
          rhop = rhop + drho
          rhomax2 = rho02 + rhop*rhop
          fphi = 2.d0*rho0*rhop
          phisum = 0.d0
!$omp simd reduction(+:phisum)
          do iphi = 1, nphi
!     Plus or minus sign doesn't matter for the value of the integral
            rho2 = rhomax2 - fphi*cos_phi(iphi)
            rho = dsqrt(rho2)
            irho = int(rho*rdrho) + 1
            phisum = fdmon(irho, jz) + phisum
          end do
!$omp end simd
          sumrhop = rhop*phisum*dphi + sumrhop
        end do
        fact = 1.d0
        if (iabs(jz - iz) .eq. ism) fact = 0.5d0
        sume = 2.d0*sumrhop*drho*fact + sume
      end do
      cdmonm(kz, iz) = 3.d0*sume*dzrfp*cdnorm*rdhs3
    end do
  end do
!$omp end parallel do
  return
end

! ============================================================================
! SUBROUTINE: AVEC
! ============================================================================
! Calculates the excess free energy arrays (ae1, ae2) and the convolution
! term (convp) from the contact density using the Carnahan-Starling equation
! of state. These quantities are used in the density functional expressions
! for the polymer chain propagators.
! ============================================================================
subroutine AVEC
  implicit none
  include 't2.inc.f90'

  ! Local variables
  integer :: iz, kz, jz
  double precision :: cdt, pcdt, xsi, rxsi, sqrxsi, flog, daex1, daex2

  ! Calculate excess free energy from contact density using Carnahan-Starling EOS
!$omp parallel do private(kz, cdt, pcdt, xsi, rxsi, sqrxsi, flog, daex1, daex2) schedule(static)
  do iz = istp1 + 2*ism, imitt
    do kz = 1, mxrho - kbl
      cdt = cdmonm(kz, iz)*dhs3
      pcdt = PIS*cdt
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
    end do
  end do
!$omp end parallel do

!$omp parallel do private(jz, iz) schedule(static)
  do kz = 1, mxrho - kbl
    jz = imitt + 1
    do iz = imitt + 1, imitt + ibl
      jz = jz - 1
      ae1(kz, iz) = ae1(kz, jz)
      ae2(kz, iz) = ae2(kz, jz)
      convp(kz, iz) = convp(kz, jz)
    end do
  end do
!$omp end parallel do
  return
end

! ============================================================================
! SUBROUTINE: EBLMNEW
! ============================================================================
! Calculates the Boltzmann weight factors for polymer end segments (ebelam)
! and middle/hinge segments (ehbclam) at each grid point. These factors
! include contributions from the excess free energy and the convolution
! integrals. Excludes regions inside colloids.
! ============================================================================
subroutine EBLMNEW
  implicit none
  include 't2.inc.f90'

  ! Local variables
  integer :: iz, kz, jstart, irho0min, krhop, krhopmax, jz, iphi, irho
  double precision :: z, zpst, diffz2, strho0, rho0, rho02, rt2, sume, zp
  double precision :: delz2, zpcsq, zpc2sq, sumrhop, rhopmax, rhop, rhomax2
  double precision :: fphi, phisum, phi, rho2, rsq, rho, fact, trams, emtrams, cmtrams

  ! Set bulk values at boundaries
!$omp parallel do private(kz)
  do iz = 1, ibl
    do kz = 1, mxrho + kbl
      ebelam(kz, iz) = bebelam
      ehbclam(kz, iz) = behbclam
    end do
  end do
!$omp end parallel do

  ! Calculate Boltzmann factors from convolution integrals
!$omp parallel do private(z, jstart, zpst, diffz2, irho0min, strho0, rho0, kz, rho02, rt2, sume, zp, jz, delz2, zpcsq, zpc2sq, sumrhop, rhopmax, krhopmax, rhop, rhomax2, fphi, phisum, iphi, rho2, rsq, rho, irho, fact, trams, emtrams, cmtrams) schedule(static)
  do iz = ibl + 1, imitt
    z = -0.5d0*dz + dble(iz)*dz
    jstart = iz - ism
    zpst = z - dhs - dz
    diffz2 = (zc1 - z)**2
    irho0min = 1
    strho0 = -0.5d0*drho

    rho0 = strho0
    ! Loop over radial positions to calculate convolution integrals
    do kz = irho0min, mxrho
      rho0 = rho0 + drho
      rho02 = rho0**2

      ! Skip points inside first colloid
      rt2 = rho02 + (z - zc1)**2
      if (rt2 .lt. Rcoll2) then
        ehbclam(kz, iz) = 0.d0
        ebelam(kz, iz) = 0.d0
        cycle
      end if
      ! Skip points inside second colloid
      rt2 = rho02 + (z - zc2)**2
      if (rt2 .lt. Rcoll2) then
        ehbclam(kz, iz) = 0.d0
        ebelam(kz, iz) = 0.d0
        cycle
      end if

      ! Compute 3D convolution integral over hard sphere volume
      ! This gives the free energy contribution from local density variations
      sume = 0.d0
      zp = zpst
      ! Loop over z' within hard sphere diameter from current point
      do jz = jstart, iz + ism
        zp = zp + dz
        delz2 = (zp - z)**2
        zpcsq = (zp - zc1)**2
        zpc2sq = (zp - zc2)**2

        sumrhop = 0.d0
        ! Maximum rho' at this z' to stay within sphere of diameter dhs
        rhopmax = dsqrt(dabs(dhs2 - delz2))
        krhopmax = nint(rhopmax*rdrho)
        rho02 = rho0**2
        rhop = -0.5d0*drho
        ! Loop over rho' (radial coordinate at integration point)
        do krhop = 1, krhopmax
          rhop = rhop + drho
          rhomax2 = rho02 + rhop*rhop
          fphi = 2.d0*rho0*rhop
          phisum = 0.d0
          ! Loop over phi (angle between rho0 and rhop vectors)
          ! This completes the cylindrical coordinate integration
!$omp simd reduction(+:phisum)
          do iphi = 1, nphi
            ! Plus or minus sign doesn't matter for the value of the integral
            ! Calculate rho at integration point using law of cosines
            rho2 = rhomax2 - fphi*cos_phi(iphi)
            ! Skip if integration point is inside first colloid
            rsq = rho2 + zpcsq
            if (rsq .lt. Rcoll2) cycle
            ! Skip if integration point is inside second colloid
            rsq = rho2 + zpc2sq
            if (rsq .lt. Rcoll2) cycle
            rho = dsqrt(rho2)
            irho = int(rho*rdrho) + 1
            ! Sum convolution term (weighted density) over angular direction
            phisum = convp(irho, jz) + phisum
          end do
!$omp end simd
          ! Integrate over rho': weight by rhop*dphi*drho (cylindrical volume element)
          sumrhop = rhop*phisum*dphi + sumrhop
        end do
        ! Trapezoidal rule: half weight at boundaries
        fact = 1.d0
        if (iabs(jz - iz) .eq. ism) fact = 0.5d0
        ! Integrate over z': weight by 2*drho*fact*dz
        sume = 2.d0*sumrhop*drho*fact + sume
      end do
      ! Normalize convolution integral to get trams (local free energy contribution)
      trams = 3.d0*sume*dzrfp*cdnorm*rdhs3

      ! Calculate Boltzmann factors from free energy contributions
      ! emtrams: end-segment contribution (half of middle segment)
      emtrams = trams + 0.5d0*ae2(kz, iz)
      ! cmtrams: middle-segment contribution with Y factor for chain connectivity
      cmtrams = trams + Y*(ae2(kz, iz) - ae1(kz, iz))
      ! exp(-beta*F): Boltzmann weights with bulk chemical potential offset
      ebelam(kz, iz) = dexp(-emtrams + emscale)
      ehbclam(kz, iz) = dexp(-0.5d0*(cmtrams) + scalem)
    end do
  end do
!$omp end parallel do
  return
end

! ============================================================================
! SUBROUTINE: EBDU
! ============================================================================
! Calculates the external potential contribution (edu) from Lennard-Jones
! interactions with the polymer solution. Computes exp(-U_LJ) where U_LJ is
! the interaction energy between a test particle at (rho,z) and the entire
! density field, using symmetry to include both colloids.
! ============================================================================
! Helper subroutine for EBDU: compute rho integration with vectorization
function compute_rho_integral(fdmon_col, bdm, hvec_slice, mxrho) result(sumrho)
  implicit none
  integer, intent(in) :: mxrho
  double precision, intent(in) :: fdmon_col(0:*), hvec_slice(0:*), bdm
  double precision :: sumrho
  integer :: kprho

  sumrho = 0.d0
  do kprho = 1, mxrho
    sumrho = sumrho + (fdmon_col(kprho) - bdm)*hvec_slice(kprho)
  end do
end function compute_rho_integral

! ============================================================================
subroutine EBDU
  implicit none
  include 't2.inc.f90'

  ! Local variables
  integer :: iz, kz, krho, kprho, ipz, itdz
  double precision :: z, sumpint, tz, tdz, sumrho

  ! Set boundary values to unity (no external potential at boundaries)
!$omp parallel do private(kz)
  do iz = 1, ibl
    do kz = 1, mxrho + kbl
      edu(kz, iz) = 1.d0
    end do
  end do
!$omp end parallel do

  ! Calculate Lennard-Jones potential energy at each grid point
!$omp parallel do private(z, krho, sumpint, tz, tdz, ipz, itdz, sumrho, kprho) schedule(static)
  do iz = ibl + 1, imitt
    z = (dble(iz) - 0.5d0)*dz
    do krho = 1, mxrho
      sumpint = 0.d0
      tz = bl - 0.5d0*dz
      ! Integrate LJ interaction with left half of density distribution
      do ipz = ibl + 1, imitt
        tz = tz + dz
        tdz = z - tz
        itdz = nint(dabs(tdz*rdz))
        sumrho = 0.d0
        do kprho = 1, mxrho
          sumrho = sumrho + (fdmon(kprho, ipz) - bdm)*hvec(kprho, krho, itdz)
        end do
        sumpint = 2.d0*sumrho*drho + sumpint
      end do

      ! Integrate LJ interaction with right half (use symmetry)
      do ipz = imitt + 1, nfack - ibl
        tz = tz + dz
        tdz = z - tz
        itdz = nint(dabs(tdz*rdz))
        sumrho = 0.d0
        do kprho = 1, mxrho
          sumrho = sumrho + (fdmon(kprho, nfack + 1 - ipz) - bdm)*hvec(kprho, krho, itdz)
        end do
        sumpint = 2.d0*sumrho*drho + sumpint
      end do
      sumpint = sumpint*dz
      edu(krho, iz) = dexp(-sumpint)
    end do
  end do
!$omp end parallel do
  return
end

