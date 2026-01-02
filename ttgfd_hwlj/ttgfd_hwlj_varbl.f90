! ============================================================================
! TTGFD_HWLJ - Generalized Flory-Dimer Theory for Polymer Solutions
! ============================================================================
! Calculates polymer density profiles near hard-wall surfaces with
! Lennard-Jones interactions using density functional theory.
! Implements the Generalized Flory-Dimer (GFD) approximation for
! inhomogeneous polymer solutions in cylindrical geometry.
! ============================================================================
program platem
  use iso_fortran_env, only: real64, int32
  use polymer_dft_data
  implicit none

  ! Constants
  integer(int32), parameter :: MAXMON = 151

  ! Additional arrays for main program - also dynamically allocated
  real(real64), allocatable :: c(:, :, :), cA(:, :), cB(:, :)
  real(real64) :: cdens(0:1000), ctvec(0:1000)

  ! Integer variables
  integer(int32) :: i, j, k, iz, jz, kz, irho, krho, kprho, iphi, imon, kmon
  integer(int32) :: itdz, ict, kct, irdz, izc1, izmin, izmax, irho0min, irhomax, jstart
  integer(int32) :: ifc, ins, iep, niter, kcm, klm, kr, npphi

  ! Real variables
  real(real64) :: add, aeta, aex1, aex2, alj, arsum, asumw, aw
  real(real64) :: baex1, baex2, bclamb, bcmtrams, bconvp, bdaex1, bdaex2
  real(real64) :: bdpol, bds, bdt, bebbe, belamb, bemtrams, bfdc, bfde, bfex
  real(real64) :: bordekoll, brsum, bsumw, bw
  real(real64) :: ccc, ccckoll, cckoll, cdt, ch2, chempp, chi, cho, chvol
  real(real64) :: ckk, ckoll, clifffi, clifffo, cmtrams, ct, ctf
  real(real64) :: ctheta, ctn, ctp, cv
  real(real64) :: daex1, daex2, ddiff, ddmax, deltazc, delz2, diffz2
  real(real64) :: dmm_adaptive, dms_adaptive, dumsum
  real(real64) :: eexc, efact, emtrams

  ! Logical variables
  logical :: use_adaptive_mixing
  integer(int32) :: oscillation_count
  real(real64) :: ddmax_prev
  real(real64) :: fact, fdc, fdcm1, fdcn, fdcp1, fde, fdm, fex, ffact, fk, flog, fphi, fsum
  real(real64) :: pb, pcdt, pdasum, phi, phisum, pint
  real(real64) :: rc, rclifffi, rclifffo, rcyl2
  real(real64) :: rho, rho0, rho02, rho2, rhoc, rhof, rhofi, rhofo, rhomax, rhon, rhosq, rhoz2, rlj
  real(real64) :: rsq, rsq1, rsq2, rt2, rxsi, rxsib, rxsibsq, s2, sqrxsi, strho0, valid
  real(real64) :: sume, sumsn, sumsp, sumw
  real(real64) :: t, t1, t2, tdmm, tdms, tdz, tdzsq, tfdm, tfem, th, tn, trams
  real(real64) :: trho, trhosq, trmix, use1, useful
  real(real64) :: x, x1, x2, x3, xsi, xsib, y1, y2, y3
  real(real64) :: z, z2, z22, zfact, zfi, zfo, zmax, zmin, zp, zpc2sq, zpcsq, zpst, zsq

  ! Local temporary variables for bulk calculations (before structs are initialized)
  real(real64) :: dhs2, dhs3, rdhs3, rnmon, rrnmon, Yfact, scalem, emscale

  ! ========================================================================
  ! ========================================================================
  ! File unit numbers (automatically assigned by runtime)

  ! Open input and output files
  ! fcdfil: unknown allows both read and write (for restart capability)
  open (newunit=ifc, file='fcdfil', form='formatted', status='unknown')
  open (newunit=ins, file='input.tsph', form='formatted', status='old')
  open (newunit=iep, file='epfil', form='formatted', status='old')
  rewind ifc
  rewind ins
  rewind iep

  ! Read simulation parameters from input file into structured input type
  read (ins, *) input%bdm         ! Monomer bulk density
  read (ins, *) input%bdtot       ! Total bulk density
  bds = input%bdtot - input%bdm   ! Solvent bulk density
  read (ins, *) input%nmon        ! Number of monomers per polymer chain
  read (ins, *) input%dz          ! Grid spacing in z direction
  read (ins, *) input%drho        ! Grid spacing in radial direction
  read (ins, *) input%dphi        ! Angular grid spacing (input in units of pi)
  input%dphi = PI*input%dphi      ! Convert to radians
  read (ins, *) input%Rcoll       ! Colloid radius
  read (ins, *) input%zc1         ! Position of first colloid center
  read (ins, *) input%collsep     ! Separation between colloid centers
  read (ins, *) input%Rcyl        ! Cylinder radius (system boundary)
  read (ins, *) input%ioimaxm     ! Maximum number of iterations
  read (ins, *) input%dmm, input%dms  ! Density mixing parameters (monomer, solvent)
  read (ins, *) input%kread       ! Read initial guess from file (0=no, 1=yes)
  read (ins, *) input%bl          ! Bond length
  read (ins, *) input%dhs         ! Hard sphere diameter (monomer)
  read (ins, *) input%dpphi       ! Angular grid spacing for potential calculation

  ! Initialize grid parameters from input
  call initialize_grid_params(input, grid)

  ! Compute local derived parameters needed for bulk thermodynamic calculations
  ! These will be recomputed and stored in structs by initialize_computed_params later
  dhs2 = input%dhs*input%dhs
  dhs3 = dhs2*input%dhs
  rdhs3 = 1.d0/dhs3
  rnmon = dble(input%nmon)
  rrnmon = 1.d0/rnmon

  ! Read Lennard-Jones energy parameter and compute LJ coefficients
  read (iep, *) input%epslj
  alj = 4.d0*input%epslj*dhs3*dhs3    ! Attractive (r^-6) coefficient
  rlj = 4.d0*input%epslj*dhs3**4      ! Repulsive (r^-12) coefficient

  ! Compute additional local parameters needed for bulk calculations
  Rcyl2 = input%Rcyl*input%Rcyl
  tdmm = 1.d0 - input%dmm
  tdms = 1.d0 - input%dms

  ! Initialize cosine lookup table for input%dphi (using grid%nphi from initialize_grid_params)
  do iphi = 1, grid%nphi
    phi = (dble(iphi) - 0.5d0)*input%dphi
    cos_phi(iphi) = dcos(phi)
  end do

  ! Bulk thermodynamic properties
  Yfact = (rnmon - 2.d0)*Y
  bdpol = input%bdm/rnmon  ! Polymer bulk density

  ! Hard sphere packing fraction and related quantities
  bdt = input%bdm*dhs3
  aeta = PIS*bdt
  xsib = 1.d0 - aeta
  rxsib = 1.d0/xsib
  rxsibsq = rxsib*rxsib

  ! Excess free energy terms (Carnahan-Starling)
  aex1 = -(C1 + 1.d0)*dlog(xsib) - &
         0.5d0*(AA1*PIS*bdt + BB1*(PIS*bdt)**2)*rxsibsq
  aex2 = -(C2 + 1.d0)*dlog(xsib) - &
         0.5d0*(AA2*PIS*bdt + BB2*(PIS*bdt)**2)*rxsibsq
  bFex = (input%bdm - 2.d0*bdpol)*Y*(aex2 - aex1) + bdpol*aex2

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
  chempp = dlog(bdpol) + Yfact*(aex2 - aex1) + aex2 + PIS*input%bdm*pdasum*dhs3
  scalem = chempp/(2.d0*rnmon)
  emscale = 2.d0*scalem

  ! Convolution terms for inhomogeneous density functional
  bconvp = (Y*(input%bdm - 2.d0*input%bdm*rrnmon)*(daex2 - daex1) + &
            input%bdm*rrnmon*daex2)*PIS*dhs3
  trams = bconvp
  emtrams = trams + 0.5d0*aex2
  cmtrams = trams + Y*(aex2 - aex1)
  bemtrams = emtrams
  bcmtrams = cmtrams

  ! Initialize computed parameters from input, grid, and bulk calculations
  ! This recomputes and stores all derived parameters in structured form
  call initialize_computed_params(input, grid, computed, chempp, emtrams, cmtrams)

  ! Print simulation parameters
  write (*, *) 'GFD POLYMER SOLUTION MODEL!'
  write (*, *) 'input%bdm,bdpol =', input%bdm, bdpol
  write (*, *) 'monomer density  = ', input%bdm
  write (*, *) 'bdt = ', bdt
  write (*, *) 'collsep,input%dz = ', input%collsep, input%dz
  write (*, *) 'input%Rcoll = ', input%Rcoll
  write (*, *) 'bond length (input%bl): ', input%bl
  write (*, *) 'monomer hs diameter (input%bl): ', input%dhs
  write (*, *) 'no. of monomers/polymer = ', input%nmon
  write (*, *) 'max no. of iterations = ', input%ioimaxm
  write (*, *) 'polymer chemical pot. (betamu) = ', chempp
  write (*, *) 'solvent chemical pot. (betamu) = ', computed%chemps
  write (*, *) 'total bulk pressure = ', Pb
  write (*, *) 'input%zc1,computed%zc2 = ', input%zc1, computed%zc2
  write (*, *) 'nfack,grid%imitt = ', grid%nfack, grid%imitt
  write (*, *) 'istp1,grid%islut = ', grid%istp1, grid%islut
  write (*, *) 'istp1s,isluts = ', grid%istp1s, grid%isluts
  write (*, *) 'ism,grid%ibl = ', grid%ism, grid%ibl
  write (*, *) 'ksm,grid%kbl = ', grid%ksm, grid%kbl
  write (*, *) 'dmm,dms (density mixing param. mon.,solv.) = ', input%dmm, input%dms
  write (*, *) 'Rcyl  = ', input%Rcyl
  write (*, *) 'bFex = ', bFex
  write (*, *) 'bebelam,computed%behbclam = ', computed%bebelam, computed%behbclam

  ! Allocate module arrays based on calculated grid dimensions
  ! Include extra space for boundary cells (grid%kbl, grid%ibl)
  ! hvec needs nfack-1 for z-dimension (used in LJ potential table)
  call allocate_arrays(grid%mxrho + grid%kbl, grid%imitt + grid%ibl, grid%nfack - 1)

  ! Allocate main program arrays
  ! Note: cB needs grid%nfack dimension because it's accessed as grid%islut + 1 - iz
  allocate (c(0:grid%mxrho + grid%kbl, 0:grid%imitt + grid%ibl, MAXMON))
  allocate (cA(0:grid%mxrho + grid%kbl, 0:grid%imitt + grid%ibl))
  allocate (cB(0:grid%mxrho + grid%kbl, 0:grid%nfack))

  ! Calculate normalization constant for contact density
  call CDFACT(input, grid, computed, cos_phi, computed%cdnorm)
  write (*, *) 'cdnorm = ', computed%cdnorm

  ! Initialize excess free energy arrays near z-boundaries (left side)
  ! These regions are near the system edge and require special treatment
  do iz = grid%istp1, grid%istp1 + 2*grid%ibl - 1
  do kz = 1, grid%mxrho + grid%kbl
    cdt = input%bdm*dhs3
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
    convp(kz, iz) = (Y*(input%bdm - 2.d0*input%bdm*rrnmon)*(daex2 - daex1) + &
                     input%bdm*rrnmon*daex2)*PIS*dhs3
  end do
  end do
  ! Initialize excess free energy arrays near radial boundaries (outer edge)
  do iz = grid%istp1 + 2*grid%ibl, grid%imitt + grid%ibl
  do kz = grid%mxrho - grid%kbl + 1, grid%mxrho + grid%kbl
    cdt = input%bdm*dhs3
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
    convp(kz, iz) = (Y*(input%bdm - 2.d0*input%bdm*rrnmon)*(daex2 - daex1) + &
                     input%bdm*rrnmon*daex2)*PIS*dhs3
  end do
  end do

  ! Initialize density fields: either from scratch or read from file
  if (input%kread .eq. 0) then
    ! Starting from bulk values with excluded volume for colloids
    z = -0.5d0*input%dz
    ! Skip boundary region at z < 0
    do iz = 1, grid%ibl
      z = z + input%dz
    end do
    ! Initialize all grid points to bulk values, then zero out colloid interiors
    do iz = grid%ibl + 1, grid%imitt
      z = z + input%dz
      z2 = (z - input%zc1)**2
      z22 = (z - computed%zc2)**2
      rho = -0.5d0*input%drho
      ! Loop over radial positions
      do kz = 1, grid%mxrho
        rho = rho + input%drho
        fdmon(kz, iz) = input%bdm
        fem(kz, iz) = 2.d0*fdmon(kz, iz)*rrnmon
        ebelam(kz, iz) = computed%bebelam
        ehbclam(kz, iz) = computed%behbclam
        ! Check if point is inside first colloid
        rt2 = rho*rho + z2
        if (rt2 .lt. computed%Rcoll2) then
          fdmon(kz, iz) = 0.d0
          fem(kz, iz) = 0.d0
          ebelam(kz, iz) = 0.d0
          ehbclam(kz, iz) = 0.d0
        end if
        ! Check if point is inside second colloid
        rt2 = rho*rho + z22
        if (rt2 .lt. computed%Rcoll2) then
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
    do iz = grid%istp1, grid%imitt
    do kz = 1, grid%mxrho
      read (ifc, *) t1, t2, fdmon(kz, iz), fem(kz, iz)
    end do
    end do
  end if

  ! Set boundary conditions at z-boundaries (left edge)
  ! All densities set to bulk values
  do iz = grid%istp1, grid%ibl
  do kz = 1, grid%mxrho + grid%kbl
    fdmon(kz, iz) = input%bdm
    fem(kz, iz) = 2.d0*fdmon(kz, iz)*rrnmon
    ebelam(kz, iz) = computed%bebelam
    ehbclam(kz, iz) = computed%behbclam
    cdmonm(kz, iz) = input%bdm
  end do
  end do
  ! Set boundary conditions at radial edge (outer cylinder boundary)
  ! All densities set to bulk values
  do iz = 1, grid%imitt
  do kz = grid%mxrho + 1, grid%mxrho + grid%kbl
    fdmon(kz, iz) = input%bdm
    fem(kz, iz) = 2.d0*fdmon(kz, iz)*rrnmon
    ebelam(kz, iz) = computed%bebelam
    ehbclam(kz, iz) = computed%behbclam
    cdmonm(kz, iz) = input%bdm
  end do
  end do

  ! Apply symmetry boundary conditions at z = grid%imitt midplane
  jz = grid%imitt + 1
  do iz = grid%imitt + 1, grid%imitt + grid%ibl
    jz = jz - 1
    do kz = 1, grid%mxrho + grid%kbl
      fdmon(kz, iz) = fdmon(kz, jz)
      fem(kz, iz) = fem(kz, jz)
      ebelam(kz, iz) = ebelam(kz, jz)
      ehbclam(kz, iz) = ehbclam(kz, jz)
      cdmonm(kz, iz) = input%bdm
    end do
  end do
  write (*, *) 'fdmon(1,1) = ', fdmon(1, 1)
  write (*, *) 'fdmon(1,11) = ', fdmon(1, 11)

  ! Precompute Lennard-Jones interaction potential on grid (hvec array)
  ! This tabulates U_LJ for all distance combinations to speed up later calculations
  write (*, *) 'dpphi = ', input%dpphi
  input%dpphi = input%dpphi*PI
  npphi = int(PI/input%dpphi + 0.01d0)
  write (*, *) 'dpphi,npphi = ', input%dpphi, npphi

  ! Initialize cosine lookup table for dpphi
  do iphi = 1, npphi
    phi = (dble(iphi) - 0.5d0)*input%dpphi
    cos_pphi(iphi) = dcos(phi)
  end do

  kcm = grid%nfack
  ! Triple loop over rho, rho', z to compute pairwise LJ interaction integrals
  ! Loop order matches F77 for numerical consistency
!$omp parallel do private(tdz, tdzsq, rho, rhosq, use1, trho, trhosq, trmix, useful, pint, iphi, s2, krho, kprho) schedule(static)
  do itdz = 0, kcm - 1
    tdz = -input%dz + dble(itdz + 1)*input%dz
    tdzsq = tdz*tdz
    rho = -0.5d0*input%drho
    do krho = 1, grid%mxrho
      rho = rho + input%drho
      rhosq = rho*rho
      use1 = tdzsq + rhosq
      trho = -0.5d0*input%drho
      do kprho = 1, grid%mxrho
        trho = trho + input%drho
        trhosq = trho*trho
        trmix = 2.d0*trho*rho
        useful = use1 + trhosq
        pint = 0.d0
        ! Angular integration for cylindrical geometry
!$omp simd reduction(+:pint)
        do iphi = 1, npphi
          s2 = useful - trmix*cos_pphi(iphi)
          if (s2 .gt. computed%dhs2) then
            ! Lennard-Jones potential: U(r) = 4*epsilon*[(sigma/r)^12 - (sigma/r)^6]
            pint = rlj/s2**6 - alj/s2**3 + pint
          end if
        end do
!$omp end simd
        hvec(kprho, krho, itdz) = trho*pint*input%dpphi
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
    if (niter .gt. input%ioimaxm) then
      write (*, *) 'NITER.GT.IOIMAXM !', niter
      exit
    end if

    ! Update fields in density functional theory calculation
    call CDCALC(input, grid, computed, fdmon, cos_phi, cdmonm)     ! Calculate contact density
    call AVEC(grid, computed, cdmonm, fdmon, fem, ae1, ae2, convp)       ! Calculate excess free energy
    call EBLMNEW(input, grid, computed, convp, ae1, ae2, cos_phi, &
                 ebelam, ehbclam)    ! Calculate end-segment Boltzmann factors
    call EBDU(input, grid, computed, fdmon, hvec, edu)       ! Calculate external potential contribution

    ! Apply boundary conditions at outer radial edge (rho > Rcyl)
    ! Set to bulk values since density should approach bulk far from colloids
    do iz = grid%istp1, grid%imitt
    do kz = grid%mxrho + 1, grid%mxrho + grid%kbl
      ebelam(kz, iz) = computed%bebelam
      ehbclam(kz, iz) = computed%behbclam
      edu(kz, iz) = 1.d0
    end do
    end do

    ! Apply symmetry boundary conditions at z = grid%imitt (midplane between colloids)
    ! The system is symmetric about the midplane, so mirror the field values
    ! This loop copies from iz = grid%imitt down to iz = 1 (reverse order)
    jz = grid%imitt + 1
    do iz = grid%imitt + 1, grid%imitt + grid%ibl
      jz = jz - 1
      do kz = 1, grid%mxrho + grid%kbl
        ebelam(kz, iz) = ebelam(kz, jz)
        ehbclam(kz, iz) = ehbclam(kz, jz)
        edu(kz, iz) = edu(kz, jz)
      end do
    end do

    ! Second symmetry application (appears redundant but ensures consistency)
    jz = grid%imitt + 1
    do iz = grid%imitt + 1, grid%imitt + grid%ibl
      jz = jz - 1
      do kz = 1, grid%mxrho + grid%kbl
        ebelam(kz, iz) = ebelam(kz, jz)
        ehbclam(kz, iz) = ehbclam(kz, jz)
        edu(kz, iz) = edu(kz, jz)
      end do
    end do

    ! Calculate cA: the propagator for polymer end segments
    ! cA(r) = exp(-beta*mu_end)*exp(-U_LJ) where:
    !   ebelam = exp(-beta*mu_end) from hard-sphere and chain connectivity
    !   edu = exp(-U_LJ) from Lennard-Jones interactions
    do iz = grid%istp1, grid%imitt + grid%ibl
    do kz = 1, grid%mxrho + grid%kbl
      cA(kz, iz) = ebelam(kz, iz)*edu(kz, iz)
    end do
    end do

    ! Propagate polymer chains segment by segment from end to end
    ! This calculates c(r,i) = propagator for segment i at position r
    imon = input%nmon
    do kmon = 1, input%nmon - 1
      imon = imon - 1

      ! Merge all chain propagation loops into single parallel region
      ! to eliminate thread synchronization barriers
!$omp parallel private(iz, z, jstart, zpst, irho0min, strho0, rho0, kz, rho02, rt2, sume, zp, jz, delz2, zpcsq, zpc2sq, phisum, zfact, rhoz2, fphi, iphi, rho2, rsq1, rsq2, valid, rho, irho, fact, efact, ffact, bebbe)

      ! Loop over all spatial grid points
!$omp do schedule(static)
      do iz = grid%istp1 + grid%ibl, grid%imitt
        z = input%bl - 0.5d0*input%dz + dble(iz - (grid%istp1 + grid%ibl) + 1)*input%dz
        jstart = iz - grid%ibl
        zpst = z - input%bl - input%dz
        irho0min = 1
        strho0 = -0.5d0*input%drho
        rho0 = strho0
        do kz = irho0min, grid%mxrho - grid%ibl
          rho0 = rho0 + input%drho

          rho02 = rho0**2
          rt2 = rho02 + (z - input%zc1)**2
          if (rt2 .lt. computed%Rcoll2) then
            c(kz, iz, imon) = 0.d0
            if (iz .gt. grid%imitt - grid%ibl - 1) cB(kz, grid%islut + 1 - iz) = 0.d0
            cB(kz, iz) = 0.d0
            cycle
          end if
          rt2 = rho02 + (z - computed%zc2)**2
          if (rt2 .lt. computed%Rcoll2) then
            c(kz, iz, imon) = 0.d0
            cB(kz, grid%islut + 1 - iz) = 0.d0
            cB(kz, iz) = 0.d0
            cycle
          end if

          ! Integrate over bond orientations: sum contributions from all points
          ! within bond length input%bl of current position (rho0, z)
          sume = 0.d0
          zp = zpst
          do jz = jstart, iz + grid%ibl
            zp = zp + input%dz
            delz2 = (zp - z)**2
            zpcsq = (zp - input%zc1)**2
            zpc2sq = (zp - computed%zc2)**2
            phisum = 0.d0
            zfact = dabs(computed%bl2 - delz2)
            rhoz2 = rho0**2 + zfact
            fphi = 2.d0*rho0*dsqrt(zfact)
!$omp simd reduction(+:phisum)
            do iphi = 1, grid%nphi
!     Plus or minus sign doesn't matter for the value of the integral
              rho2 = rhoz2 - fphi*cos_phi(iphi)
              rsq1 = rho2 + zpcsq
              rsq2 = rho2 + zpc2sq

              ! Mask-based approach: 1.0 if outside both colloids, 0.0 if inside either
              ! This eliminates conditional exits (cycle) for full SIMD vectorization
              valid = merge(1.0d0, 0.0d0, rsq1 >= computed%Rcoll2 .and. rsq2 >= computed%Rcoll2)

              rho = dsqrt(rho2)
              irho = int(rho*computed%rdrho) + 1
              phisum = phisum + valid * cA(irho, jz)
            end do
!$omp end simd
            fact = 1.d0
            if (iabs(jz - iz) .eq. grid%ibl) fact = 0.5d0
            sume = 2.d0*phisum*input%dphi*fact + sume
          end do
          efact = dsqrt(edu(kz, iz))
          ffact = sume*computed%dzrfp*ehbclam(kz, iz)/input%bl*efact
          c(kz, iz, imon) = ffact
          if (iz .gt. grid%imitt - grid%ibl - 1) cB(kz, grid%islut + 1 - iz) = ffact*ehbclam(kz, iz)*efact
          cB(kz, iz) = ffact*ehbclam(kz, iz)*efact
        end do
      end do
!$omp end do
      ! Implicit barrier: Loop 3 reads cB(kz,iz) which Loop 1 writes

      ! Handle boundary regions: propagators at z-boundaries
!$omp do schedule(static)
      do iz = grid%istp1, grid%ibl
      do kz = 1, grid%mxrho + grid%kbl
        bebbe = computed%behbclam*cA(kz, iz)
        c(kz, iz, imon) = bebbe
        cA(kz, iz) = computed%behbclam*bebbe
      end do
      end do
!$omp end do nowait

      ! Handle radial boundaries and update propagators
!$omp do schedule(static)
      do iz = grid%ibl + 1, grid%imitt
        ! Outer radial boundary: use bulk propagators
        do kz = grid%mxrho - grid%kbl, grid%mxrho + grid%kbl
          bebbe = computed%behbclam*cA(kz, iz)
          c(kz, iz, imon) = bebbe
          cA(kz, iz) = computed%behbclam*bebbe
        end do
        ! Interior region: use backward propagator for next iteration
        do kz = 1, grid%mxrho - grid%ibl - 1
          cA(kz, iz) = cB(kz, iz)
        end do
      end do
!$omp end do
      ! Implicit barrier: Loop 4 reads cA which Loops 2 & 3 write

      ! Apply symmetry to propagators at midplane
!$omp do schedule(static)
      do iz = grid%imitt + 1, grid%imitt + grid%ibl
        jz = grid%imitt + 1 - (iz - grid%imitt)
        do kz = 1, grid%mxrho + grid%kbl
          cA(kz, iz) = cA(kz, jz)
        end do
      end do
!$omp end do

!$omp end parallel

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
      dmm_adaptive = input%dmm
      dms_adaptive = input%dms
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
          write (*, '(A,E12.5,A,F5.3,A,I3)') 'Adaptive mixing (osc): ddmax=', ddmax, ', input%dmm=', &
                dmm_adaptive, ', osc_count=', oscillation_count
        else
          write (*, '(A,E12.5,A,F5.3)') 'Adaptive mixing: ddmax=', ddmax, ', input%dmm=', dmm_adaptive
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
    z = -0.5d0*input%dz
    do i = grid%istp1, grid%imitt
      z = z + input%dz
      diffz2 = (z - input%zc1)**2
      rho = -0.5d0*input%drho
      do j = 1, grid%mxrho
        rho = rho + input%drho
        rsq = rho*rho + diffz2
        if (rsq .lt. computed%Rcoll2) then
          fem(j, i) = 0.d0
          fdmon(j, i) = 0.d0
        else
          ! Calculate total monomer density from chain propagators
          ! Sum over all internal segments (convolution of forward and backward propagators)
          dumsum = 0.d0
          do k = 2, input%nmon - 1
            dumsum = c(j, i, k)*c(j, i, input%nmon + 1 - k) + dumsum
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
    ! Mirror fdmon and fem across z = grid%imitt to maintain symmetry
    jz = grid%imitt + 1
    do iz = grid%imitt + 1, grid%imitt + grid%ibl
      jz = jz - 1
      do kz = 1, grid%mxrho + grid%kbl
        fdmon(kz, iz) = fdmon(kz, jz)
        fem(kz, iz) = fem(kz, jz)
      end do
    end do

  end do  ! End of main iteration loop

  ! Check if maximum iterations exceeded
  if (niter .gt. input%ioimaxm) then
    stop
  end if

  ! ===== Output converged results =====
  ! Write density profiles and calculate thermodynamic properties
  rewind 89
  rewind 78
  rewind 85
  sumW = 0.d0
  z = -0.5d0*input%dz

  ! Write density profiles along z-axis at rho=0 (centerline)
  do iz = grid%istp1, grid%imitt
    z = z + input%dz
    ! File 85: monomer and end-segment densities at centerline
    write (85, *) z, fdmon(1, iz), fem(1, iz)
    ! File 89: propagators for segments 1,3,5,9 and ehbclam at centerline
    write (89, '(6f14.7)') z, c(1, iz, 1), c(1, iz, 3), c(1, iz, 5), &
      ehbclam(1, iz), c(1, iz, 9)

    ! Radially integrate density within radius 1.0 to get average
    fsum = 0.d0
    klm = nint(1.d0/input%drho)
    rho = -0.5d0*input%drho
    do i = 1, klm
      rho = rho + input%drho
      fsum = fsum + fdmon(i, iz)*2.d0*PI*rho
    end do
    ! File 78: z-position and radially averaged density
    write (78, *) z, fsum*input%drho/(PI*1.d0**2)
  end do

  ! Write radial profiles at z = input%zc1 (first colloid center position)
  rewind 83
  rewind 87
  iz = int(input%zc1*computed%rdz) + 1
  rho = -0.5d0*input%drho
  do kr = 1, grid%mxrho
    rho = rho + input%drho
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
  bfdc = input%bdm - bfde  ! Bulk internal-segment density
  asumW = 0.d0
  sumsp = 0.d0
  sumsn = 0.d0
  bsumW = 0.d0
  chvol = 0.d0
  z = -0.5d0*input%dz

  ! Integrate grand potential density over system volume
  do iz = grid%istp1, grid%imitt
    z = z + input%dz
    arsum = 0.d0
    brsum = 0.d0
    cv = 0.d0
    diffz2 = (z - input%zc1)**2
    rho = -0.5d0*input%drho
    do kz = 1, grid%mxrho
      rho = rho + input%drho
      rsq = rho*rho + diffz2
      fdm = fdmon(kz, iz)

      ! Only integrate outside colloid volume
      if (rsq .ge. computed%Rcoll2) then
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
      arsum = arsum - 0.5d0*rho*eexc*(fdm + input%bdm)
      brsum = brsum + rho*(0.5d0*(fdm - input%bdm)*eexc - fdm*eexc)
    end do
    ! Integrate radially: multiply by 2*pi*rho*input%drho
    asumW = 2.d0*PI*arsum*input%drho + asumW
    bsumW = 2.d0*PI*brsum*input%drho + bsumW
  end do
  ! Integrate along z-axis: multiply by input%dz
  asumW = asumW*input%dz
  bsumW = bsumW*input%dz
  ! Factor of 2 accounts for both halves of symmetric system
  aW = 2.d0*asumW
  bW = 2.d0*bsumW
  write (*, *)
  write (*, *) 'aW = ', aW
  write (*, *) 'bW = ', bW

  ! ===== Calculate forces on colloid from contact density =====
  ! Integrate contact density over colloid surface to get net force

  ! Determine integration limits for first colloid (centered at input%zc1)
  izmin = nint((input%zc1 + 0.5d0*input%dz - input%Rcoll)*computed%rdz + 0.5d0)
  zmin = (dfloat(izmin) - 0.5d0)*input%dz
  izmax = nint((input%zc1 - 0.5d0*input%dz + input%Rcoll)*computed%rdz + 0.5d0)
  zmax = (dfloat(izmax) - 0.5d0)*input%dz
  izc1 = nint(input%zc1*computed%rdz + 0.5d0)
  write (*, *) 'zmin,input%zc1,zmax = ', zmin, input%zc1, zmax
  write (*, *) 'izmin,izc1,izmax = ', izmin, izc1, izmax
  write (*, *) dfloat(izmin)*input%dz - 0.5d0*input%dz, dfloat(izmax)*input%dz - 0.5d0*input%dz
  write (*, *) dfloat(izc1)*input%dz - 0.5d0*input%dz
  ict = 0

  ! Calculate force on outer hemisphere (z < input%zc1) of first colloid
  rhoFo = 0.d0
  rcliffFo = 0.d0
  z = zmin - input%dz
  do iz = izmin, izc1 - 1
    z = z + input%dz
    zsq = (z - input%zc1)**2
    ! Only process z-slices that intersect the colloid
    if (zsq .le. computed%Rcoll2) then
      rho = -0.5d0*input%drho
      irho = 0
      ! Find first grid point outside colloid at this z
      do
        rho = rho + input%drho
        irho = irho + 1
        if ((rho*rho + zsq) .gt. computed%Rcoll2) exit
      end do
      Rc = dsqrt(rho*rho + zsq)
      rhoc = dsqrt(computed%Rcoll2 - zsq)

      ! Quadratic interpolation to get density at exact colloid surface
      ! Use 3 points near boundary (irho, irho+1, irho+2)
      if (dabs(fdmon(irho, iz)) .gt. 0.00000001d0) then
        y3 = fdmon(irho, iz)
        y2 = fdmon(irho + 1, iz)
        y1 = fdmon(irho + 2, iz)
        x3 = rho
        x2 = rho + input%drho
        x1 = rho + 2.d0*input%drho
      else
        y3 = fdmon(irho + 1, iz)
        y2 = fdmon(irho + 2, iz)
        y1 = fdmon(irho + 3, iz)
        x3 = rho + input%drho
        x2 = rho + 2.d0*input%drho
        x1 = rho + 3.d0*input%drho
        write (*, *) 'TJOHO!'
      end if

      ! Lagrange interpolation to get density at colloid surface
      x = rhoc
      fdc = y1*(x - x2)*(x - x3)/((x1 - x2)*(x1 - x3)) + &
            y2*(x - x1)*(x - x3)/((x2 - x1)*(x2 - x3)) + &
            y3*(x - x1)*(x - x2)/((x3 - x1)*(x3 - x2))
      ! cos(theta) = (z - input%zc1)/input%Rcoll for surface normal direction
      ctheta = (z - input%zc1)/input%Rcoll
      rhoFo = 2.d0*PI*rhoc*ctheta*fdc + rhoFo
      rcliffFo = 2.d0*PI*ctheta*fdc + rcliffFo
      ict = ict + 1
      ctvec(ict) = ctheta
      cdens(ict) = fdc
    end if
  end do
  write (*, *) 'rcliffFo = ', input%Rcoll*rcliffFo*input%dz
  write (*, *) 'z = ', z

  ! Calculate force on inner hemisphere (z > input%zc1) of first colloid
  rhoFi = 0.d0
  rcliffFi = 0.d0
  z = input%zc1 - 0.5d0*input%dz
  do iz = izc1, izmax
    z = z + input%dz
    zsq = (z - input%zc1)**2
    ! Only process z-slices that intersect the colloid
    if (zsq .le. computed%Rcoll2) then
      rho = -0.5d0*input%drho
      irho = 0
      ! Find first grid point outside colloid at this z
      do
        rho = rho + input%drho
        irho = irho + 1
        if ((rho*rho + zsq) .gt. computed%Rcoll2) exit
      end do
      Rc = dsqrt(rho*rho + zsq)
      rhoc = dsqrt(computed%Rcoll2 - zsq)

      if (dabs(fdmon(irho, iz)) .gt. 0.00000001d0) then
        y3 = fdmon(irho, iz)
        y2 = fdmon(irho + 1, iz)
        y1 = fdmon(irho + 2, iz)
        x3 = rho
        x2 = rho + input%drho
        x1 = rho + 2.d0*input%drho
      else
        y3 = fdmon(irho + 1, iz)
        y2 = fdmon(irho + 2, iz)
        y1 = fdmon(irho + 3, iz)
        x3 = rho + input%drho
        x2 = rho + 2.d0*input%drho
        x1 = rho + 3.d0*input%drho
        write (*, *) 'TJOHO!!!!', fdmon(irho, iz), rho
      end if

      x = rhoc
      fdc = y1*(x - x2)*(x - x3)/((x1 - x2)*(x1 - x3)) + &
            y2*(x - x1)*(x - x3)/((x2 - x1)*(x2 - x3)) + &
            y3*(x - x1)*(x - x2)/((x3 - x1)*(x3 - x2))
      ctheta = (z - input%zc1)/input%Rcoll
      rhoFi = 2.d0*PI*rhoc*ctheta*fdc + rhoFi
      rcliffFi = 2.d0*PI*ctheta*fdc + rcliffFi
      ict = ict + 1
      ctvec(ict) = ctheta
      cdens(ict) = fdc
    end if
  end do
  write (*, *) 'rcliffFi = ', input%Rcoll*rcliffFi*input%dz
  rhoF = (rhoFi + rhoFo)*input%dz
  write (*, *)
  write (*, *) 'rcliffF = ', input%Rcoll*(rcliffFi + rcliffFo)*input%dz
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
  ctF = 2.d0*PI*computed%Rcoll2*ctF
  ch2 = 2.d0*PI*computed%Rcoll2*ch2
  write (*, *)
  write (*, *) 'ctF = ', ctF
  write (*, *)
  write (*, *) 'ch2 = ', ch2
  write (*, *)

  ! Alternative force calculation: integrate over rho slices at constant z

  ! Determine maximum radial index inside sphere
  irhomax = nint(input%Rcoll*computed%rdrho + 1.d0)
  rhomax = (dfloat(irhomax) - 0.5d0)*input%drho
  ! Stay inside the sphere to avoid boundary issues
  irhomax = irhomax - 1
  rhomax = rhomax - input%drho
  write (*, *) 'rhomax,irhomax = ', rhomax, irhomax
  ict = 0
  zFo = 0.d0
  cho = 0.d0
  cliffFo = 0.d0
  rho = -0.5d0*input%dz

  ! Loop over radial slices from center outward (outer hemisphere in z)
  do irho = 1, irhomax
    rho = rho + input%drho
    rhosq = rho*rho
    ! Only process rho values that intersect the colloid
    if (rhosq .le. computed%Rcoll2) then
      z = input%zc1 + 0.5d0*input%dz
      iz = izc1
      ! Find first grid point outside colloid at this rho (moving down in z)
      do
        z = z - input%dz
        iz = iz - 1
        zsq = (z - input%zc1)**2
        if ((rhosq + zsq) .gt. computed%Rcoll2) exit
      end do
      Rc = dsqrt(rhosq + zsq)
      deltazc = dsqrt(computed%Rcoll2 - rhosq)

      ! Quadratic interpolation in z-direction to get surface density
      if (dabs(fdmon(irho, iz)) .gt. 0.00000001d0) then
        y3 = fdmon(irho, iz)
        y2 = fdmon(irho, iz - 1)
        y1 = fdmon(irho, iz - 2)
        x3 = dabs(z - input%zc1)
        x2 = x3 + input%dz
        x1 = x3 + 2.d0*input%dz
      else
        y3 = fdmon(irho, iz - 1)
        y2 = fdmon(irho, iz - 2)
        y1 = fdmon(irho, iz - 3)
        x3 = dabs(z - input%zc1) + input%dz
        x2 = x3 + input%dz
        x1 = x3 + 2.d0*input%dz
        write (*, *) 'TJOHO1!'
      end if

      x = deltazc
      fdc = y1*(x - x2)*(x - x3)/((x1 - x2)*(x1 - x3)) + &
            y2*(x - x1)*(x - x3)/((x2 - x1)*(x2 - x3)) + &
            y3*(x - x1)*(x - x2)/((x3 - x1)*(x3 - x2))
      ctheta = -deltazc/input%Rcoll
      ict = ict + 1
      ctvec(ict) = ctheta
      cdens(ict) = fdc
    end if
  end do

  ! Loop over radial slices in reverse (inner hemisphere in z)
  rho = rho + input%drho
  irho = irhomax + 1
  zFi = 0.d0
  chi = 0.d0
  cliffFi = 0.d0
  do krho = 1, irhomax
    irho = irho - 1
    rho = rho - input%drho
    rhosq = rho*rho
    ! Only process rho values that intersect the colloid
    if (rhosq .le. computed%Rcoll2) then
      z = input%zc1 - 0.5d0*input%dz
      iz = izc1 - 1
      ! Find first grid point outside colloid at this rho (moving up in z)
      do
        z = z + input%dz
        iz = iz + 1
        zsq = (z - input%zc1)**2
        if ((rhosq + zsq) .gt. computed%Rcoll2) exit
      end do
      Rc = dsqrt(rhosq + zsq)
      deltazc = dsqrt(computed%Rcoll2 - rhosq)

      if (dabs(fdmon(irho, iz)) .gt. 0.00000001d0) then
        y3 = fdmon(irho, iz)
        y2 = fdmon(irho, iz + 1)
        y1 = fdmon(irho, iz + 2)
        x3 = dabs(z - input%zc1)
        x2 = x3 + input%dz
        x1 = x3 + 2.d0*input%dz
      else
        y3 = fdmon(irho, iz + 1)
        y2 = fdmon(irho, iz + 2)
        y1 = fdmon(irho, iz + 3)
        x3 = dabs(z - input%zc1) + input%dz
        x2 = x3 + input%dz
        x1 = x3 + 2.d0*input%dz
        write (*, *) 'TJOHO2!'
      end if

      x = deltazc
      fdc = y1*(x - x2)*(x - x3)/((x1 - x2)*(x1 - x3)) + &
            y2*(x - x1)*(x - x3)/((x2 - x1)*(x2 - x3)) + &
            y3*(x - x1)*(x - x2)/((x3 - x1)*(x3 - x2))
      ctheta = deltazc/input%Rcoll
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
    rho = input%Rcoll*dsqrt(1.d0 - ct*ct)
    rhon = input%Rcoll*dsqrt(1.d0 - ctn*ctn)
    if (kct .eq. 0) rho = input%Rcoll
    if (kct .eq. ict) rhon = input%Rcoll
    if (dabs(rhon - rho) .lt. 0.00000001d0) then
      fk = 0.d0
    else
      fk = (fdcn - fdc)/(rhon - rho)
    end if
    ccc = (fdc - fk*rho)*(rhon - rho) + 0.5d0*fk*(rhon*rhon - rho*rho) + ccc
  end do

  ctF = 2.d0*PI*computed%Rcoll2*ctF
  ch2 = 2.d0*PI*computed%Rcoll2*ch2
  cckoll = 2.d0*PI*cckoll*input%drho
  ckoll = 2.d0*PI*ckoll*input%Rcoll*input%drho
  ckk = 2.d0*PI*ckk*input%Rcoll*input%drho
  ccc = 2.d0*PI*ccc
  write (*, *)
  write (*, *) 'ctF = ', ctF
  write (*, *)
  write (*, *) 'ch2 = ', ch2

  rewind ifc
  z = -0.5d0*input%dz
  do iz = grid%istp1, grid%imitt
    z = z + input%dz
    rho = -0.5d0*input%drho
    do kz = 1, grid%mxrho
      rho = rho + input%drho
      write (ifc, '(2f12.5,2f21.12)') &
        z, rho, fdmon(kz, iz), fem(kz, iz)
    end do
  end do

  ! Close files to ensure buffers are flushed
  close (ifc)
  close (ins)
  close (iep)

  ! Deallocate arrays before exit
  deallocate (c, cA, cB)
  call deallocate_arrays()

  STOP
END
