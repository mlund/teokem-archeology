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

  ! Bulk thermodynamic properties
  type(bulk_properties_t) :: bulk

  ! Integer variables
  integer(int32) :: i, j, k, iz, jz, kz, irho, iphi, imon, kmon
  integer(int32) :: irho0min, jstart
  integer(int32) :: ifc, ins, iep, niter, npphi
  integer(int32) :: iout_zdens, iout_zprop, iout_zavg, iout_rdens, iout_rprop

  ! Real variables
  real(real64) :: alj, aw
  real(real64) :: bds, bebbe
  real(real64) :: bw
  real(real64) :: ch2, ctf
  real(real64) :: ddiff, ddmax, delz2, diffz2
  real(real64) :: dmm_adaptive, dms_adaptive, dumsum
  real(real64) :: efact

  ! Logical variables
  logical :: use_adaptive_mixing
  integer(int32) :: oscillation_count
  real(real64) :: ddmax_prev
  real(real64) :: fact, ffact, fphi
  real(real64) :: phi, phisum
  real(real64) :: rcliffF, rcyl2
  real(real64) :: rho, rho0, rho02, rho2, rhoz2, rlj
  real(real64) :: rsq, rsq1, rsq2, rt2, strho0, valid
  real(real64) :: sume
  real(real64) :: tdmm, tdms, tfdm, tfem
  real(real64) :: z, zfact, zp, zpc2sq, zpcsq, zpst

  ! ========================================================================
  ! ========================================================================
  ! File unit numbers (automatically assigned by runtime)

  ! Open input and output files
  ! fcdfil: unknown allows both read and write (for restart capability)
  open (newunit=ifc, file='fcdfil', form='formatted', status='unknown')
  open (newunit=ins, file='input.tsph', form='formatted', status='old')
  open (newunit=iep, file='epfil', form='formatted', status='old')

  ! Open output files for density profiles (use traditional fort.* names for compatibility)
  open (newunit=iout_zdens, file='fort.85', form='formatted', status='replace')
  open (newunit=iout_zprop, file='fort.89', form='formatted', status='replace')
  open (newunit=iout_zavg, file='fort.78', form='formatted', status='replace')
  open (newunit=iout_rdens, file='fort.83', form='formatted', status='replace')
  open (newunit=iout_rprop, file='fort.87', form='formatted', status='replace')

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
  computed%dhs2 = input%dhs*input%dhs
  computed%dhs3 = computed%dhs2*input%dhs
  computed%rdhs3 = 1.d0/computed%dhs3
  computed%rnmon = dble(input%nmon)
  computed%rrnmon = 1.d0/computed%rnmon

  ! Read Lennard-Jones energy parameter and compute LJ coefficients
  read (iep, *) input%epslj
  alj = 4.d0*input%epslj*computed%dhs3*computed%dhs3    ! Attractive (r^-6) coefficient
  rlj = 4.d0*input%epslj*computed%dhs3**4      ! Repulsive (r^-12) coefficient

  ! Compute additional local parameters needed for bulk calculations
  Rcyl2 = input%Rcyl*input%Rcyl
  tdmm = 1.d0 - input%dmm
  tdms = 1.d0 - input%dms

  ! Initialize cosine lookup table for input%dphi (using grid%nphi from initialize_grid_params)
  do iphi = 1, grid%nphi
    phi = (dble(iphi) - 0.5d0)*input%dphi
    cos_phi(iphi) = dcos(phi)
  end do

  ! Calculate bulk thermodynamic properties
  computed%Yfact = (computed%rnmon - 2.d0)*Y
  call calculate_bulk_properties(input, computed, bulk)

  ! Initialize computed parameters from input, grid, and bulk calculations
  ! This recomputes and stores all derived parameters in structured form
  call initialize_computed_params(input, grid, computed, bulk%chempp, bulk%emtrams, bulk%cmtrams)

  ! Print simulation parameters
  write (*, *) 'GFD POLYMER SOLUTION MODEL!'
  write (*, *) 'input%bdm,bdpol =', input%bdm, bulk%bdpol
  write (*, *) 'monomer density  = ', input%bdm
  write (*, *) 'bdt = ', bulk%bdt
  write (*, *) 'collsep,input%dz = ', input%collsep, input%dz
  write (*, *) 'input%Rcoll = ', input%Rcoll
  write (*, *) 'bond length (input%bl): ', input%bl
  write (*, *) 'monomer hs diameter (input%bl): ', input%dhs
  write (*, *) 'no. of monomers/polymer = ', input%nmon
  write (*, *) 'max no. of iterations = ', input%ioimaxm
  write (*, *) 'polymer chemical pot. (betamu) = ', bulk%chempp
  write (*, *) 'solvent chemical pot. (betamu) = ', computed%chemps
  write (*, *) 'total bulk pressure = ', bulk%Pb
  write (*, *) 'input%zc1,computed%zc2 = ', input%zc1, computed%zc2
  write (*, *) 'grid%nfack,grid%imitt = ', grid%nfack, grid%imitt
  write (*, *) 'grid%istp1,grid%islut = ', grid%istp1, grid%islut
  write (*, *) 'grid%istp1s,grid%isluts = ', grid%istp1s, grid%isluts
  write (*, *) 'grid%ism,grid%ibl = ', grid%ism, grid%ibl
  write (*, *) 'grid%ksm,grid%kbl = ', grid%ksm, grid%kbl
  write (*, *) 'dmm,dms (density mixing param. mon.,solv.) = ', input%dmm, input%dms
  write (*, *) 'Rcyl  = ', input%Rcyl
  write (*, *) 'bFex = ', bulk%bFex
  write (*, *) 'computed%bebelam,computed%behbclam = ', computed%bebelam, computed%behbclam

  ! Allocate module arrays based on calculated grid dimensions
  ! Include extra space for boundary cells (grid%kbl, grid%ibl)
  ! hvec needs grid%nfack-1 for z-dimension (used in LJ potential table)
  call allocate_arrays(fields, grid%mxrho + grid%kbl, grid%imitt + grid%ibl, grid%nfack - 1)

  ! Allocate main program arrays
  ! Note: cB needs grid%nfack dimension because it's accessed as grid%islut + 1 - iz
  allocate (c(0:grid%mxrho + grid%kbl, 0:grid%imitt + grid%ibl, MAXMON))
  allocate (cA(0:grid%mxrho + grid%kbl, 0:grid%imitt + grid%ibl))
  allocate (cB(0:grid%mxrho + grid%kbl, 0:grid%nfack))

  ! Calculate normalization constant for contact density
  call CDFACT(input, grid, computed, cos_phi, computed%cdnorm)
  write (*, *) 'computed%cdnorm = ', computed%cdnorm

  ! Initialize excess free energy arrays at system boundaries
  call initialize_boundary_excess_free_energy(input, grid, computed, fields)

  ! Initialize density fields: either from scratch or read from file
  call initialize_density_fields(input, grid, computed, fields, ifc)
  write (*, *) 'fields%fdmon(1,1) = ', fields%fdmon(1, 1)
  write (*, *) 'fields%fdmon(1,11) = ', fields%fdmon(1, 11)

  ! Precompute Lennard-Jones interaction potential on grid (hvec array)
  ! This tabulates U_LJ for all distance combinations to speed up later calculations
  write (*, *) 'dpphi = ', input%dpphi
  call calculate_lj_potential_table(input, grid, computed, alj, rlj, fields, cos_pphi, npphi)
  write (*, *) 'dpphi,npphi = ', input%dpphi*PI, npphi
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
    call CDCALC(input, grid, computed, fields, cos_phi)     ! Calculate contact density
    call AVEC(grid, computed, fields)       ! Calculate excess free energy
    call EBLMNEW(input, grid, computed, fields, cos_phi)    ! Calculate end-segment Boltzmann factors
    call EBDU(input, grid, computed, fields)       ! Calculate external potential contribution

    ! Apply boundary conditions at outer radial edge (rho > Rcyl)
    ! Set to bulk values since density should approach bulk far from colloids
    do iz = grid%istp1, grid%imitt
    do kz = grid%mxrho + 1, grid%mxrho + grid%kbl
      fields%ebelam(kz, iz) = computed%bebelam
      fields%ehbclam(kz, iz) = computed%behbclam
      fields%edu(kz, iz) = 1.d0
    end do
    end do

    ! Apply symmetry boundary conditions at z = grid%imitt (midplane between colloids)
    ! The system is symmetric about the midplane, so mirror the field values
    ! This loop copies from iz = grid%imitt down to iz = 1 (reverse order)
    jz = grid%imitt + 1
    do iz = grid%imitt + 1, grid%imitt + grid%ibl
      jz = jz - 1
      do kz = 1, grid%mxrho + grid%kbl
        fields%ebelam(kz, iz) = fields%ebelam(kz, jz)
        fields%ehbclam(kz, iz) = fields%ehbclam(kz, jz)
        fields%edu(kz, iz) = fields%edu(kz, jz)
      end do
    end do

    ! Second symmetry application (appears redundant but ensures consistency)
    jz = grid%imitt + 1
    do iz = grid%imitt + 1, grid%imitt + grid%ibl
      jz = jz - 1
      do kz = 1, grid%mxrho + grid%kbl
        fields%ebelam(kz, iz) = fields%ebelam(kz, jz)
        fields%ehbclam(kz, iz) = fields%ehbclam(kz, jz)
        fields%edu(kz, iz) = fields%edu(kz, jz)
      end do
    end do

    ! Calculate cA: the propagator for polymer end segments
    ! cA(r) = exp(-beta*mu_end)*exp(-U_LJ) where:
    !   ebelam = exp(-beta*mu_end) from hard-sphere and chain connectivity
    !   edu = exp(-U_LJ) from Lennard-Jones interactions
    do iz = grid%istp1, grid%imitt + grid%ibl
    do kz = 1, grid%mxrho + grid%kbl
      cA(kz, iz) = fields%ebelam(kz, iz)*fields%edu(kz, iz)
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
          efact = dsqrt(fields%edu(kz, iz))
          ffact = sume*computed%dzrfp*fields%ehbclam(kz, iz)/input%bl*efact
          c(kz, iz, imon) = ffact
          if (iz .gt. grid%imitt - grid%ibl - 1) cB(kz, grid%islut + 1 - iz) = ffact*fields%ehbclam(kz, iz)*efact
          cB(kz, iz) = ffact*fields%ehbclam(kz, iz)*efact
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
          fields%fem(j, i) = 0.d0
          fields%fdmon(j, i) = 0.d0
        else
          ! Calculate total monomer density from chain propagators
          ! Sum over all internal segments (convolution of forward and backward propagators)
          dumsum = 0.d0
          do k = 2, input%nmon - 1
            dumsum = c(j, i, k)*c(j, i, input%nmon + 1 - k) + dumsum
          end do
          tfem = 2.d0*c(j, i, 1)*fields%ebelam(j, i)*dsqrt(fields%edu(j, i))/fields%ehbclam(j, i)
          tfdm = dumsum + tfem
          if (dabs(tfdm) .gt. 1.0d-14) then
            ddiff = abs(tfdm - fields%fdmon(j, i))/tfdm
            if (ddiff .gt. ddmax) ddmax = ddiff
          end if
          fields%fem(j, i) = fields%fem(j, i)*dmm_adaptive + tdmm*tfem
          fields%fdmon(j, i) = fields%fdmon(j, i)*dmm_adaptive + tdmm*tfdm
        end if
      end do
    end do

    ! Apply symmetry to updated densities at midplane
    ! Mirror fdmon and fem across z = grid%imitt to maintain symmetry
    jz = grid%imitt + 1
    do iz = grid%imitt + 1, grid%imitt + grid%ibl
      jz = jz - 1
      do kz = 1, grid%mxrho + grid%kbl
        fields%fdmon(kz, iz) = fields%fdmon(kz, jz)
        fields%fem(kz, iz) = fields%fem(kz, jz)
      end do
    end do

  end do  ! End of main iteration loop

  ! Check if maximum iterations exceeded
  if (niter .gt. input%ioimaxm) then
    stop
  end if

  ! ===== Output converged results =====
  ! Write density profiles and calculate thermodynamic properties
  call output_density_profiles(input, grid, computed, fields, c, &
                                iout_zdens, iout_zprop, iout_zavg, iout_rdens, iout_rprop)

  ! ===== Calculate grand potential (thermodynamic potential) =====
  call calculate_grand_potential(input, grid, computed, fields, bulk, aW, bW)

  ! ===== Calculate forces on colloid from contact density =====
  call calculate_colloid_forces(input, computed, fields, rcliffF, ctF, ch2)

  rewind ifc
  z = -0.5d0*input%dz
  do iz = grid%istp1, grid%imitt
    z = z + input%dz
    rho = -0.5d0*input%drho
    do kz = 1, grid%mxrho
      rho = rho + input%drho
      write (ifc, '(2f12.5,2f21.12)') &
        z, rho, fields%fdmon(kz, iz), fields%fem(kz, iz)
    end do
  end do

  ! Close files to ensure buffers are flushed
  close (ifc)
  close (ins)
  close (iep)
  close (iout_zdens)
  close (iout_zprop)
  close (iout_zavg)
  close (iout_rdens)
  close (iout_rprop)

  ! Deallocate arrays before exit
  deallocate (c, cA, cB)
  call deallocate_arrays(fields)

  STOP
END
