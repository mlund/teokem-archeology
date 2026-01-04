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
  integer(int32) :: iz, kz, iphi
  integer(int32) :: ifc, niter

  ! Real variables
  real(real64) :: aw, bw
  real(real64) :: bds
  real(real64) :: ch2, ctf
  real(real64) :: ddmax
  real(real64) :: dmm_adaptive, dms_adaptive

  ! Adaptive mixing state variables
  logical :: use_adaptive_mixing
  integer(int32) :: oscillation_count
  real(real64) :: ddmax_prev
  real(real64) :: phi, rcliffF, rho, z

  ! ========================================================================
  ! ========================================================================
  ! File unit numbers (automatically assigned by runtime)

  ! Read all simulation parameters from input files
  call read_input_parameters(input, bds)

  ! Initialize grid parameters from input
  call initialize_grid_params(input, grid)

  ! Compute local derived parameters needed for bulk thermodynamic calculations
  ! These will be recomputed and stored in structs by initialize_computed_params later
  computed%dhs2 = input%dhs*input%dhs
  computed%dhs3 = computed%dhs2*input%dhs
  computed%rdhs3 = 1.d0/computed%dhs3
  computed%rnmon = dble(input%nmon)
  computed%rrnmon = 1.d0/computed%rnmon

  ! Open files for I/O
  ! fcdfil: unknown allows both read and write (for restart capability)
  open (newunit=ifc, file='fcdfil', form='formatted', status='unknown')
  rewind ifc

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
  call initialize_computed_params(input, grid, computed, bulk)

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
  ! LJ coefficients (alj, rlj) and angular grid (cos_pphi, npphi) are computed internally
  write (*, *) 'dpphi = ', input%dpphi
  call calculate_lj_potential_table(input, grid, computed, fields)
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

    ! Apply boundary conditions (radial + z-symmetry)
    call apply_boundary_conditions(grid, computed, fields)

    ! Propagate polymer chains segment by segment
    call propagate_polymer_chain(input, grid, computed, fields, cos_phi, c, cA, cB)

    if (ddmax .lt. CONV_TOL) exit  ! Converged

    ! Calculate adaptive mixing parameters
    call calculate_adaptive_mixing(niter, ddmax, input, use_adaptive_mixing, &
                                    oscillation_count, ddmax_prev, dmm_adaptive, dms_adaptive)

    ! Update densities and check convergence
    call update_densities_and_check_convergence(input, grid, computed, fields, c, &
                                                 dmm_adaptive, ddmax)

  end do  ! End of main iteration loop

  ! Check if maximum iterations exceeded
  if (niter .gt. input%ioimaxm) then
    stop
  end if

  ! ===== Output converged results =====
  ! Write density profiles and calculate thermodynamic properties
  call output_density_profiles(input, grid, computed, fields, c)

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

  ! Deallocate arrays before exit
  deallocate (c, cA, cB)
  call deallocate_arrays(fields)

  STOP
END
