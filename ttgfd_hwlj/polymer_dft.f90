! MIT License
!
! Copyright (c) 2025 Jan Forsman and Mikael Lund
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
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

  ! Additional arrays for main program - dynamically allocated
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

  call read_input_parameters(input, bds)
  call initialize_grid_params(input, grid)

  ! Compute local derived parameters needed for bulk thermodynamic calculations
  ! These will be recomputed and stored in structs by initialize_computed_params later
  computed%dhs2 = input%dhs*input%dhs
  computed%dhs3 = computed%dhs2*input%dhs
  computed%rdhs3 = 1.d0/computed%dhs3
  computed%rnmon = dble(input%nmon)
  computed%rrnmon = 1.d0/computed%rnmon
  computed%Yfact = (computed%rnmon - 2.d0)*Y

  ! Initialize cosine lookup table for input%dphi
  allocate (cos_phi(grid%nphi))
  do iphi = 1, grid%nphi
    phi = (dble(iphi) - 0.5d0)*input%dphi
    cos_phi(iphi) = dcos(phi)
  end do

  call calculate_bulk_properties(input, computed, bulk)
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
  ! Dimensions calculated internally: nrho = mxrho + kbl, nz = imitt + ibl, nz_hvec = nfack - 1
  call allocate_arrays(fields, grid)

  ! Allocate main program arrays
  ! Note: cB needs grid%nfack dimension because it's accessed as grid%islut + 1 - iz
  ! c array sized based on actual number of monomers from input file
  allocate (c(0:grid%mxrho + grid%kbl, 0:grid%imitt + grid%ibl, input%nmon))
  allocate (cA(0:grid%mxrho + grid%kbl, 0:grid%imitt + grid%ibl))
  allocate (cB(0:grid%mxrho + grid%kbl, 0:grid%nfack))

  ! Calculate normalization constant for contact density
  call calculate_contact_density_normalization(input, grid, computed, cos_phi, computed%cdnorm)
  write (*, *) 'computed%cdnorm = ', computed%cdnorm

  call initialize_boundary_excess_free_energy(input, grid, computed, fields)

  open (newunit=ifc, file='fcdfil', form='formatted', status='unknown')
  rewind ifc
  call initialize_density_fields(input, grid, computed, fields, ifc)
  write (*, *) 'fields%fdmon(1,1) = ', fields%fdmon(1, 1)
  write (*, *) 'fields%fdmon(1,11) = ', fields%fdmon(1, 11)

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

    call calculate_contact_density(input, grid, computed, fields, cos_phi)
    call calculate_excess_free_energy(grid, computed, fields)
    call calculate_boltzmann_factors(input, grid, computed, fields, cos_phi)
    call calculate_external_potential(input, grid, computed, fields)
    call apply_boundary_conditions(grid, computed, fields)
    call propagate_polymer_chain(input, grid, computed, fields, cos_phi, c, cA, cB)

    if (ddmax .lt. CONV_TOL) exit ! Exit if converged

    call calculate_adaptive_mixing(niter, ddmax, input, use_adaptive_mixing, &
                                   oscillation_count, ddmax_prev, dmm_adaptive, dms_adaptive)

    call update_densities_and_check_convergence(input, grid, computed, fields, c, &
                                                dmm_adaptive, ddmax)

  end do

  ! Check if maximum iterations exceeded
  if (niter .gt. input%ioimaxm) then
    stop
  end if

  ! Output converged results
  call output_density_profiles(input, grid, computed, fields, c)
  call calculate_grand_potential(input, grid, computed, fields, bulk, aW, bW)
  call calculate_colloid_forces(input, computed, fields, rcliffF, ctF, ch2)

  ! Write fcdfil
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
  close (ifc)

  ! Deallocate arrays before exit
  deallocate (cos_phi, c, cA, cB)
  call deallocate_arrays(fields)

  STOP
END
