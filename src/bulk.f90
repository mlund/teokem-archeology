! =============================================================================
! =============================================================================
!
! MONTE CARLO SIMULATION OF MULTICOMPONENT ISOTROPIC IONIC SOLUTION
!
! =============================================================================
! =============================================================================


program bulk

   implicit none
   real :: ran2
   include 'bulk_f90.inc'

   ! Explicit interfaces for subroutines
   interface
      subroutine calculate_statistics(data_array, standard_deviation, mean_value, array_size)
         integer :: array_size
         double precision :: data_array(25), standard_deviation, mean_value
      end subroutine calculate_statistics

      subroutine calculate_statistics_per_species(a, bbbwww, xb, nnn, num)
         integer, parameter :: max_species = 10
         integer :: nnn, num
         double precision :: a(25, max_species), xb(max_species), bbbwww(max_species)
      end subroutine calculate_statistics_per_species

      subroutine initialize_random_configuration()
      end subroutine initialize_random_configuration

      subroutine evaluate_trial_move()
      end subroutine evaluate_trial_move

      subroutine recalculate_total_energy()
      end subroutine recalculate_total_energy

      subroutine calculate_collision_pressure()
      end subroutine calculate_collision_pressure

      subroutine calculate_collision_pressure1()
      end subroutine calculate_collision_pressure1

      subroutine calculate_collision_pressure2()
      end subroutine calculate_collision_pressure2

      subroutine calculate_collision_pressure3()
      end subroutine calculate_collision_pressure3

      subroutine calculate_widom_insertion()
      end subroutine calculate_widom_insertion

      subroutine calculate_widom_insertion1()
      end subroutine calculate_widom_insertion1

      subroutine calculate_widom_insertion2()
      end subroutine calculate_widom_insertion2

      subroutine calculate_widom_insertion3()
      end subroutine calculate_widom_insertion3
   end interface

   ! Note: Variable names have been updated to be more descriptive
   ! See bulk_f90.inc for the complete variable declarations

   ! Local variables for Monte Carlo simulation
   integer :: step_inner, step_middle                        ! Loop counters for inner and middle MC loops
   integer :: total_hardcore_rejections                      ! Total hard core overlap rejections
   integer :: total_energy_rejections                        ! Total energy-based rejections
   integer :: total_accepted_moves                           ! Total accepted Monte Carlo moves
   integer :: total_mc_steps                                 ! Total Monte Carlo steps attempted
   integer :: species_index, particle_index                  ! Loop indices for species and particles
   integer :: loop_i, loop_j, loop_k, loop_l                 ! Generic loop counters
   integer :: particle_count_temp                            ! Temporary particle count
   integer :: total_configurations                           ! Total number of configurations
   integer :: total_widom_tests                              ! Total Widom test particle insertions
   integer :: moves_per_species_total(max_species)           ! Total moves attempted per species
   integer :: moves_per_species_accepted(max_species)        ! Accepted moves per species
   integer :: moves_per_species_energy_rejected(max_species) ! Energy rejected moves per species
   integer :: moves_per_species_hardcore_rejected(max_species) ! Hard core rejected moves per species

   double precision :: energy_accumulator_current(2)         ! Energy accumulator for current macro step
   double precision :: energy_accumulator_total(2)           ! Total energy accumulator across all macro steps
   double precision :: energy_per_macrostep(25)              ! Energy stored for each macro step
   double precision :: time_array(10)                        ! Time measurements array (currently unused)
   double precision :: energy_difference                     ! Energy difference for MC acceptance
   double precision :: boltzmann_factor                      ! Boltzmann factor for MC acceptance
   double precision :: energy_variance                       ! Variance of energy measurements
   double precision :: temp_kelvin                           ! Temperature in Kelvin
   double precision :: pressure_excess_energy                ! Excess pressure from energy
   double precision :: pressure_excess_variance              ! Variance of excess pressure
   double precision :: pressure_ideal                        ! Ideal gas pressure contribution
   double precision :: pressure_total                        ! Total system pressure
   double precision :: pressure_total_variance               ! Variance of total pressure
   double precision :: energy_consistency_check              ! Relative energy consistency check
   double precision :: macro_step_inverse                    ! Inverse of number of macro steps (1/N)
   double precision :: initial_energy_scaled                 ! Initial configuration energy (scaled to kJ/mol)

   data energy_accumulator_total  / 2*0.0 /
   data energy_per_macrostep   / 25*0.0 /
   data moves_per_species_total, moves_per_species_accepted, &
        moves_per_species_energy_rejected, moves_per_species_hardcore_rejected / 40*0 /
   data total_hardcore_rejections, total_energy_rejections, &
        total_accepted_moves, total_mc_steps / 4*0 /
   data time_array  / 10*0.0 /

43 format( 'MONTE CARLO SIMULATION OF MULTICOMPONENT ISOTROPIC    ', &
            'IONIC SYSTEM  ', /, /, 'Written by Peter Bolhuis, February', &
            '1992', / )

44 format( '************************************************************', &
            '************************', / )

   unit_macro = 12
   unit_conf = 13
   unit_input = 15
   unit_output = 6

   write(unit_output, 44)
   write(unit_output, 43)
   write(unit_output, 44)


   ! ==========================================================================
   ! Physical constants are now defined as parameters in bulk_f90.inc
   ! ==========================================================================


   ! ==========================================================================
   ! Sample units are angstrom and elementary charge
   ! Kilojoules are used for a mole
   !
   ! MEANING OF THE INPUT VARIABLES
   !
   ! initial_config_type        = 0   from file
   ! initial_config_type        = 1   from slump
   ! num_species      = number of ionic components
   ! hion        = array contains the input info on the ions
   !               (1=nr, 2=number, 3=radius, 4=charge, 5=displacement)
   ! ny(1,2,3)   = number of MC loops
   ! temperature      = temperature
   ! dielectric_constant        = relative dielectric constant
   ! box_size        = box_sizesize
   ! widom_interval      = a widom/collision calc is performed every widom_intervalcycle
   ! num_widom_insertions      = number of widom particles inserted each calc.
   ! measurement_location       = where is the chem. pot. measured
   !               0 = averaged over the cell
   !
   ! ==========================================================================
   ! ==========================================================================

   open(unit_conf, file='bulk.conf', form='formatted')
   open(unit_input, file='bulk.inp',  form='formatted')
   open(unit_macro, file='bulk.macro', form='formatted')

   read(unit_input, *)
   read(unit_input, *) num_species
   read(unit_input, *)
   read(unit_input, *)
   read(unit_input, *) (species_properties(loop_l, 1), species_properties(loop_l, 2), species_properties(loop_l, 3), species_properties(loop_l, 4), species_properties(loop_l, 5), loop_l = 1, num_species)
   read(unit_input, *)
   read(unit_input, *)
   read(unit_input, *) mc_steps_inner, mc_steps_middle, mc_steps_outer
   read(unit_input, *)
   read(unit_input, *)
   read(unit_input, *) initial_config_type, random_seed
   random_seed = -random_seed
   read(unit_input, *)
   read(unit_input, *)
   read(unit_input, *) box_size
   read(unit_input, *)
   read(unit_input, *)
   read(unit_input, *) temperature, dielectric_constant
   read(unit_input, *)
   read(unit_input, *)
   read(unit_input, *) num_widom_insertions, widom_interval, measurement_location

   if (measurement_location> 2)  measurement_location= 2
   if (random_seed == 0) random_seed = 7
   if (widom_interval== 0) widom_interval= 2 * mc_steps_inner


   ! ==========================================================================
   ! Boxsize variables
   ! ==========================================================================

   box_half = box_size/ 2.0
   box_half_inverse = 1.0 / box_half


   ! ==========================================================================
   ! Set up hard core and charge vectors
   !
   ! hc2v contains in the first column the particle type given
   ! by species_properties(i,1), whereas in the following columns the squared
   ! hard core distances between the particles are given.
   ! chv gives for every particle its charge
   ! num_particlesis the total number of particles
   ! ==========================================================================

   num_particles= 0

   do species_index = 1, num_species
      particle_count_temp = int(species_properties(species_index, 2))

      do particle_index = 1, particle_count_temp
         num_particles= num_particles+ 1
         particle_species(num_particles) = species_index

         do loop_k = 1, num_species
            hard_core_distance_sq(num_particles, loop_k) = (species_properties(species_index, 3) + species_properties(loop_k, 3)) ** 2
         end do

         particle_charge(num_particles) = species_properties(species_index, 4)
         displacement_max(num_particles)  = species_properties(species_index, 5)
      end do
   end do

   do species_index = 1, num_species
      species_concentration(species_index) = (species_properties(species_index, 2)) / (box_size** 3)
   end do

   do particle_index = 1, num_particles
      trial_energy(particle_index) = 0.0
   end do

   total_configurations = mc_steps_inner * mc_steps_middle * mc_steps_outer

   if (mc_steps_outer > 25) mc_steps_outer = 25

   temp_kelvin = temperature
   beta_inverse_temp= -1.0 / (gas_constant * temperature)


   ! ==========================================================================
   ! Read configuration from file
   ! ==========================================================================

   if (initial_config_type== 0) then
      read(unit_conf, *) (particle_x(loop_l), particle_y(loop_l), particle_z(loop_l), loop_l = 1, num_particles)
      write(unit_output, *) 'Configuration read from file'
   end if

   if (initial_config_type== 1) call initialize_random_configuration
   call calculate_collision_pressure
   call calculate_widom_insertion


   ! ==========================================================================
   ! Write input data
   ! ==========================================================================

   total_widom_tests = num_widom_insertions* (int(mc_steps_inner / widom_interval) * mc_steps_middle * mc_steps_outer)

   write(unit_output, 800) num_particles
   write(unit_output, 801) (loop_l, int(species_properties(loop_l, 2)), species_properties(loop_l, 3), species_properties(loop_l, 4), &
                    species_properties(loop_l, 5), species_concentration(loop_l) * 1.0d27 / avogadro_number, loop_l = 1, num_species)
   write(unit_output, 802) dielectric_constant, temperature, box_size, mc_steps_inner, mc_steps_middle, mc_steps_outer, total_configurations, &
                   dble(total_configurations / num_particles), total_widom_tests

800 format (/, 'VARIABLES FROM INPUT FILE', /, &
            /, ' number of particles        =', i10, /, &
            /, ' species   number   radius    charge   dp   average conc', /)
801 format (i6, i10, f10.1, f10.1, f10.2, f10.6)
802 format (/, ' dielectric constant         =', f10.1, /, &
            ' temperature                 =', f10.1, /, &
            ' box_sizesize (AA)               =', f10.1, /, &
            ' number of configurations    =', i6, ' *', i6, ' *', i6, ' =', i8, /, &
            ' number of conf. per part    =', f10.1, /, &
            ' widom test per species      =', i10, /)

   energy_conversion_factor= elementary_charge** 2 * 1.0d10 * avogadro_number/ (dielectric_constant * vacuum_permittivity* 4.0 * pi)


   ! ==========================================================================
   ! Energy evaluation initial configuration
   ! ==========================================================================

   call recalculate_total_energy

   initial_energy_scaled = total_coulomb_energy* energy_conversion_factor* 1.0d-3

   write(unit_output, '(a)') 'Initial configuration'
   write(unit_output, '(a, e12.5)') 'Coulomb energy (kJ/mol)  =', initial_energy_scaled


   ! ==========================================================================
   ! Start of Monte Carlo loop
   !
   ! Loop over mc_steps_inner*mc_steps_middle*mc_steps_outer particle moves.
   ! current_particle   = current particle
   ! ei    = new coulomb energy of particle il
   ! esa   = all coulomb interaction energies of the particle
   ! total_coulomb_energy = total energy (also ulaa(1))
   ! displacement_max(1..10)  displacement parameter for every ion
   !
   ! ==========================================================================

   do current_macro_step = 1, mc_steps_outer
      energy_accumulator_current(1) = 0.0

      do step_middle = 1, mc_steps_middle
         do step_inner = 1, mc_steps_inner
            current_particle   = 1 + mod(total_mc_steps, num_particles)
            current_species= particle_species(current_particle)
            total_mc_steps = total_mc_steps + 1
            moves_per_species_total(current_species) = moves_per_species_total(current_species) + 1

            trial_x= particle_x(current_particle) + displacement_max(current_particle) * (ran2(random_seed) - 0.5)
            trial_y= particle_y(current_particle) + displacement_max(current_particle) * (ran2(random_seed) - 0.5)
            trial_z= particle_z(current_particle) + displacement_max(current_particle) * (ran2(random_seed) - 0.5)

            if (trial_x>  box_half) trial_x= trial_x- box_size
            if (trial_x< -box_half) trial_x= trial_x+ box_size
            if (trial_y>  box_half) trial_y= trial_y- box_size
            if (trial_y< -box_half) trial_y= trial_y+ box_size
            if (trial_z>  box_half) trial_z= trial_z- box_size
            if (trial_z< -box_half) trial_z= trial_z+ box_size

            overlap_status= 99
            call evaluate_trial_move

            if (overlap_status/= 0) go to 63

            energy_difference = 0.0
            do loop_j = 1, num_particles
               energy_difference = energy_difference + trial_energy(loop_j) - energy_matrix(loop_j, current_particle)
            end do

            boltzmann_factor = beta_inverse_temp* (energy_difference * energy_conversion_factor)

            if (boltzmann_factor < -80.0) go to 64
            if (boltzmann_factor >   0.0) go to 62
            if (exp(boltzmann_factor) < ran2(random_seed)) go to 64

62          moves_per_species_accepted(current_species) = moves_per_species_accepted(current_species) + 1

            ! Trial config accepted
            particle_x(current_particle) = trial_x
            particle_y(current_particle) = trial_y
            particle_z(current_particle) = trial_z

            do loop_j = 1, num_particles
               energy_matrix(current_particle, loop_j) = trial_energy(loop_j)
               energy_matrix(loop_j, current_particle) = trial_energy(loop_j)
            end do

            total_coulomb_energy= total_coulomb_energy+ energy_difference
            go to 65

63          moves_per_species_hardcore_rejected(current_species) = moves_per_species_hardcore_rejected(current_species) + 1
            go to 65

64          moves_per_species_energy_rejected(current_species) = moves_per_species_energy_rejected(current_species) + 1

65          energy_accumulator_current(1) = energy_accumulator_current(1) + total_coulomb_energy

            if (mod(step_inner, widom_interval) == 0) then
               call calculate_collision_pressure1
               call calculate_widom_insertion1
            end if

         end do
      end do

      energy_consistency_check = total_coulomb_energy
      call recalculate_total_energy

      if (total_coulomb_energy/= 0) energy_consistency_check = (energy_consistency_check - total_coulomb_energy) / total_coulomb_energy

      write(unit_macro, *)
      write(unit_macro, 44)
      write(unit_macro, 73) current_macro_step, energy_consistency_check

      if (abs(energy_consistency_check) > 0.001) stop

      call recalculate_total_energy

      if (widom_interval<= mc_steps_inner) then
         call calculate_collision_pressure2
         call calculate_widom_insertion2
      end if

      energy_accumulator_current(1) = energy_accumulator_current(1) / dble(mc_steps_inner * mc_steps_middle)
      energy_accumulator_total(1) = energy_accumulator_total(1) + energy_accumulator_current(1)
      energy_accumulator_current(1) = energy_accumulator_current(1) * energy_conversion_factor* 1.0d-3
      energy_per_macrostep(current_macro_step) = energy_accumulator_current(1)

      ! Calculate total acceptance and rejection
      total_accepted_moves = 0
      total_hardcore_rejections = 0
      total_energy_rejections = 0

      do loop_i = 1, num_species
         total_accepted_moves = total_accepted_moves + moves_per_species_accepted(loop_i)
         total_hardcore_rejections = total_hardcore_rejections + moves_per_species_hardcore_rejected(loop_i)
         total_energy_rejections = total_energy_rejections + moves_per_species_energy_rejected(loop_i)
      end do

      write(unit_macro, 810) total_mc_steps, total_accepted_moves, total_energy_rejections, total_hardcore_rejections
      write(unit_macro, '(a, e12.5)') 'Coulomb energy (kj/mol)  =', energy_accumulator_current(1)
   end do

73  format(/ ' MACROSTEP  ', i10, /, /, 'Checkparameters : ', e10.3)
810 format(/, ' total number of configurations    =', i8, /, &
           ' accepted configurations           =', i8, /, &
           ' energy rejected configurations    =', i8, /, &
           ' hard core rejected configurations =', i8, /)
811 format(/, 'Divided per species', /, 'species      accepted   ', &
           ' energy rej.   hard core rej.', /, 10(i6, 3(6x, f8.4), /), /)


   ! ==========================================================================
   ! End of Monte Carlo loop
   !
   ! Calculation of the average energies and the pressures
   !
   ! ==========================================================================

   macro_step_inverse = 1.0 / float(mc_steps_outer)
   energy_accumulator_total(1) = macro_step_inverse * energy_accumulator_total(1)
   energy_accumulator_total(1) = energy_accumulator_total(1) * energy_conversion_factor* 1.0d-3

   call calculate_statistics(energy_per_macrostep, energy_variance, energy_accumulator_total(1), mc_steps_outer)

   write(unit_output, '(/)')
   write(unit_output, 44)
   write(unit_output, '(/, a, i3, a, /)') 'FINAL RESULTS AFTER ', mc_steps_outer, &
                                  ' MACROSTEPS'
   write(unit_output, 810) total_mc_steps, total_accepted_moves, total_energy_rejections, total_hardcore_rejections
   write(unit_output, 811) (loop_l, dble(moves_per_species_accepted(loop_l)) / moves_per_species_total(loop_l), &
                    dble(moves_per_species_energy_rejected(loop_l)) / moves_per_species_total(loop_l), &
                    dble(moves_per_species_hardcore_rejected(loop_l)) / moves_per_species_total(loop_l), loop_l = 1, num_species)
   write(unit_output, '(a, 2e12.5)') 'Coulomb energy (kj/mol)  =', energy_accumulator_total(1), energy_variance
   write(unit_output, *)

   if (widom_interval<= mc_steps_inner) then
      call calculate_collision_pressure3
      call calculate_widom_insertion3
   end if

   rewind unit_conf
   write(unit_conf, 771) (particle_x(loop_l), particle_y(loop_l), particle_z(loop_l), loop_l = 1, num_particles)
771 format(5e16.8)

   close(unit_conf)


   ! ==========================================================================
   !
   ! Calculation of pressures from the densities
   !
   ! ==========================================================================

   pressure_ideal = 0.0
   do loop_k = 1, num_species
      pressure_ideal = pressure_ideal + species_concentration(loop_k) * 1.0d27 / avogadro_number
   end do

   pressure_excess_energy  = energy_accumulator_total(1) * 1.0d30 / (3 * gas_constant* temperature* box_size** 3 * avogadro_number)
   pressure_excess_variance = energy_variance   * 1.0d30 / (3 * gas_constant* temperature* box_size** 3 * avogadro_number)
   pressure_total   = pressure_ideal + pressure_collision_average+ pressure_excess_energy
   pressure_total_variance  = sqrt(pressure_collision_variance* pressure_collision_variance+ pressure_excess_variance * pressure_excess_variance)

   write(unit_output, '(/, a, /)') 'TOTAL BULK PRESSURE '
   write(unit_output, 729) 'Ideal pressure       ', pressure_ideal, 0.0
   write(unit_output, 729) 'Energy con. <E/3V>   ', pressure_excess_energy, pressure_excess_variance
   write(unit_output, 729) 'Collision press      ', pressure_collision_average, pressure_collision_variance
   write(unit_output, 731)
   write(unit_output, 729) 'Bulk pressure        ' , pressure_total, pressure_total_variance

729 format(a, f12.4, f10.4)
731 format(70('_'))

   stop

end program bulk



! =============================================================================
! =============================================================================
!
! Subroutine: earth
!
! Calculates the average and deviation of nnn values stored in
! array a(1..25)
!
! =============================================================================
! =============================================================================

subroutine calculate_statistics(data_array, standard_deviation, mean_value, array_size)

   implicit none
   integer :: loop_index, array_size
   double precision :: data_array(25), standard_deviation, mean_value
   double precision :: sum_squared_deviations, normalization_factor

   normalization_factor = 1.0 / (array_size * (array_size - 1))
   sum_squared_deviations   = 0.0
   mean_value  = 0.0

   do loop_index = 1, array_size
      mean_value = mean_value + data_array(loop_index)
   end do

   mean_value = mean_value / array_size

   do loop_index = 1, array_size
      sum_squared_deviations = sum_squared_deviations + (data_array(loop_index) - mean_value) * (data_array(loop_index) - mean_value)
   end do

   sum_squared_deviations = sqrt(normalization_factor * sum_squared_deviations)
   standard_deviation = sum_squared_deviations

   return

end subroutine calculate_statistics


! =============================================================================
! =============================================================================
!
! Subroutine: earth2
!
! Calculates the average and deviation of nnn values stored in
! matrix a(1..25,1..num) and stores them in arrays xb and bbbwww
!
! =============================================================================
! =============================================================================

subroutine calculate_statistics_per_species(a, bbbwww, xb, nnn, num)

   implicit none
   integer, parameter :: max_species = 10
   integer :: i, k, nnn, num
   double precision :: a(25, max_species), xb(max_species), bbbwww(max_species)
   double precision :: b, yak

   yak = 1.0 / (nnn * (nnn - 1))

   do i = 1, num
      b = 0.0
      xb(i) = 0.0

      do k = 1, nnn
         xb(i) = xb(i) + a(k, i)
      end do

      xb(i) = xb(i) / nnn

      do k = 1, nnn
         b = b + (a(k, i) - xb(i)) * (a(k, i) - xb(i))
      end do

      b = sqrt(yak * b)
      bbbwww(i) = b
   end do

   return

end subroutine calculate_statistics_per_species


! =============================================================================
! =============================================================================
!
! Subroutine: slump
!
! Set up of the random initial configuration.
!
! =============================================================================
! =============================================================================

subroutine initialize_random_configuration

   implicit none
   real :: ran2
   include 'bulk_f90.inc'

   ! Local variables
   integer :: i, k12, nsl, nslmax
   double precision :: ddx, ddy, ddz, r2, x6tt, y6tt, z6tt

   ! Note: Variable names have been updated to be more descriptive
   ! See bulk_f90.inc for the complete variable declarations

   write(unit_output, *) 'Random configuration generated'

   nsl    = 0
   nslmax = 100000
   k12    = 0

1  nsl = nsl + 1

   if (nsl > nslmax) then
      write(unit_output, *)' too dense system'
      stop
   end if

   x6tt  = (ran2(random_seed) - 0.5) * box_size
   y6tt  = (ran2(random_seed) - 0.5) * box_size
   z6tt  = (ran2(random_seed) - 0.5) * box_size
   current_species= particle_species(k12 + 1)

   do i = 1, k12
      ddx = x6tt - particle_x(i)
      ddy = y6tt - particle_y(i)
      ddz = z6tt - particle_z(i)
      ddx = ddx - aint(ddx * box_half_inverse) * box_size
      ddy = ddy - aint(ddy * box_half_inverse) * box_size
      ddz = ddz - aint(ddz * box_half_inverse) * box_size
      r2  = ddx ** 2 + ddy ** 2 + ddz ** 2

      if (r2 < hard_core_distance_sq(i, current_species)) go to 1
   end do

   k12 = k12 + 1
   particle_x(k12) = x6tt
   particle_y(k12) = y6tt
   particle_z(k12) = z6tt

   if (k12 == num_particles) return

   go to 1

end subroutine initialize_random_configuration


! =============================================================================
! =============================================================================
!
! Subroutine: qin
!
! Essential part of the program. Checks for overlap and calculate
! the energy change after a MC step.
!
! =============================================================================
! =============================================================================

subroutine evaluate_trial_move

   implicit none
   include 'bulk_f90.inc'

   ! Local variables
   integer :: k
   double precision :: ddx, ddy, ddz, chil

   ! Note: Variable names have been updated to be more descriptive
   ! See bulk_f90.inc for the complete variable declarations

   do k = 1, num_particles
      ddx = abs(trial_x- particle_x(k))
      ddy = abs(trial_y- particle_y(k))
      ddz = abs(trial_z- particle_z(k))

      if (ddx > box_half) ddx = ddx - box_size
      if (ddy > box_half) ddy = ddy - box_size
      if (ddz > box_half) ddz = ddz - box_size

      distance_squared(k) = ddx * ddx + ddy * ddy + ddz * ddz

      ! Cancel out the hard core overlap with itself
      distance_squared(current_particle) = 1000000

      if (distance_squared(k) < hard_core_distance_sq(k, current_species)) return
   end do

   overlap_status= 0
   chil = particle_charge(current_particle)

   do k = 1, num_particles
      trial_energy(k) = chil * particle_charge(k) / sqrt(distance_squared(k))
   end do

   trial_energy(current_particle) = 0

   return

end subroutine evaluate_trial_move


! =============================================================================
! =============================================================================
!
! Subroutine: liv
!
! Gives the total coulombic energy and the force by recalculating
! every particleinteraction. Used to check QIN
!
! =============================================================================
! =============================================================================

subroutine recalculate_total_energy

   implicit none
   include 'bulk_f90.inc'

   ! Local variables
   integer :: i, k, ip1
   double precision :: ddx, ddy, ddz, uj1

   ! Note: Variable names have been updated to be more descriptive
   ! See bulk_f90.inc for the complete variable declarations

   total_coulomb_energy= 0.0

   do i = 1, num_particles- 1
      trial_x= particle_x(i)
      trial_y= particle_y(i)
      trial_z= particle_z(i)
      ip1 = i + 1

      do k = ip1, num_particles
         ddx = trial_x- particle_x(k)
         ddy = trial_y- particle_y(k)
         ddz = trial_z- particle_z(k)
         ddx = ddx - aint(ddx * box_half_inverse) * box_size
         ddy = ddy - aint(ddy * box_half_inverse) * box_size
         ddz = ddz - aint(ddz * box_half_inverse) * box_size
         distance_squared(k) = ddx * ddx + ddy * ddy + ddz * ddz
      end do

      do k = ip1, num_particles
         uj1 = particle_charge(i) * particle_charge(k) / sqrt(distance_squared(k))
         total_coulomb_energy= total_coulomb_energy+ uj1
         energy_matrix(i, k) = uj1
      end do
   end do

   do i = 1, num_particles- 1
      do k = i + 1, num_particles
         energy_matrix(k, i) = energy_matrix(i, k)
      end do
   end do

   do i = 1, num_particles
      energy_matrix(i, i) = 0.0
   end do

   return

end subroutine recalculate_total_energy


! =============================================================================
!
! Subroutine: collision
!
! This routine checkes if a particle of type i in one halfcel is
! able to collide with a patricle of type j in the other. The
! contribution to the pressure is calculation by intergration
! At a random position around the real particle of type i a ghost
! particle of type j is inserted and the potential with every
! other particle in the system is calculated. This gives the
! collision probability
!
! =============================================================================
! =============================================================================

subroutine calculate_collision_pressure

   implicit none
   real :: ran2
   include 'bulk_f90.inc'

   ! Local variables
   integer :: i, j, k, isp, ksp, nisp, nnn, num, nwtot
   double precision :: rel(max_species)
   double precision :: scoll(25, max_species, max_species), coll(max_species, max_species)
   double precision :: ddx, ddy, ddz, aa, dis, v, g2, wx6, wy6, wz6, urej, wtot2
   double precision :: collav, collv

   ! Note: Variable names have been updated to be more descriptive
   ! See bulk_f90.inc for the complete variable declarations

   do i = 1, max_species
      do k = 1, max_species
         coll(i, k) = 0.0

         do j = 1, mc_steps_outer
            scoll(j, i, k) = 0.0
         end do
      end do
   end do

   return


   ! -------------------------------------------------------------------------

   entry calculate_collision_pressure1

   do i = 1, num_widom_insertions
      num = 1

      do isp = 1, num_species
         nisp = int(species_properties(isp, 2))

         do ksp = 1, num_species
            nnn = int(ran2(random_seed) * nisp) + num
            dis = species_properties(isp, 3) + species_properties(ksp, 3)
            aa  = (2 * ran2(random_seed) - 1) * dis
            v   = 2 * pi * ran2(random_seed)
            wz6 = particle_z(nnn) + aa
            g2  = sqrt(dis * dis - aa * aa) + 0.00001
            wx6 = particle_x(nnn) + g2 * cos(v)
            wy6 = particle_y(nnn) + g2 * sin(v)
            wx6 = wx6 - aint(wx6 * box_half_inverse) * box_size
            wy6 = wy6 - aint(wy6 * box_half_inverse) * box_size
            wz6 = wz6 - aint(wz6 * box_half_inverse) * box_size
            urej = 0

            do k = 1, num_particles
               ddx = dabs(wx6 - particle_x(k))
               ddy = dabs(wy6 - particle_y(k))
               ddz = dabs(wz6 - particle_z(k))

               if (ddx > box_half) ddx = ddx - box_size
               if (ddy > box_half) ddy = ddy - box_size
               if (ddz > box_half) ddz = ddz - box_size

               distance_squared(k) = ddx * ddx + ddy * ddy + ddz * ddz

               if (distance_squared(k) < hard_core_distance_sq(k, ksp)) go to 121
            end do

            do k = 1, num_particles
               distance_inverse(k) = 1.0 / sqrt(distance_squared(k))
            end do

            wtot2 = 0

            do k = 1, num_particles
               wtot2 = wtot2 + distance_inverse(k) * particle_charge(k)
            end do

            coll(isp, ksp) = exp(beta_inverse_temp* wtot2 * species_properties(ksp, 4) * energy_conversion_factor) + coll(isp, ksp)

121         continue
         end do

         num = num + nisp
      end do
   end do

   return


   ! -------------------------------------------------------------------------

   entry calculate_collision_pressure2

   nwtot = num_widom_insertions* int(mc_steps_inner / widom_interval) * mc_steps_middle

   write (unit_macro, '(/, /, a)') 'COLLISION PRESSURE MATRIX'
   write (unit_macro, '(/, a, i6)') 'Total collision trials per species ', nwtot
   write (unit_macro, 2010) (i, i = 1, num_species)

   do k = 1, num_species
      do i = 1, num_species
         scoll(current_macro_step, i, k) = coll(i, k) / nwtot
         coll(i, k) = 0
      end do

      write(unit_macro, 2012) (scoll(current_macro_step, i, k), i = 1, num_species)
   end do

   return

2010 format('Species          ', 10(12x, i6), /)
2012 format('        pressure    ', 10e18.6, /)
2013 format(i4, '        ', 20e12.5)


   ! -------------------------------------------------------------------------

   entry calculate_collision_pressure3

   write(unit_output, '(/, a, /)') 'COLLISION MATRIX  AND RELATIVE ERROR'
   write (unit_output, 2010) (i, i = 1, num_species)

   do i = 1, num_species
      do k = 1, num_species
         collav = 0

         do j = 1, mc_steps_outer
            contact_correlation(j, i, k) = scoll(j, i, k)
            temp_array(j) = scoll(j, i, k)
         end do

         call calculate_statistics(temp_array, collv, collav, mc_steps_outer)

         rel(k) = 0
         if (collav /= 0) rel(k) = collv / collav
         coll(i, k) = species_concentration(i) * collav
      end do

      write(unit_output, 2013) i, (coll(i, k), rel(k), k = 1, num_species)
   end do

   write(unit_output, '(/)')

   return

end subroutine calculate_collision_pressure


! =============================================================================
! =============================================================================
!
! Subroutine: widom
!
! Puts a ghost particle anywhere in the cell and calculates the
! energy with the rest of the system. Averaging of the exponent
! results in the chemical potential. A rescaling of the ion-
! charges is neccessary to assure neutrality. The integration
! is performed in entry 2 with a simpson's rule.
! The total chemical potential is divided in to the ideal (chid),
! hard core (chhc) and the electric contribution (chel). For
! each species the chem.pot. is given cellaverage, at the wall,
! and in the midplan.
!
! =============================================================================
! =============================================================================

subroutine calculate_widom_insertion

   implicit none
   real :: ran2
   include 'bulk_f90.inc'

   ! Local variables
   ! max_widom_species accounts for species at multiple measurement locations (bulk, wall, midplane)
   integer, parameter :: max_widom_species = max_species * 3  ! max_species * (max_measurement_locations + 1)
   integer :: i, j, k, mp, ntel, ntocp, nwtot, nh, nl, jm, k1, irsum
   integer :: ihc(max_widom_species), irej(max_widom_species), mwcn(25), ihcall(0:5)
   double precision :: chel(25, max_widom_species), chhc(25, max_widom_species), chex(25, max_widom_species)
   double precision :: chto(25, max_widom_species), chexw(25, max_widom_species), dch1(25, max_widom_species)
   double precision :: dch2(25, max_widom_species), chint(max_widom_species, 11), ewnom(max_widom_species, 11)
   double precision :: ewden(max_widom_species, 11), chinta(max_widom_species, 11)
   double precision :: expuw(max_widom_species), chid(max_widom_species)
   double precision :: chelav(max_widom_species), chelv(max_widom_species), chhcav(max_widom_species), chhcv(max_widom_species)
   double precision :: chexav(max_widom_species), chexv(max_widom_species), chtoav(max_widom_species), chtov(max_widom_species)
   double precision :: chexwa(max_widom_species), chexwv(max_widom_species)
   double precision :: ddx, ddy, ddz, x, y, z, ew, ewd, ewla, wtot2, wtot3, aint1, aint2, aint4, wsum
   double precision :: uj1, pcollav, pcollv, cwiav, cwiv
   character(len=80) :: str(max_species)

   ! Note: Variable names have been updated to be more descriptive
   ! See bulk_f90.inc for the complete variable declarations

   do i = 0, measurement_location
      do j = 1, num_species
         if (species_concentration(j) /= 0) then
            chid(j + num_species* i) = dlog( species_concentration(j) / avogadro_number* 1.0d27)
         else
            chid(j + num_species* i) = -77
         end if
      end do
   end do

   do j = 1, max_species
      do i = 1, 25
         chel(i, j) = 0
         chhc(i, j) = 0
         chex(i, j) = 0
         chto(i, j) = 0
      end do

      chelav(j) = 0
      chhcav(j) = 0
      chexav(j) = 0
      chtoav(j) = 0
      expuw(j)  = 0
      ihc(j)    = 0

      do k = 1, 11
         ewden(j, k)  = 0
         ewnom(j, k)  = 0
         chinta(j, k) = 0
      end do
   end do

   do j = 1, 25
      mwcn(j) = 0
   end do

   do j = 0, 5
      ihcall(j) = 0
   end do

   ntel = 0

   str(1) = 'Chemical potentials averaged over the whole cell'
   str(2) = 'Chemical potentials at the wall'
   str(3) = 'Chemical potentials at the midplane'

   return


   entry calculate_widom_insertion1

   mwcn(current_macro_step) = mwcn(current_macro_step) + 1

   do mp = 0, measurement_location
      do i = 1, num_widom_insertions
         x = box_size* (ran2(random_seed) - 0.5)
         y = box_size* (ran2(random_seed) - 0.5)
         z = 0

         if (mp <= 1) then
            z = box_size* (ran2(random_seed) - 0.5)
            if (mp == 1) z = sign(box_half, z)
         end if

         do j = 1, num_particles
            ddx = abs(x - particle_x(j))
            if (ddx > box_half) ddx = ddx - box_size
            ddy = abs(y - particle_y(j))
            if (ddy > box_half) ddy = ddy - box_size
            ddz = abs(z - particle_z(j))
            if (ddz > box_half) ddz = ddz - box_size
            distance_squared(j) = ddx * ddx + ddy * ddy + ddz * ddz
         end do

         irsum = 0

         do j = 1, num_species
            jm = mp * num_species+ j
            irej(jm) = 0

            do k = 1, num_particles
               if (distance_squared(k) < hard_core_distance_sq(k, j)) irej(jm) = 1
            end do

            irsum = irsum + irej(jm)
         end do

         if (irsum == num_species) then
            ihcall(mp) = ihcall(mp) + 1
            go to 110
         end if

         do k = 1, num_particles
            distance_inverse(k) = 1.0 / sqrt(distance_squared(k))
         end do

         wtot2 = 0
         wtot3 = 0
         nl    = 1

         do j = 1, num_species
            nh   = int(species_properties(j, 2))
            wsum = 0

            do k = nl, nl + nh - 1
               wsum = wsum + distance_inverse(k)
            end do

            nl = nl + nh
            wtot2 = wtot2 + wsum * species_properties(j, 4)
            wtot3 = wtot3 + wsum
         end do

         wtot2 = wtot2 + uj1

         do j = 1, num_species
            jm = mp * num_species+ j

            if (irej(jm) == 1) then
               ihc(jm) = ihc(jm) + 1
               go to 160
            end if

            expuw(jm) = expuw(jm) + exp(beta_inverse_temp* wtot2 * species_properties(j, 4) * energy_conversion_factor)

            do k1 = 0, 10
               k   = k1 + 1
               ew  = species_properties(j, 4) * (wtot2 - k1 * 0.1 * species_properties(j, 4) * wtot3 / num_particles)
               ewla = ew * k1 * 0.1
               ewd  = exp(beta_inverse_temp* energy_conversion_factor* ewla)
               ewden(jm, k) = ewden(jm, k) + ewd
               ewnom(jm, k) = ewnom(jm, k) - ew * beta_inverse_temp* energy_conversion_factor* ewd
            end do

160         continue
         end do

110      continue
      end do
   end do

   return


   entry calculate_widom_insertion2

   ntocp = num_species* (measurement_location+ 1)

   do i = 1, ntocp
      do j = 1, 11
         if (ewden(i, j) == 0) then
            write(unit_output, *) ' WIDOM DENOMINATOR EQUALS ZERO', i, j
         else
            chint(i, j)  = ewnom(i, j) / ewden(i, j)
            ewnom(i, j)  = 0
            ewden(i, j)  = 0
            chinta(i, j) = chinta(i, j) + 1.0 / mc_steps_outer * chint(i, j)
         end if
      end do

      aint4 = chint(i, 2) + chint(i, 4) + chint(i, 6) + chint(i, 8) + chint(i, 10)
      aint2 = chint(i, 3) + chint(i, 5) + chint(i, 7) + chint(i, 9)
      aint1 = chint(i, 1) + chint(i, 11)
      chel(current_macro_step, i) = 1.0 / 30.0 * (aint1 + 2 * aint2 + 4 * aint4)
   end do

   nwtot = mwcn(current_macro_step) * num_widom_insertions
   nwtot = num_widom_insertions* int(mc_steps_inner / widom_interval) * mc_steps_middle

   do i = 1, ntocp
      ihc(i) = ihc(i) + ihcall(int((i - 1) / num_species))
      chhc(current_macro_step, i)  = -dlog(dble(nwtot - ihc(i)) / nwtot)
      chexw(current_macro_step, i) = -dlog(expuw(i) / nwtot)
      chex(current_macro_step, i)  = chel(current_macro_step, i) + chhc(current_macro_step, i)
      chto(current_macro_step, i)  = chex(current_macro_step, i) + chid(i)
      expuw(i) = 0
      ihc(i)   = 0
   end do

   do j = 1, num_species
      dch1(current_macro_step, j) = chto(current_macro_step, j) - chex(current_macro_step, num_species+ j)
      dch2(current_macro_step, j) = chto(current_macro_step, j) - chex(current_macro_step, 2 * num_species+ j)
   end do

   write(unit_macro, 2010) nwtot

   do j = 0, measurement_location
      ihcall(j) = 0
      write(unit_macro, '(/, a)') str(j + 1)
      write(unit_macro, 2015)
      write(unit_macro, 2020) (i, chid(i), chhc(current_macro_step, i), chel(current_macro_step, i), &
                        chex(current_macro_step, i), chto(current_macro_step, i), chexw(current_macro_step, i), &
                        i = j * num_species+ 1, (j + 1) * num_species)
   end do

   return


   entry calculate_widom_insertion3

   nwtot = 0

   do i = 1, mc_steps_outer
      nwtot = nwtot + mwcn(i) * num_widom_insertions
   end do

   ntocp = num_species* (measurement_location+ 1)
   i = 1

   call calculate_statistics_per_species(chel,  chelv,  chelav,  mc_steps_outer, ntocp)
   call calculate_statistics_per_species(chexw, chexwv, chexwa,  mc_steps_outer, ntocp)
   call calculate_statistics_per_species(chhc,  chhcv,  chhcav,  mc_steps_outer, ntocp)
   call calculate_statistics_per_species(chex,  chexv,  chexav,  mc_steps_outer, ntocp)
   call calculate_statistics_per_species(chto,  chtov,  chtoav,  mc_steps_outer, ntocp)

   write (unit_output, '(a, /)') 'CONTACT CORRELATION g(r)'
   write (unit_output, 2001) (i, i = 1, num_species)

   pressure_collision_average= 0.0
   pressure_collision_variance = 0.0

   do i = 1, num_species
      do k = 1, num_species
         do j = 1, mc_steps_outer
            temp_array(j) = contact_correlation(j, i, k) * exp(chexw(j, k))
         end do

         call calculate_statistics(temp_array, cwiv, cwiav, mc_steps_outer)
         contact_correlation(11, i, k) = cwiv
         contact_correlation(12, i, k) = cwiav
      end do

      write(unit_output, 2002) i, (contact_correlation(12, i, k), contact_correlation(11, i, k), k = 1, num_species)
   end do

   write (unit_output, '(/, a, /)') 'CONTACT PRESSURE MATRIX'
   write (unit_output, 2001) (i, i = 1, num_species)

   do i = 1, num_species
      do k = 1, num_species
         do j = 1, mc_steps_outer
            temp_array(j) = 2.0 / 3.0 * pi * species_concentration(i) * contact_correlation(j, i, k) * &
                     exp(chexw(j, k) + chid(k)) * (species_properties(i, 3) + species_properties(k, 3)) ** 3
         end do

         call calculate_statistics(temp_array, cwiv, cwiav, mc_steps_outer)
         contact_correlation(11, i, k) = cwiv
         contact_correlation(12, i, k) = cwiav
         pressure_collision_average= pressure_collision_average+ cwiav
         pressure_collision_variance = pressure_collision_variance + cwiv * cwiv
      end do

      write(unit_output, 2002) i, (contact_correlation(12, i, k), contact_correlation(11, i, k), k = 1, num_species)
   end do

   pressure_collision_variance= sqrt(pcollv)

   write(unit_output, '(/, a, f12.5, f10.5, /)') 'Total Collision pressure   =        &
                                          &', pcollav, pcollv


   write(unit_output, 44)
   write(unit_output, 2034)
   write(unit_output, 2035) ((i - 1.0) / 10.0, i = 1, 11)

   do i = 1, ntocp
      write(unit_output, 2040) i, (chinta(i, j), j = 1, 11)
   end do

   write(unit_output, 2010) nwtot

   do j = 0, measurement_location
      write(unit_output, '(/, a)') str(j + 1)
      write(unit_output, 2015)
      write(unit_output, 2030) (mod(i - 1, num_species) + 1, chid(i), chhcav(i), chhcv(i) &
                        , chelav(i), chelv(i), chexav(i), chexv(i), chtoav(i), chtov(i) &
                        , chexwa(i), chexwv(i), i = j * num_species+ 1, num_species* (j + 1))
      write(unit_output, 2033) ( exp((2 * chexav(i - 2) + chexav(i - 1)) / 3) )
   end do


2001 format('Species      ', 10(i12, 12x))
2002 format(i4, '        ', 10(f12.5, f10.5))
2010 format(/, /, 'CHEMICAL POTENTIALS IN UNITS OF KT, MEASURED', i8, ' TIMES')
2015 format(/, 'Type  Ideal     Hard Core       El. Static       Excess', &
           '         Total         Uncorrected Excess')
2020 format(10(i3, f8.4, 5(f9.4, '       '), /))
2030 format(10(i3, f8.4, 5(f9.4, f7.4), /))
2033 format('2:1 ', 10(f8.4), /)
2034 format(/, /, 'INTEGRAND IN WIDOM CORRECTION FOR DIFFERENT CHARGE PARAMETERS')
2035 format('TYPE', 11('  ', f4.2, '  '), /)
2040 format(10(i3, 11f8.3), /)
44   format(/, '************************************************************', &
           '**********************************************')

   return

end subroutine calculate_widom_insertion
