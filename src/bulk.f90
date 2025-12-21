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

   ! Note: Variable names have been updated to be more descriptive
   ! See bulk_f90.inc for the complete variable declarations

   ! Local variables for Monte Carlo simulation
   integer :: step_inner, step_middle         ! Loop counters for inner and middle MC loops
   integer :: kn1, kn2, key, kntot            ! Counters for accepted/rejected moves
   integer :: i, j, k, l, num                 ! Loop counters and temporary variables
   integer :: iclk, lclk, nkonf, nwtot        ! Clock and configuration counters
   integer :: mtot(max_species), macc(max_species), menrj(max_species), mhcrj(max_species)

   double precision :: uula(2), uuta(2), suu(25), tijd(10)
   double precision :: deltae, dzxp, uuvar     ! Energy and variance variables
   double precision :: t                       ! Temperature variable
   double precision :: pexen, pexenv, pid      ! Pressure variables
   double precision :: ptot, ptotv             ! Total pressure variables
   double precision :: qww1, yn3, ytot1        ! Energy check and averaging variables

   data uuta  / 2*0.0 /
   data suu   / 25*0.0 /
   data mtot, macc, menrj, mhcrj / 40*0 /
   data kn1, kn2, key, kntot / 4*0 /
   data tijd  / 10*0.0 /

43 format( 'MONTE CARLO SIMULATION OF MULTICOMPONENT ISOTROPIC    ', &
            'IONIC SYSTEM  ', /, /, 'Written by Peter Bolhuis, February', &
            '1992', / )

44 format( '************************************************************', &
            '************************', / )

   iclk = 0
   lclk = 0

   unit_macro = 12
   unit_conf = 13
   unit_input = 15
   unit_output = 6

   write(unit_output, 44)
   write(unit_output, 43)
   write(unit_output, 44)


   ! ==========================================================================
   ! Physical constants and conversion factors
   ! ==========================================================================

   pi   = acos(-1.0d0)
   vacuum_permittivity= 8.85418782d-12
   elementary_charge = 1.60219d-19
   avogadro_number= 6.0223d23
   gas_constant  = 8.314


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
   read(unit_input, *) (species_properties(l, 1), species_properties(l, 2), species_properties(l, 3), species_properties(l, 4), species_properties(l, 5), l = 1, num_species)
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

   do i = 1, num_species
      num = int(species_properties(i, 2))

      do j = 1, num
         num_particles= num_particles+ 1
         particle_species(num_particles) = i

         do k = 1, num_species
            hard_core_distance_sq(num_particles, k) = (species_properties(i, 3) + species_properties(k, 3)) ** 2
         end do

         particle_charge(num_particles) = species_properties(i, 4)
         displacement_max(num_particles)  = species_properties(i, 5)
      end do
   end do

   do i = 1, num_species
      species_concentration(i) = (species_properties(i, 2)) / (box_size** 3)
   end do

   do i = 1, num_particles
      trial_energy(i) = 0.0
   end do

   nkonf = mc_steps_inner * mc_steps_middle * mc_steps_outer

   if (mc_steps_outer > 25) mc_steps_outer = 25

   t     = temperature
   beta_inverse_temp= -1.0 / (gas_constant * temperature)


   ! ==========================================================================
   ! Read configuration from file
   ! ==========================================================================

   if (initial_config_type== 0) then
      read(unit_conf, *) (particle_x(l), particle_y(l), particle_z(l), l = 1, num_particles)
      write(unit_output, *) 'Configuration read from file'
   end if

   if (initial_config_type== 1) call slump
   call collision
   call widom


   ! ==========================================================================
   ! Write input data
   ! ==========================================================================

   nwtot = num_widom_insertions* (int(mc_steps_inner / widom_interval) * mc_steps_middle * mc_steps_outer)

   write(unit_output, 800) num_particles
   write(unit_output, 801) (l, int(species_properties(l, 2)), species_properties(l, 3), species_properties(l, 4), &
                    species_properties(l, 5), species_concentration(l) * 1.0d27 / avogadro_number, l = 1, num_species)
   write(unit_output, 802) dielectric_constant, temperature, box_size, mc_steps_inner, mc_steps_middle, mc_steps_outer, nkonf, &
                   dble(nkonf / num_particles), nwtot

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

   call liv

   ytot1 = total_coulomb_energy* energy_conversion_factor* 1.0d-3

   write(unit_output, '(a)') 'Initial configuration'
   write(unit_output, '(a, e12.5)') 'Coulomb energy (kJ/mol)  =', ytot1


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
      uula(1) = 0.0

      do step_middle = 1, mc_steps_middle
         do step_inner = 1, mc_steps_inner
            current_particle   = 1 + mod(kntot, num_particles)
            current_species= particle_species(current_particle)
            kntot = kntot + 1
            mtot(current_species) = mtot(current_species) + 1

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
            call qin

            if (overlap_status/= 0) go to 63

            deltae = 0.0
            do j = 1, num_particles
               deltae = deltae + trial_energy(j) - energy_matrix(j, current_particle)
            end do

            dzxp = beta_inverse_temp* (deltae * energy_conversion_factor)

            if (dzxp < -80.0) go to 64
            if (dzxp >   0.0) go to 62
            if (exp(dzxp) < ran2(random_seed)) go to 64

62          macc(current_species) = macc(current_species) + 1

            ! Trial config accepted
            particle_x(current_particle) = trial_x
            particle_y(current_particle) = trial_y
            particle_z(current_particle) = trial_z

            do j = 1, num_particles
               energy_matrix(current_particle, j) = trial_energy(j)
               energy_matrix(j, current_particle) = trial_energy(j)
            end do

            total_coulomb_energy= total_coulomb_energy+ deltae
            go to 65

63          mhcrj(current_species) = mhcrj(current_species) + 1
            go to 65

64          menrj(current_species) = menrj(current_species) + 1

65          uula(1) = uula(1) + total_coulomb_energy

            if (mod(step_inner, widom_interval) == 0) then
               call collision1
               call widom1
            end if

         end do
      end do

      qww1 = total_coulomb_energy
      call liv

      if (total_coulomb_energy/= 0) qww1 = (qww1 - total_coulomb_energy) / total_coulomb_energy

      write(unit_macro, *)
      write(unit_macro, 44)
      write(unit_macro, 73) current_macro_step, qww1, energy_check_parameter

      if (abs(qww1) > 0.001) stop

      call liv

      if (widom_interval<= mc_steps_inner) then
         call collision2
         call widom2
      end if

      uula(1) = uula(1) / dble(mc_steps_inner * mc_steps_middle)
      uuta(1) = uuta(1) + uula(1)
      uula(1) = uula(1) * energy_conversion_factor* 1.0d-3
      suu(current_macro_step) = uula(1)

      ! Calculate total acceptance and rejection
      key = 0
      kn1 = 0
      kn2 = 0

      do i = 1, num_species
         key = key + macc(i)
         kn1 = kn1 + mhcrj(i)
         kn2 = kn2 + menrj(i)
      end do

      write(unit_macro, 810) kntot, key, kn2, kn1
      write(unit_macro, '(a, e12.5)') 'Coulomb energy (kj/mol)  =', uula(1)
   end do

73  format(/ ' MACROSTEP  ', i10, /, /, 'Checkparameters : ', 3e10.3)
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

   yn3     = 1.0 / float(mc_steps_outer)
   uuta(1) = yn3 * uuta(1)
   uuta(1) = uuta(1) * energy_conversion_factor* 1.0d-3

   call earth(suu, uuvar, uuta(1), mc_steps_outer)

   write(unit_output, '(/)')
   write(unit_output, 44)
   write(unit_output, '(/, a, i3, a, /)') 'FINAL RESULTS AFTER ', mc_steps_outer, &
                                  ' MACROSTEPS'
   write(unit_output, 810) kntot, key, kn2, kn1
   write(unit_output, 811) (l, dble(macc(l)) / mtot(l), dble(menrj(l)) &
                    / mtot(l), dble(mhcrj(l)) / mtot(l), l = 1, num_species)
   write(unit_output, '(a, 2e12.5)') 'Coulomb energy (kj/mol)  =', uuta(1), uuvar
   write(unit_output, *)

   if (widom_interval<= mc_steps_inner) then
      call collision3
      call widom3
   end if

   rewind unit_conf
   write(unit_conf, 771) (particle_x(l), particle_y(l), particle_z(l), l = 1, num_particles)
771 format(5e16.8)

   close(unit_conf)


   ! ==========================================================================
   !
   ! Calculation of pressures from the densities
   !
   ! ==========================================================================

   pid = 0.0
   do k = 1, num_species
      pid = pid + species_concentration(k) * 1.0d27 / avogadro_number
   end do

   pexen  = uuta(1) * 1.0d30 / (3 * gas_constant* temperature* box_size** 3 * avogadro_number)
   pexenv = uuvar   * 1.0d30 / (3 * gas_constant* temperature* box_size** 3 * avogadro_number)
   ptot   = pid + pressure_collision_average+ pexen
   ptotv  = sqrt(pressure_collision_variance* pressure_collision_variance+ pexenv * pexenv)

   write(unit_output, '(/, a, /)') 'TOTAL BULK PRESSURE '
   write(unit_output, 729) 'Ideal pressure       ', pid, 0.0
   write(unit_output, 729) 'Energy con. <E/3V>   ', pexen, pexenv
   write(unit_output, 729) 'Collision press      ', pressure_collision_average, pressure_collision_variance
   write(unit_output, 731)
   write(unit_output, 729) 'Bulk pressure        ' , ptot, ptotv

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

subroutine earth(a, bbbwww, xb, nnn)

   implicit none
   integer :: k, nnn
   double precision :: a(25), bbbwww, xb
   double precision :: b, yak

   yak = 1.0 / (nnn * (nnn - 1))
   b   = 0.0
   xb  = 0.0

   do k = 1, nnn
      xb = xb + a(k)
   end do

   xb = xb / nnn

   do k = 1, nnn
      b = b + (a(k) - xb) * (a(k) - xb)
   end do

   b = sqrt(yak * b)
   bbbwww = b

   return

end subroutine earth


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

subroutine earth2(a, bbbwww, xb, nnn, num)

   implicit none
   integer, parameter :: max_species = 10
   integer :: i, j, k, nnn, num
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

end subroutine earth2


! =============================================================================
! =============================================================================
!
! Subroutine: slump
!
! Set up of the random initial configuration.
!
! =============================================================================
! =============================================================================

subroutine slump

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

end subroutine slump


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

subroutine qin

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

end subroutine qin


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

subroutine liv

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

end subroutine liv


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

subroutine collision

   implicit none
   real :: ran2
   include 'bulk_f90.inc'

   ! Local variables
   integer :: i, j, k, isp, ksp, nisp, nnn, num, nwtot
   double precision :: rel(max_species)
   double precision :: scoll(25, max_species, max_species), coll(max_species, max_species)
   double precision :: ddx, ddy, ddz, aa, dis, v, g2, wx6, wy6, wz6, urej, wtot2
   double precision :: collav, collv, cwiav, cwiv

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

   entry collision1

   do i = 1, num_widom_insertions
      num = 1

      do isp = 1, num_species
         nisp = int(species_properties(isp, 2))

         do ksp = 1, num_species
            nnn = ran2(random_seed) * nisp + num
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

   entry collision2

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
2011 format(/, i4, '    Col. samples', 10('            ', i6))
2012 format('        pressure    ', 10e18.6, /)
2013 format(i4, '        ', 20e12.5)


   ! -------------------------------------------------------------------------

   entry collision3

   write(unit_output, '(/, a, /)') 'COLLISION MATRIX  AND RELATIVE ERROR'
   write (unit_output, 2010) (i, i = 1, num_species)

   do i = 1, num_species
      do k = 1, num_species
         collav = 0

         do j = 1, mc_steps_outer
            contact_correlation(j, i, k) = scoll(j, i, k)
            temp_array(j) = scoll(j, i, k)
         end do

         call earth(temp_array, collv, collav, mc_steps_outer)

         rel(k) = 0
         if (collav /= 0) rel(k) = collv / collav
         coll(i, k) = species_concentration(i) * collav
      end do

      write(unit_output, 2013) i, (coll(i, k), rel(k), k = 1, num_species)
   end do

   write(unit_output, '(/)')

   return

end subroutine collision


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

subroutine widom

   implicit none
   real :: ran2
   include 'bulk_f90.inc'

   ! Local variables
   integer :: i, j, k, mp, ntel, ntocp, nwtot, nh, nl, jm, k1, irsum
   integer :: ihc(max_species), irej(max_species), mwcn(25), ihcall(0:5)
   double precision :: chel(25, max_species), chhc(25, max_species), chex(25, max_species)
   double precision :: chto(25, max_species), chexw(25, max_species), dch1(25, max_species)
   double precision :: dch2(25, max_species), chint(max_species, 11), ewnom(max_species, 11)
   double precision :: ewden(max_species, 11), chinta(max_species, 11)
   double precision :: expuw(max_species), chid(max_species)
   double precision :: chelav(max_species), chelv(max_species), chhcav(max_species), chhcv(max_species)
   double precision :: chexav(max_species), chexv(max_species), chtoav(max_species), chtov(max_species)
   double precision :: chexwa(max_species), chexwv(max_species)
   double precision :: ddx, ddy, ddz, x, y, z, ew, ewd, ewla, wtot2, wtot3, aint1, aint2, aint4, wsum
   double precision :: uj1, pcollav, pcollv, cwiav, cwiv
   character*80 str(max_species)

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


   entry widom1

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


   entry widom2

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


   entry widom3

   nwtot = 0

   do i = 1, mc_steps_outer
      nwtot = nwtot + mwcn(i) * num_widom_insertions
   end do

   ntocp = num_species* (measurement_location+ 1)
   i = 1

   call earth2(chel,  chelv,  chelav,  mc_steps_outer, ntocp)
   call earth2(chexw, chexwv, chexwa,  mc_steps_outer, ntocp)
   call earth2(chhc,  chhcv,  chhcav,  mc_steps_outer, ntocp)
   call earth2(chex,  chexv,  chexav,  mc_steps_outer, ntocp)
   call earth2(chto,  chtov,  chtoav,  mc_steps_outer, ntocp)

   write (unit_output, '(a, /)') 'CONTACT CORRELATION g(r)'
   write (unit_output, 2001) (i, i = 1, num_species)

   pressure_collision_average= 0.0
   pressure_collision_variance = 0.0

   do i = 1, num_species
      do k = 1, num_species
         do j = 1, mc_steps_outer
            temp_array(j) = contact_correlation(j, i, k) * exp(chexw(j, k))
         end do

         call earth(temp_array, cwiv, cwiav, mc_steps_outer)
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

         call earth(temp_array, cwiv, cwiav, mc_steps_outer)
         contact_correlation(11, i, k) = cwiv
         contact_correlation(12, i, k) = cwiav
         pressure_collision_average= pressure_collision_average+ cwiav
         pressure_collision_variance = pressure_collision_variance + cwiv * cwiv
      end do

      write(unit_output, 2002) i, (contact_correlation(12, i, k), contact_correlation(11, i, k), k = 1, num_species)
   end do

   pressure_collision_variance= sqrt(pcollv)

   write(unit_output, '(/, a, f12.5, f10.5, /)') 'Total Collision pressure   =        &
                                          ', pcollav, pcollv


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

end subroutine widom
