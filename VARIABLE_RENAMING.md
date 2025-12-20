# Variable Renaming Map for F90 Refactoring

## Parameters
- `mion` → `max_ions` - Maximum number of ions
- `mxspec` → `max_species` - Maximum number of species

## Physical Constants
- `pi` → `pi` - Mathematical constant (keep)
- `avno` → `avogadro_number` - Avogadro's number
- `bk` → `gas_constant` - Gas constant (R)
- `eps` → `dielectric_constant` - Relative dielectric constant
- `epsx` → `vacuum_permittivity` - Vacuum permittivity (ε₀)
- `ecf` → `energy_conversion_factor` - Energy conversion factor
- `ech` → `elementary_charge` - Elementary charge
- `abeta` → `beta_inverse_temp` - -1/(kT) inverse temperature parameter
- `dtemp` → `temperature` - Temperature

## Particle Arrays
- `x6` → `particle_x` - X coordinates of particles
- `y6` → `particle_y` - Y coordinates of particles
- `z6` → `particle_z` - Z coordinates of particles
- `chv` → `particle_charge` - Charge of each particle
- `dp` → `displacement_max` - Maximum displacement for MC moves
- `ispc` → `particle_species` - Species index for each particle

## Energy Arrays
- `esa` → `energy_matrix` - Pair interaction energy matrix
- `ei` → `trial_energy` - Energy for trial configuration
- `rw2` → `distance_squared` - Squared distances
- `rwi` → `distance_inverse` - Inverse distances (1/r)

## Species Properties
- `hion` → `species_properties` - Properties of each species (number, radius, charge, displacement)
- `hc2v` → `hard_core_distance_sq` - Squared hard core distances
- `caver` → `species_concentration` - Average concentration per species
- `hmik` → `species_interaction_matrix` - Species interaction matrix
- `nspec` → `num_species` - Number of species
- `ispec` → `current_species` - Current species index

## Trial Move Variables
- `tx6` → `trial_x` - Trial X position
- `ty6` → `trial_y` - Trial Y position
- `tz6` → `trial_z` - Trial Z position
- `il` → `current_particle` - Index of current particle
- `isos` → `overlap_status` - Overlap status flag

## Simulation Box
- `box` → `box_size` - Simulation box size
- `box2` → `box_half` - Half box size
- `box2i` → `box_half_inverse` - Inverse of half box size

## Monte Carlo Loop Counters
- `ny1` → `mc_steps_inner` - Inner MC loop steps
- `ny2` → `mc_steps_middle` - Middle MC loop steps
- `ny3` → `mc_steps_outer` - Outer MC loop steps (macro steps)
- `my3` → `current_macro_step` - Current macro step counter
- `npart` → `num_particles` - Total number of particles

## Widom Insertion Variables
- `nwins` → `num_widom_insertions` - Number of Widom particle insertions
- `nwint` → `widom_interval` - Interval between Widom calculations
- `nfix` → `measurement_location` - Location for chemical potential measurement
- `cwi` → `contact_correlation` - Contact correlation function

## Pressure Variables
- `pwiwal` → `pressure_widom_wall` - Widom pressure at wall
- `pwimid` → `pressure_widom_midplane` - Widom pressure at midplane
- `pwiwalv` → `pressure_widom_wall_variance` - Variance of wall pressure
- `pwimidv` → `pressure_widom_midplane_variance` - Variance of midplane pressure
- `pcollav` → `pressure_collision_average` - Average collision pressure
- `pcollv` → `pressure_collision_variance` - Collision pressure variance

## Energy Tracking
- `xww1` → `total_coulomb_energy` - Total Coulomb energy
- `qfww1` → `energy_check_parameter` - Energy check parameter

## Random Number Generator
- `islu` → `random_seed` - Random number generator seed
- `ran2` → `random_uniform` - Function name for RNG

## Configuration Type
- `ink` → `initial_config_type` - Initial configuration type (0=file, 1=random)

## File Units
- `lll` → `unit_conf` - Configuration file unit
- `mmm` → `unit_input` - Input file unit
- `jjj` → `unit_output` - Output file unit
- `iii` → `unit_misc` - Miscellaneous file unit
- `kkk` → `unit_macro` - Macro output file unit

## Temporary/Helper Variables
- `dum` → `temp_array` - Temporary array for calculations
