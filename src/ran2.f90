! =============================================================================
! Random Number Generator: Long Period Uniform Distribution
!
! This is the ran2 algorithm from Numerical Recipes, providing a long period
! (> 2×10¹⁸) random number generator using Schrage's method to avoid overflow
! in modular arithmetic.
!
! The algorithm combines two linear congruential generators and uses a shuffle
! table to break up sequential correlations.
! =============================================================================

function ran2(random_seed)

    implicit none

    ! Arguments
    integer, intent(inout) :: random_seed      ! Seed for random number generator (modified on each call)
    real                   :: ran2             ! Returns uniform random number in (0,1)

    ! Linear congruential generator 1 parameters
    integer, parameter :: modulus_1         = 2147483563  ! IM1: First modulus
    integer, parameter :: multiplier_1      = 40014       ! IA1: First multiplier
    integer, parameter :: quotient_1        = 53668       ! IQ1: modulus_1 / multiplier_1
    integer, parameter :: remainder_1       = 12211       ! IR1: modulus_1 mod multiplier_1

    ! Linear congruential generator 2 parameters
    integer, parameter :: modulus_2         = 2147483399  ! IM2: Second modulus
    integer, parameter :: multiplier_2      = 40692       ! IA2: Second multiplier
    integer, parameter :: quotient_2        = 52774       ! IQ2: modulus_2 / multiplier_2
    integer, parameter :: remainder_2       = 3791        ! IR2: modulus_2 mod multiplier_2

    ! Shuffle table parameters
    integer, parameter :: shuffle_table_size = 32            ! NTAB: Size of shuffle table
    integer, parameter :: modulus_1_minus_1  = modulus_1 - 1 ! IMM1
    integer, parameter :: table_divisor      = 1 + modulus_1_minus_1 / shuffle_table_size  ! NDIV

    ! Scaling and safety parameters
    real, parameter :: scaling_factor     = 1.0 / modulus_1  ! AM: Converts to (0,1)
    real, parameter :: epsilon            = 1.2e-7           ! EPS: Small value to avoid exact endpoints
    real, parameter :: random_max         = 1.0 - epsilon    ! RNMX: Maximum return value

    ! Saved state variables (persistent across function calls)
    integer :: generator_2_state                      ! State of second generator
    integer :: shuffle_table(shuffle_table_size)      ! Shuffle table for decorrelation
    integer :: previous_output                        ! Previous output from shuffle table

    save shuffle_table, previous_output, generator_2_state

    ! Initial values for saved variables
    data generator_2_state  / 123456789 /
    data shuffle_table      / shuffle_table_size * 0 /
    data previous_output    / 0 /

    ! Local variables
    integer :: table_index  ! Index into shuffle table
    integer :: temp_value   ! Temporary value for Schrage's method


    ! =========================================================================
    ! Initialization: First call or when seed is reset (negative value)
    ! =========================================================================
    if (random_seed <= 0) then
        ! Ensure seed is positive
        random_seed = max(-random_seed, 1)
        generator_2_state = random_seed

        ! Fill shuffle table using generator 1
        ! Extra iterations (table_size + 8) ensure table is well shuffled
        do table_index = shuffle_table_size + 8, 1, -1
            ! Schrage's method: compute (multiplier_1 * seed) mod modulus_1
            ! without overflow by using seed = quotient * Q + remainder
            temp_value  = random_seed / quotient_1
            random_seed = multiplier_1 * (random_seed - temp_value * quotient_1) - temp_value * remainder_1

            ! Ensure result is positive
            if (random_seed < 0) random_seed = random_seed + modulus_1

            ! Store in shuffle table (only the first shuffle_table_size values)
            if (table_index <= shuffle_table_size) shuffle_table(table_index) = random_seed
        end do

        ! Initialize previous output with first table entry
        previous_output = shuffle_table(1)
    end if


    ! =========================================================================
    ! Generate random number using both generators and shuffle table
    ! =========================================================================

    ! Update generator 1 using Schrage's method
    temp_value  = random_seed / quotient_1
    random_seed = multiplier_1 * (random_seed - temp_value * quotient_1) - temp_value * remainder_1
    if (random_seed < 0) random_seed = random_seed + modulus_1

    ! Update generator 2 using Schrage's method
    temp_value        = generator_2_state / quotient_2
    generator_2_state = multiplier_2 * (generator_2_state - temp_value * quotient_2) - temp_value * remainder_2
    if (generator_2_state < 0) generator_2_state = generator_2_state + modulus_2

    ! Use previous output to select index in shuffle table
    ! This breaks up sequential correlations
    table_index = 1 + previous_output / table_divisor

    ! Combine generator 2 output with shuffled generator 1 output
    previous_output            = shuffle_table(table_index) - generator_2_state
    shuffle_table(table_index) = random_seed

    ! Ensure output is positive
    if (previous_output < 1) previous_output = previous_output + modulus_1_minus_1

    ! Scale to (0,1) and ensure we never return exact 1.0
    ran2 = min(scaling_factor * previous_output, random_max)

    return

end function ran2
