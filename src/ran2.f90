! =============================================================================
! Random number generator: ran2
! Long period (> 2×10^18) random number generator
! =============================================================================

function ran2(idum)

    implicit none

    ! Arguments
    integer, intent(inout) :: idum
    real                   :: ran2

    ! Constants
    integer, parameter :: IM1   = 2147483563
    integer, parameter :: IM2   = 2147483399
    integer, parameter :: IMM1  = IM1 - 1
    integer, parameter :: IA1   = 40014
    integer, parameter :: IA2   = 40692
    integer, parameter :: IQ1   = 53668
    integer, parameter :: IQ2   = 52774
    integer, parameter :: IR1   = 12211
    integer, parameter :: IR2   = 3791
    integer, parameter :: NTAB  = 32
    integer, parameter :: NDIV  = 1 + IMM1 / NTAB

    real, parameter :: AM   = 1.0 / IM1
    real, parameter :: EPS  = 1.2e-7
    real, parameter :: RNMX = 1.0 - EPS

    ! Local variables
    integer            :: idum2, j, k
    integer            :: iv(NTAB), iy

    save iv, iy, idum2

    data idum2 / 123456789 /
    data iv    / NTAB*0 /
    data iy    / 0 /

    ! Initialize
    if (idum <= 0) then
        idum  = max(-idum, 1)
        idum2 = idum

        do j = NTAB + 8, 1, -1
            k    = idum / IQ1
            idum = IA1 * (idum - k * IQ1) - k * IR1

            if (idum < 0) idum = idum + IM1
            if (j <= NTAB) iv(j) = idum
        end do

        iy = iv(1)
    end if

    ! Generate random number
    k    = idum / IQ1
    idum = IA1 * (idum - k * IQ1) - k * IR1
    if (idum < 0) idum = idum + IM1

    k     = idum2 / IQ2
    idum2 = IA2 * (idum2 - k * IQ2) - k * IR2
    if (idum2 < 0) idum2 = idum2 + IM2

    j     = 1 + iy / NDIV
    iy    = iv(j) - idum2
    iv(j) = idum

    if (iy < 1) iy = iy + IMM1

    ran2 = min(AM * iy, RNMX)

    return

end function ran2
