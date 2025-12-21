# Software archeology from Division of Theoretical Chemistry, Lund University 

Collection of old fortran codes for computational chemistry.

## Programs

Name        | Description
--------    |:-----------------
`bulk`      | Metropolis Monte Carlo simulation of multicomponent electrolyte solutions. Ions are modelled as charged, hard particles using the primitive model of electrolytes. Uses [scaled Widom analysis](https://doi.org/10.1080/00268978800100203) for calculating single ion activity coefficients. Both the original Fortran 77 code is available (`make bulk_f77`) as well as a modernized version (`make bulk`). They give identical results. The code has been used to study [excess chemical potentials of seasalt](https://doi.org/10.1016/S0304-4203(02)00039-7).

## Building

The build uses `gfortran`, but it should be easy to modify to other compilers.

```sh
make clean
make
```

## Contributors

The codes were usually shared from a common source and modified to suit special
purposes. While there's no single author, here's a list of likely contributors,
mainly from / visitors at Division of Theoretical Chemistry, Lund University, Sweden:

_Bo Jönsson, Cliff Woodward, Bo Svensson, Peter Bolhuis, Jan Forsman,
Torbjörn Åkesson, Magnus Ullner, Mikael Lund._ 

