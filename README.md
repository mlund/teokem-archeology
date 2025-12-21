# Software archeology from Division of Theoretical Chemistry, Lund University 

Collection of old fortran codes for computational chemistry.

## Programs

- `bulk`. Metropolis Monte Carlo simulation of multicomponent electrolyte solutions.
  Ions are modelled as charged, hard particles using the primitive model of electrolytes.
  Scaled Widom analysis for calculating single ion activity coefficients.
  Both the original Fortran 77 code is available (`make bulk_f77`) as well as a modernized
  version (`make bulk`). They give identical results.

## Building

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

