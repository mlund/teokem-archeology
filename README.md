# Software archeology from Division of Theoretical Chemistry, Lund University 

Collection of old fortran codes for computational chemistry.

## Programs

Name        | Description
--------    |:-----------------
`bulk`      | Metropolis Monte Carlo simulation of multicomponent electrolyte solutions. Ions are modelled as charged, hard particles using the primitive model of electrolytes. Uses [scaled Widom analysis](https://doi.org/10.1080/00268978800100203) for calculating single ion activity coefficients. Calculates radial distribution functions g(r) for all species pairs. Both the original Fortran 77 code is available (`make bulk_f77`) as well as a modernized version (`make bulk`). They give identical results. The code has been used to study [excess chemical potentials of seasalt](https://doi.org/10.1016/S0304-4203(02)00039-7) and the [validity of primitive models](https://doi.org/10.1021/jp808427f).

## Building

The build uses `gfortran`, but it should be easy to modify to other compilers.

```sh
make clean
make
```

## Code Formatting

Fortran 90 source files (.f90) can be automatically formatted using `fprettify`:

```sh
make format
```

This target:
- Formats all `.f90` files using settings from `.fprettify.yaml`
- Uses 4-space indentation
- Adds whitespace around operators
- Enforces lowercase for keywords
- Creates backup files with `.bak` extension before formatting
- Preserves `.f` (Fortran 77) files unchanged

Formatting settings can be customized by editing `.fprettify.yaml`.

To remove backup files:
```sh
make clean-backups
```

**Installation:**
```sh
pip install fprettify
```

## Analysis Tools

### Radial Distribution Function (RDF) Analysis

The `bulk` program saves radial distribution functions g(r) for all species pairs to `rdf.csv`.
To visualize the RDF data, use the provided Python plotting script:

```sh
./src/plot_rdf.py rdf.csv
```

For interactive exploration of RDF data in Jupyter notebooks:

```python
from plot_rdf import plot_interactive
plot_interactive('rdf.csv')
```

**Requirements:**
- Python 3
- pandas
- matplotlib
- numpy
- ipywidgets (for interactive mode only)

## Contributors

The codes were usually shared from a common source and modified to suit special
purposes. While there's no single author, here's a list of likely contributors,
mainly from / visitors at Division of Theoretical Chemistry, Lund University, Sweden:

_Bo Jönsson, Cliff Woodward, Bo Svensson, Peter Bolhuis, Jan Forsman,
Torbjörn Åkesson, Magnus Ullner, Mikael Lund._ 

