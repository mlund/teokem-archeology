# Dev Container Setup

This directory contains the Docker-based development container configuration for the Polymer DFT project.

## What's Included

- **GCC/GFortran 15.2.0** - Latest compiler with full Fortran 2008/2018 support
- **Python 3** - For running polymer_dft.py driver script
- **fprettify** - Fortran code formatter
- **matplotlib** - Python plotting library
- **numpy** - Numerical computing
- **Development tools** - make, cmake, gdb, valgrind, git

## Quick Start

### Using VS Code

1. Install the "Dev Containers" extension in VS Code
2. Open this project folder in VS Code
3. Press `F1` and select "Dev Containers: Reopen in Container"
4. Wait for the container to build (first time only, ~2-5 minutes)
5. Once built, you'll have a complete development environment

### Using Command Line

Build the container:
```bash
docker build -t polymer-dft-dev .devcontainer
```

Run the container:
```bash
docker run -it --rm -v $(pwd):/workspace polymer-dft-dev
```

## Verification

After the container starts, verify the installation:

```bash
# Check compiler version
gfortran --version
# Should show: GNU Fortran (GCC) 15.x.x

# Check Python tools
python3 --version
fprettify --version
python3 -c "import matplotlib; import numpy; print('OK')"

# Build and test the project
make clean && make
make test
```

## Using the Environment

### Compile the project
```bash
make clean
make          # Build optimized version
make debug    # Build debug version
make legacy   # Build F77 version
```

### Run simulations
```bash
make run        # Run optimized version
make run_debug  # Run debug version
make run_legacy # Run legacy version
```

### Format code
```bash
make format
# or
fprettify -i 2 --enable-decl --c-relations polymer_dft.f90
```

### Debug with GDB
```bash
make debug
gdb ./polymer_dft_debug
```

## Configuration

### OpenMP Threads

The default is set to 4 threads via `OMP_NUM_THREADS=4` environment variable. You can override this:

```bash
export OMP_NUM_THREADS=8
make run
```

### VS Code Extensions

The following extensions are automatically installed:
- `fortran-lang.linter-gfortran` - Fortran linting
- `ekibun.fortranbreaker` - Fortran debugging support
- `ms-python.python` - Python support
- `ms-vscode.cpptools` - C/C++ tools (for gdb integration)

## Troubleshooting

### Container fails to build

Check your Docker installation:
```bash
docker --version
docker info
```

### Slow build times

The first build downloads ~500MB of base images. Subsequent rebuilds are much faster due to layer caching.

### Permission issues

The container runs as root by default. If you need to change ownership of generated files:
```bash
sudo chown -R $USER:$USER .
```

## Updating

To rebuild with the latest GCC 15 updates:

```bash
# In VS Code
F1 -> "Dev Containers: Rebuild Container"

# Or from command line
docker build --no-cache -t polymer-dft-dev .devcontainer
```

## Notes

- The container uses Debian Trixie as the base (via gcc:15)
- All project files are mounted at `/workspace`
- Changes to files are immediately reflected on both host and container
- The container has full ptrace capabilities for debugging
