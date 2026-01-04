# Makefile for TTGFD_HWLJ - Polymer DFT Code
# Generalized Flory-Dimer theory with Hard-Wall and Lennard-Jones interactions

# Compiler
FC = gfortran

# Program name
PROGRAM = polymer_dft

# Source files
MODULE_SRC = polymer_dft_data.f90
SRC = polymer_dft.f90

# Compiler flags
# Note: Large COMMON blocks (~6GB) require special handling on ARM64
# Successfully using -O3 and -march=native after optimization work
# -fopenmp: Enable OpenMP parallelization and SIMD directives
# Removed -fno-automatic and -finit-local-zero: not needed (verified by testing)
FFLAGS = -O3 -ftree-vectorize -march=native -std=gnu -fopenmp -ffast-math -flto -funroll-loops

# Strong warning flags
WARNFLAGS = -Wall \
            -Wextra \
            -Wpedantic \
            -Wimplicit-interface \
            -Wunused \
            -Wmaybe-uninitialized \
            -Wuninitialized \
            -Wconversion \
            -Wconversion-extra \
            -Wshadow \
            -Wsurprising \
            -Waliasing \
            -Wampersand \
            -Warray-bounds \
            -Wcharacter-truncation \
            -Wline-truncation \
            -Wintrinsics-std \
            -Wno-tabs \
            -Wreal-q-constant \
            -Wunderflow \
            -Wunused-parameter \
            -Wrealloc-lhs \
            -Wrealloc-lhs-all \
            -Wtarget-lifetime \
            -Wcompare-reals

# Runtime checking flags (for debugging)
DEBUGFLAGS = -g \
             -fcheck=all \
             -fbacktrace \
             -ffpe-trap=invalid,zero,overflow \
             -finit-real=snan \
             -finit-integer=-2147483648 \
             -finit-logical=true

# Default target
.PHONY: all
all: $(PROGRAM)

# Build the program with warnings
# Module must be compiled first to generate .mod file
$(PROGRAM): $(MODULE_SRC) $(SRC)
	$(FC) $(FFLAGS) $(WARNFLAGS) -o $(PROGRAM) $(MODULE_SRC) $(SRC)

# Debug build with runtime checks
.PHONY: debug
debug: FFLAGS = -O0 -g
debug: $(MODULE_SRC) $(SRC)
	$(FC) $(FFLAGS) $(WARNFLAGS) $(DEBUGFLAGS) -o $(PROGRAM)_debug $(MODULE_SRC) $(SRC)

# Build original Fortran 77 version with legacy standard (no optimization)
.PHONY: legacy
legacy: legacy.f legacy.inc
	$(FC) $(WARNFLAGS) -Wno-implicit-interface -Wno-unused-parameter -O3 -ftree-vectorize -march=native -std=legacy -fno-automatic -ffast-math -finit-local-zero -o legacy legacy.f -lm

# Clean build artifacts
.PHONY: clean
clean:
	rm -f $(PROGRAM) $(PROGRAM)_debug
	rm -f legacy
	rm -f *.o *.mod *.dat
	rm -f fort.*
	rm -f gmon.out cpu.prof profile.pdf
	rm -f fcdfil epfil
	rm -fR *.dSYM

.PHONY: run
run: $(PROGRAM)
	python3 polymer_dft.py --bdm 0.1 --epslj 0.5 --dz 1.0 --nmon 151 --drho 0.3 --executable ./$(PROGRAM)

.PHONY: run_debug
run_debug: debug
	python3 polymer_dft.py --bdm 0.1 --epslj 0.5 --dz 1.0 --nmon 151 --drho 0.3 --executable ./$(PROGRAM)_debug

.PHONY: run_legacy
run_legacy: legacy
	python3 polymer_dft.py --bdm 0.1 --epslj 0.5 --dz 1.0 --nmon 151 --drho 0.3 --executable ./legacy

# Format Fortran 90 source files using fprettify
.PHONY: format
format:
	fprettify -i 2 --enable-decl --c-relations $(MODULE_SRC) $(SRC)

# Run regression tests
.PHONY: test
test: $(PROGRAM)
	@echo "Running regression tests..."
	@python3 test_regression.py

