# Makefile for TTGFD_HWLJ - Polymer DFT Code
# Generalized Flory-Dimer theory with Hard-Wall and Lennard-Jones interactions

# Compiler
FC = gfortran

# Program name
PROGRAM = ttgfd_hwlj

# Source files
MODULE_SRC = polymer_dft_data.f90
SRC = ttgfd_hwlj_varbl.f90
INC = t2.inc.f90  # Kept for legacy F77 build

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
            -Wno-maybe-uninitialized \
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
            -Wno-unused-parameter \
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
debug: $(SRC) $(INC)
	$(FC) $(FFLAGS) $(WARNFLAGS) $(DEBUGFLAGS) -o $(PROGRAM)_debug $(SRC)

# Build with profiling support (gprof)
.PHONY: profile
profile: FFLAGS = -O1 -pg
profile: $(SRC) $(INC)
	$(FC) $(FFLAGS) -o $(PROGRAM)_prof $(SRC)

# Build with gperftools profiling support (macOS compatible)
.PHONY: profile-gperftools
profile-gperftools: FFLAGS = -O3 -ftree-vectorize -march=native -std=gnu -fno-automatic -finit-local-zero -fopenmp -ffast-math -flto -funroll-loops
profile-gperftools: $(SRC) $(INC)
	$(FC) $(FFLAGS) $(WARNFLAGS) -o $(PROGRAM)_gperf $(SRC)

# Run with gperftools CPU profiler
.PHONY: run-profile
run-profile: profile-gperftools
	@echo "Running with gperftools profiler..."
	@echo "Note: Install gperftools first with: brew install gperftools"
	CPUPROFILE=cpu.prof DYLD_INSERT_LIBRARIES=/opt/homebrew/lib/libprofiler.dylib ./$(PROGRAM)_gperf < input.tsph
	@echo "Profile data saved to cpu.prof"
	@echo "Analyze with: pprof --text ./$(PROGRAM)_gperf cpu.prof"
	@echo "Or generate PDF: pprof --pdf ./$(PROGRAM)_gperf cpu.prof > profile.pdf"

# Build original Fortran 77 version with legacy standard (no optimization)
.PHONY: legacy
legacy: ttgfd_hwlj_varbl.f t2.inc
	$(FC) $(WARNFLAGS) -O3 -ftree-vectorize -march=native -std=legacy -fno-automatic -ffast-math -finit-local-zero -o $(PROGRAM)_f77 ttgfd_hwlj_varbl.f -lm

# Clean build artifacts
.PHONY: clean
clean:
	rm -f $(PROGRAM) $(PROGRAM)_debug $(PROGRAM)_prof $(PROGRAM)_f77 $(PROGRAM)_gperf
	rm -f *.o *.mod
	rm -f gmon.out cpu.prof profile.pdf
	rm -f fcdfil epfil
	rm -fR *.dSYM

.PHONY: run
run: $(PROGRAM)
	python3 run.py --bdm 0.1 --epslj 0.5 --dz 1.0 --nmon 151 --drho 0.3 --executable ./$(PROGRAM)

.PHONY: run_legacy
run_legacy: legacy
	python3 run.py --bdm 0.1 --epslj 0.5 --dz 1.0 --nmon 151 --drho 0.3 --executable ./$(PROGRAM)_f77

# Format Fortran 90 source files using fprettify
.PHONY: format
format:
	fprettify -i 2 --enable-decl --c-relations $(SRC) $(INC)

# Run regression tests
.PHONY: test
test: $(PROGRAM)
	@echo "Running regression tests..."
	@python3 test_regression.py

