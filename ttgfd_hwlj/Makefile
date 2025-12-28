# Makefile for TTGFD_HWLJ - Polymer DFT Code
# Generalized Flory-Dimer theory with Hard-Wall and Lennard-Jones interactions

# Compiler
FC = gfortran

# Program name
PROGRAM = ttgfd_hwlj

# Source files
SRC = ttgfd_hwlj_varbl.f90
INC = t2.inc.f90

# Compiler flags
# Note: Large COMMON blocks (~6GB) require special handling on ARM64
# Using -O2 instead of -O3 to avoid linker fixup errors
# Removed -march=native which can cause ADRP range issues with large data
FFLAGS = -O3

# Strong warning flags
WARNFLAGS = -Wall \
            -Wextra \
            -Wpedantic \
            -Wno-implicit-interface \
            -Wunused \
            -Wuninitialized \
            -Wconversion \
            -Wconversion-extra \
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
$(PROGRAM): $(SRC) $(INC)
	$(FC) $(FFLAGS) $(WARNFLAGS) -o $(PROGRAM) $(SRC)

# Debug build with runtime checks
.PHONY: debug
debug: FFLAGS = -O0 -g
debug: $(SRC) $(INC)
	$(FC) $(FFLAGS) $(WARNFLAGS) $(DEBUGFLAGS) -o $(PROGRAM)_debug $(SRC)

# Optimized build without warnings (for production)
.PHONY: fast
fast: $(SRC) $(INC)
	$(FC) -O3 -ffast-math -funroll-loops -o $(PROGRAM) $(SRC)

# Build with profiling support
.PHONY: profile
profile: FFLAGS = -O1 -pg
profile: $(SRC) $(INC)
	$(FC) $(FFLAGS) -o $(PROGRAM)_prof $(SRC)

# Clean build artifacts
.PHONY: clean
clean:
	rm -f $(PROGRAM) $(PROGRAM)_debug $(PROGRAM)_prof
	rm -f *.o *.mod
	rm -f gmon.out

# Clean output files
.PHONY: cleanout
cleanout:
	rm -f fcdfil epfil

# Full clean (build artifacts and output)
.PHONY: distclean
distclean: clean cleanout

# Run the program
.PHONY: run
run: $(PROGRAM)
	./$(PROGRAM)

# Check if input files exist
.PHONY: check
check:
	@echo "Checking for required input files..."
	@test -f input.tsph || echo "WARNING: input.tsph not found"
	@test -f epfil || echo "WARNING: epfil not found"
	@test -f input.tsph && echo "✓ input.tsph found"
	@test -f epfil && echo "✓ epfil found"

# Format Fortran 90 source files using fprettify
.PHONY: format
format:
	fprettify -i 2 --enable-decl --c-relations $(SRC) $(INC)

# Help target
.PHONY: help
help:
	@echo "TTGFD_HWLJ Makefile"
	@echo "=================="
	@echo ""
	@echo "Targets:"
	@echo "  make          - Build with optimizations and warnings (default)"
	@echo "  make debug    - Build with debugging symbols and runtime checks"
	@echo "  make fast     - Build with aggressive optimizations (no warnings)"
	@echo "  make profile  - Build with profiling support (use with gprof)"
	@echo "  make run      - Build and run the program"
	@echo "  make check    - Check for required input files"
	@echo "  make clean    - Remove compiled binaries and objects"
	@echo "  make cleanout - Remove output files"
	@echo "  make distclean- Remove all generated files"
	@echo "  make format   - Format Fortran 90 source files"
	@echo "  make help     - Show this help message"
	@echo ""
	@echo "Compiler: $(FC)"
	@echo "Source:   $(SRC)"
	@echo ""
