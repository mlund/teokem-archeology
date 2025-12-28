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
FFLAGS = -O3 -ftree-vectorize -march=native -std=gnu 

# Strong warning flags
WARNFLAGS = -Wall \
            -Wextra \
            -Wpedantic \
            -Wno-implicit-interface \
            -Wunused \
            -Wno-maybe-uninitialized \
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

# Build with profiling support
.PHONY: profile
profile: FFLAGS = -O1 -pg
profile: $(SRC) $(INC)
	$(FC) $(FFLAGS) -o $(PROGRAM)_prof $(SRC)

# Build original Fortran 77 version with legacy standard
.PHONY: legacy
legacy: ttgfd_hwlj_varbl.f t2.inc
	$(FC) -O3 -std=legacy -o $(PROGRAM)_f77 ttgfd_hwlj_varbl.f -lm

# Clean build artifacts
.PHONY: clean
clean:
	rm -f $(PROGRAM) $(PROGRAM)_debug $(PROGRAM)_prof $(PROGRAM)_f77
	rm -f *.o *.mod
	rm -f gmon.out
	rm -f fcdfil epfil
	rm -fR *.dsym

# Run the program
.PHONY: run
run: $(PROGRAM)
	./$(PROGRAM)

# Format Fortran 90 source files using fprettify
.PHONY: format
format:
	fprettify -i 2 --enable-decl --c-relations $(SRC) $(INC)

