# Fortran compiler
FC = gfortran

# Compiler flags for Fortran 90 free-form
FFLAGS = -O3 -g -ftree-vectorize -march=native -std=gnu -funroll-loops -fno-automatic -finit-local-zero -fimplicit-none -Wall -Wpedantic -Wshadow -Wextra -Wno-compare-reals -Wno-implicit-interface -Wno-implicit-procedure -Wno-unused-parameter -Waliasing -Wampersand -Wconversion -Wsurprising -Wline-truncation -Wintrinsics-std -Wtabs -Wmaybe-uninitialized -Wuninitialized -fopenmp

# Compiler flags for legacy Fortran 77
FFLAGS_F77 = -std=legacy -ftree-vectorize -march=native -O3 -g -funroll-loops -fno-automatic -finit-local-zero -Wall -Wshadow
RAN2FLAGS_F77 = -std=legacy -O3 -g

# Target executables
TARGET = bulk
TARGET_F77 = bulk_f77

# Default target
all: $(TARGET)

# Compile the program (Fortran 90 free-form)
$(TARGET): bulk.f90
	$(FC) $(FFLAGS) -c bulk.f90 -o bulk.o
	$(FC) bulk.o -o $(TARGET) -lm -fopenmp

# Compile the original Fortran 77 version
$(TARGET_F77): ran2.f bulk.f
	$(FC) $(RAN2FLAGS_F77) -c ran2.f -o ran2_f77.o
	$(FC) $(FFLAGS_F77) -c bulk.f -o bulk_f77.o
	$(FC) ran2_f77.o bulk_f77.o -o $(TARGET_F77) -lm

# Clean up compiled files
clean:
	rm -fR $(TARGET) $(TARGET_F77) *.o *.dSYM

# Format Fortran 90 source files using fprettify
# Settings are read from .fprettify.yaml config file
# Formats .f90 and .inc files (legacy .f files are not formatted)
format:
	@echo "Formatting Fortran 90 source files with fprettify..."
	@for file in *.f90 *_f90.inc; do \
		if [ -f "$$file" ]; then \
			echo "Formatting $$file..."; \
			fprettify "$$file"; \
		fi \
	done
	@echo "Formatting complete!"
	@echo "Note: Legacy .f files are not formatted to preserve historical code"

# Phony targets
.PHONY: all clean format
