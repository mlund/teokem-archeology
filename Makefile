# Fortran compiler
FC = gfortran

# Compiler flags for Fortran 90 free-form
FFLAGS = -O3 -g -std=gnu -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -funroll-loops -fno-automatic -finit-local-zero -fimplicit-none -Wall -Wpedantic -Wshadow -Wextra -Wimplicit-interface -Wimplicit-procedure -Waliasing -Wampersand -Wconversion -Wsurprising -Wline-truncation -Wintrinsics-std -Wtabs
RAN2FLAGS = -O3 -g -fimplicit-none -Wshadow

# Compiler flags for legacy Fortran 77
FFLAGS_F77 = -std=legacy -O3 -g -funroll-loops -fno-automatic -finit-local-zero
RAN2FLAGS_F77 = -std=legacy -O3 -g

# Source files
SRCDIR = src

# Target executables
TARGET = bulk
TARGET_F77 = bulk_f77

# Default target
all: $(TARGET)

# Compile the program (Fortran 90 free-form)
$(TARGET): $(SRCDIR)/ran2.f90 $(SRCDIR)/bulk.f90
	$(FC) $(RAN2FLAGS) -c $(SRCDIR)/ran2.f90 -o ran2.o
	$(FC) $(FFLAGS) -c $(SRCDIR)/bulk.f90 -o bulk.o
	$(FC) ran2.o bulk.o -o $(TARGET) -lm

# Compile the original Fortran 77 version
$(TARGET_F77): $(SRCDIR)/ran2.f $(SRCDIR)/bulk.f
	$(FC) $(RAN2FLAGS_F77) -c $(SRCDIR)/ran2.f -o ran2_f77.o
	$(FC) $(FFLAGS_F77) -c $(SRCDIR)/bulk.f -o bulk_f77.o
	$(FC) ran2_f77.o bulk_f77.o -o $(TARGET_F77) -lm

# Clean up compiled files
clean:
	rm -f $(TARGET) $(TARGET_F77) *.o

# Run bulk test script
test_bulk: $(TARGET)
	cd examples/bulk && bash script.sh

# Phony targets
.PHONY: all clean test_bulk bulk_f77
