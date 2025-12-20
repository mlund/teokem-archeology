# Fortran compiler
FC = gfortran

# Compiler flags for main program
FFLAGS = -std=legacy -O3 -g -funroll-loops -fno-automatic -finit-local-zero
# Flags for ran2 (needs standard variable handling for SAVE to work)
RAN2FLAGS = -std=legacy -O3 -g

# Source files
SRCDIR = src

# Target executable
TARGET = bulk

# Default target
all: $(TARGET)

# Compile the program
$(TARGET): $(SRCDIR)/ran2.f $(SRCDIR)/bulk.f
	$(FC) $(RAN2FLAGS) -c $(SRCDIR)/ran2.f -o ran2.o
	$(FC) $(FFLAGS) -c $(SRCDIR)/bulk.f -o bulk.o
	$(FC) ran2.o bulk.o -o $(TARGET) -lm

# Clean up compiled files
clean:
	rm -f $(TARGET) *.o

# Run bulk test script
test_bulk:
	cd examples/bulk && bash script.sh

# Phony targets
.PHONY: all clean test_bulk
