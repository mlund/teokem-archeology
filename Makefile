# Fortran compiler
FC = gfortran

# Compiler flags
FFLAGS = -std=legacy -O0 -g -funroll-loops -finit-local-zero -fno-automatic

# Source files
SRCDIR = src
SOURCES = $(SRCDIR)/ran2.f $(SRCDIR)/bulk.f

# Target executable
TARGET = bulk

# Default target
all: $(TARGET)

# Compile the program
$(TARGET): $(SOURCES)
	$(FC) $(FFLAGS) $(SOURCES) -o $(TARGET) -lm

# Clean up compiled files
clean:
	rm -f $(TARGET)

# Phony targets
.PHONY: all clean
