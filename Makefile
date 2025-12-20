# Simple Makefile for Fortran programs using gfortran
# This Makefile can compile Fortran programs (.f, .f90, .F, .F90 files)

# Compiler
FC = gfortran

# Compiler flags
FFLAGS = -O2 -Wall -Wextra
DEBUGFLAGS = -g -fbacktrace -fcheck=all -Wall -Wextra

# Directories
SRCDIR = .
OBJDIR = obj
BINDIR = bin

# Default executable name (can be overridden: make TARGET=myprogram)
TARGET = program

# Source files to compile (can be overridden: make SOURCES="file1.f90 file2.f90")
# If not specified, tries to find a single program file matching TARGET
ifndef SOURCES
  # Look for a source file with the same name as TARGET
  ifneq ($(wildcard $(TARGET).f90),)
    SOURCES = $(TARGET).f90
  else ifneq ($(wildcard $(TARGET).f),)
    SOURCES = $(TARGET).f
  else ifneq ($(wildcard $(TARGET).F90),)
    SOURCES = $(TARGET).F90
  else ifneq ($(wildcard $(TARGET).F),)
    SOURCES = $(TARGET).F
  else
    # If no matching file found, compile all Fortran files in directory
    SOURCES = $(wildcard $(SRCDIR)/*.f90 $(SRCDIR)/*.f $(SRCDIR)/*.F90 $(SRCDIR)/*.F)
  endif
endif

# Generate object file names from sources
OBJECTS = $(patsubst %,$(OBJDIR)/%.o,$(basename $(notdir $(SOURCES))))

# Default rule
all: directories $(BINDIR)/$(TARGET)

# Create necessary directories
directories:
	@mkdir -p $(OBJDIR) $(BINDIR)

# Link object files to create executable
$(BINDIR)/$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $(OBJECTS)
	@echo "Built $(TARGET) successfully!"

# Compile .f90 files
$(OBJDIR)/%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@ -J$(OBJDIR)

# Compile .f files
$(OBJDIR)/%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@ -J$(OBJDIR)

# Compile .F90 files (with preprocessing)
$(OBJDIR)/%.o: %.F90
	$(FC) $(FFLAGS) -c $< -o $@ -J$(OBJDIR)

# Compile .F files (with preprocessing)
$(OBJDIR)/%.o: %.F
	$(FC) $(FFLAGS) -c $< -o $@ -J$(OBJDIR)

# Debug build
debug: FFLAGS = $(DEBUGFLAGS)
debug: clean all

# Clean build artifacts
clean:
	rm -rf $(OBJDIR) $(BINDIR)
	@echo "Cleaned build artifacts"

# Clean and rebuild
rebuild: clean all

# Help target
help:
	@echo "Makefile for Fortran programs using gfortran"
	@echo ""
	@echo "Usage:"
	@echo "  make                      - Build the program (default target)"
	@echo "  make all                  - Same as 'make'"
	@echo "  make debug                - Build with debug flags"
	@echo "  make clean                - Remove build artifacts"
	@echo "  make rebuild              - Clean and rebuild"
	@echo "  make help                 - Show this help message"
	@echo ""
	@echo "Variables:"
	@echo "  TARGET=name               - Set executable name (default: program)"
	@echo "  SOURCES='file1.f90 ...'   - Specify source files to compile"
	@echo "  FFLAGS='flags'            - Set custom compiler flags"
	@echo ""
	@echo "Examples:"
	@echo "  make TARGET=myprogram"
	@echo "  make TARGET=test_modules SOURCES='mathops.f90 test_modules.f90'"
	@echo "  make debug TARGET=myprogram"

.PHONY: all directories debug clean rebuild help
