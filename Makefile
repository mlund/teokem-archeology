# Root Makefile - delegates to subdirectories

# Default target
all:
	$(MAKE) -C bulk

# Build ttgfd_hwlj polymer DFT code
ttgfd_hwlj:
	$(MAKE) -C ttgfd_hwlj

# Build Fortran 77 version
bulk_f77:
	$(MAKE) -C bulk bulk_f77

# Clean up compiled files
clean:
	$(MAKE) -C bulk clean
	$(MAKE) -C ttgfd_hwlj clean
	rm -fR *.o *.dSYM

# Format Fortran 90 source files using fprettify
format:
	$(MAKE) -C bulk format
	$(MAKE) -C ttgfd_hwlj format

# Phony targets
.PHONY: all clean bulk_f77 format ttgfd_hwlj
