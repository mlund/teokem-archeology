FROM gcc:15

# Install Python and system dependencies
RUN apt-get update && apt-get install -y \
    python3 python3-pip python3-venv \
    make gdb valgrind \
    vim git curl \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages for Fortran development
RUN pip3 install --break-system-packages \
    fprettify \
    matplotlib \
    numpy

# Verify installations
RUN gfortran --version && \
    python3 --version && \
    fprettify --version

# Set working directory
WORKDIR /workspace

# Set OpenMP threads default
ENV OMP_NUM_THREADS=4
