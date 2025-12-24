# Makefile for RAMGEO A100 GPU Optimized Version

# Compilers
FC = nvfortran
CC = nvcc

# Compiler flags
FFLAGS = -acc -ta=tesla:cc80 -Minfo=all,accel -fast -Mlarge_arrays
CFLAGS = -O3 -arch=sm_80 --compiler-options -Wall
LDFLAGS = -acc -ta=tesla:cc80 -lcudart -lstdc++

# Source files
F_SOURCES = ramgeo_a100_optimized.f epade.f
CU_SOURCES = cuda_a100_module.cuf
OBJECTS = $(F_SOURCES:.f=.o) $(CU_SOURCES:.cuf=.o)

# Target executable
TARGET = ramgeo_a100

# Default target
all: $(TARGET)

# Link all objects
$(TARGET): $(OBJECTS)
	$(FC) $(LDFLAGS) -o $@ $(OBJECTS)

# Compile Fortran sources
%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

# Compile CUDA Fortran sources
%.o: %.cuf
	$(FC) $(FFLAGS) -c $< -o $@

# Compile CUDA C sources
%.o: %.cu
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up
clean:
	rm -f $(OBJECTS) $(TARGET) *.mod

# Run test
test: $(TARGET)
	@echo "Running test..."
	./$(TARGET) < test_input.in

# Profile with nvprof
profile: $(TARGET)
	nvprof ./$(TARGET) < large_input.in

# Memory check
memcheck: $(TARGET)
	compute-sanitizer ./$(TARGET) < test_input.in

# Performance check
perf: $(TARGET)
	nsys profile ./$(TARGET) < large_input.in
	nsys stats report1.nsys-rep

# Help message
help:
	@echo "Available targets:"
	@echo "  all     - Build the executable (default)"
	@echo "  clean   - Remove object files and executable"
	@echo "  test    - Run with test input"
	@echo "  profile - Profile with nvprof"
	@echo "  memcheck- Check memory with compute-sanitizer"
	@echo "  perf    - Detailed performance analysis with nsys"
	@echo "  help    - Show this help message"

.PHONY: all clean test profile memcheck perf help