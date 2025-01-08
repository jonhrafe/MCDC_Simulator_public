# Makefile

# Compiler Definitions
CXX = g++
NVCC = nvcc

# Host Include Paths
HOST_INCLUDES = -I. -IEigen

# Device Include Paths
DEVICE_INCLUDES = -I.

# Common C++ Flags for CPU Build
CXXFLAGS = -std=c++17 -pthread  -w -O3 $(HOST_INCLUDES)

# CUDA Compiler Flags for CUDA Build
#NVCCFLAGS = -std=c++11 -gencode arch=compute_61,code=sm_61 -lcudart -lstdc++ -lm \
            --compiler-bindir=/usr/bin/gcc-10

# CUDA Compiler Flags for CUDA Build
NVCCFLAGS = -std=c++14 -lcudart -lstdc++ -lm -w \
            --compiler-bindir=/usr/bin/gcc-8
# Source Files
SOURCES_CPP = main.cpp \
              simulablesequence.cpp \
              vertex.cpp \
              obstacle.cpp \
              collision.cpp \
              scheme.cpp \
              voxel.cpp \
              cylinder.cpp \
              walker.cpp \
              mcsimulation.cpp \
              parallelmcsimulation.cpp \
              trajectory.cpp \
              triangle.cpp \
              parameters.cpp \
              plyobstacle.cpp \
              pgsesequence.cpp \
              dynamicsSimulation.cpp \
              simerrno.cpp \
              collisionsphere.cpp \
              cylindergammadistribution.cpp \
              sentinel.cpp \
              subdivision.cpp \
              gradientwaveform.cpp \
              propagator.cpp \
              sphere.cpp \
              spheregammadistribution.cpp \
              benchmark.cpp \
              FixedGrid.cpp

SRC_CU = CUDAInt.cu TriangleInt.cu

# Object Files
OBJECTS_CPP = $(SOURCES_CPP:.cpp=.o)
OBJECTS_CU = $(SRC_CU:.cu=.o)

# Static Library
LIB_HOST = libhost.a

# Executables
TARGET = ../MC-DC_Simulator
CUDA_TARGET = ../CUDAInt

# Default Target: Build CPU Executable
all: $(TARGET)

# CPU Compilation
$(TARGET): $(OBJECTS_CPP)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJECTS_CPP)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Create Static Library
$(LIB_HOST): $(OBJECTS_CPP)
	ar rcs $@ $(OBJECTS_CPP)

# CUDA Target
cuda: $(CUDA_TARGET)

# CUDA Compilation
$(CUDA_TARGET): $(OBJECTS_CU) $(LIB_HOST)
	$(NVCC) $(NVCCFLAGS) -o $@ $(OBJECTS_CU) -L. -lhost

%.o: %.cu
	$(NVCC) -c $< -o $@ $(NVCCFLAGS) $(DEVICE_INCLUDES)

# Clean Target
clean:
	rm -f $(OBJECTS_CPP) $(OBJECTS_CU) $(LIB_HOST) $(TARGET) $(CUDA_TARGET) *.o
