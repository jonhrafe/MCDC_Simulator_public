# Makefile

# Compiler Definitions
CXX = g++
NVCC = nvcc

# Host Include Paths
HOST_INCLUDES = -I. -IEigen

# Device Include Paths
DEVICE_INCLUDES = -I.

# Common C++ Flags for CPU Build
CXXFLAGS_RELEASE = -std=c++17 -pthread -w -O3 $(HOST_INCLUDES)

# CUDA Compiler Flags for CUDA Build (Release)
NVCCFLAGS_RELEASE = -std=c++11 -gencode arch=compute_61,code=sm_61 -lcudart -lstdc++ -lm -Xcompiler -w \
                    --compiler-bindir=/usr/bin/gcc-10

# Debug Flags
CXXFLAGS_DEBUG = -std=c++17 -pthread -g -O0 $(HOST_INCLUDES)
NVCCFLAGS_DEBUG = -std=c++11 -gencode arch=compute_61,code=sm_61 -g -G -lcudart -lstdc++ -lm -Xcompiler -g \
                  --compiler-bindir=/usr/bin/gcc-10

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

# Debug Object Files
DEBUG_OBJECTS_CPP = $(SOURCES_CPP:.cpp=.debug_o)
DEBUG_OBJECTS_CU = $(SRC_CU:.cu=.debug_o)

# Static Library
LIB_HOST = libhost.a

# Executables
TARGET = ../MC-DC_Simulator
CUDA_TARGET = ../CUDAIntersections
PLYInt_TARGET = ../PLYIntersections

DEBUG_TARGET = ../MC-DC_Simulator_debug
CPUInt_DEBUG_TARGET = ../PLYIntersections_debug
GPUInt_DEBUG_TARGET = ../CUDAIntersections_debug

# Default Target: Build CPU Executable (Release)
all: $(TARGET)

# CPU Compilation for Main Simulation (Release)
$(TARGET): $(OBJECTS_CPP)
	$(CXX) $(CXXFLAGS_RELEASE) -o $@ $(OBJECTS_CPP)

# Release Pattern Rules
%.o: %.cpp
	$(CXX) $(CXXFLAGS_RELEASE) -c $< -o $@

%.o: %.cu
	$(NVCC) -c $< -o $@ $(NVCCFLAGS_RELEASE) $(DEVICE_INCLUDES)

# Debug Pattern Rules
%.debug_o: %.cpp
	$(CXX) $(CXXFLAGS_DEBUG) -c $< -o $@

%.debug_o: %.cu
	$(NVCC) -c $< -o $@ $(NVCCFLAGS_DEBUG) $(DEVICE_INCLUDES)

# Debug build for main.cpp
debug: $(DEBUG_OBJECTS_CPP)
	$(CXX) $(CXXFLAGS_DEBUG) -o $(DEBUG_TARGET) $^

# Debug build for CPUPLYInt.cpp
CPUInt_debug: CPUPLYInt.debug_o $(PLY_OBJECTS:.o=.debug_o)
	$(CXX) $(CXXFLAGS_DEBUG) -o $(CPUInt_DEBUG_TARGET) $(PLY_OBJECTS:.o=.debug_o)

# Debug build for TriangleInt.cu
GPUInt_debug: $(DEBUG_OBJECTS_CU)
	$(NVCC) $(NVCCFLAGS_DEBUG) -o $(GPUInt_DEBUG_TARGET) $^ $(LIB_HOST)

# Create Static Library
$(LIB_HOST): $(OBJECTS_CPP)
	ar rcs $@ $(OBJECTS_CPP)

# New Target for PLY Intersections
PLY_OBJECTS = simulablesequence.o \
              vertex.o \
              obstacle.o \
              collision.o \
              scheme.o \
              voxel.o \
              cylinder.o \
              walker.o \
              mcsimulation.o \
              parallelmcsimulation.o \
              trajectory.o \
              triangle.o \
              parameters.o \
              plyobstacle.o \
              pgsesequence.o \
              dynamicsSimulation.o \
              simerrno.o \
              collisionsphere.o \
              cylindergammadistribution.o \
              sentinel.o \
              subdivision.o \
              gradientwaveform.o \
              propagator.o \
              sphere.o \
              spheregammadistribution.o \
              benchmark.o \
              FixedGrid.o \
              CPUPLYInt.o

plyInt: $(PLYInt_TARGET)

$(PLYInt_TARGET): $(PLY_OBJECTS)
	$(CXX) $(CXXFLAGS_RELEASE) -o $@ $(PLY_OBJECTS)

# CUDA Target
cudaInt: $(CUDA_TARGET)
# Datasynth Compilation
$(CUDA_TARGET): $(OBJECTS_CU) $(LIB_HOST)
	$(NVCC) $(NVCCFLAGS) -o $@ $(OBJECTS_CU) -L. -lhost -w
    
# Clean Target
clean:
	rm -f $(OBJECTS_CPP) $(OBJECTS_CU) $(DEBUG_OBJECTS_CPP) $(DEBUG_OBJECTS_CU) $(LIB_HOST) $(TARGET) $(CUDA_TARGET) $(PLYInt_TARGET) $(DEBUG_TARGET) $(CPUInt_DEBUG_TARGET) $(GPUInt_DEBUG_TARGET) *.o *.debug_o
