# Compiler
CXX = g++

# Common flags (for release build)
CXXFLAGS = -O3 -std=c++17 -march=native -flto -lpthread -pthread  -I.

# Debug flags (additional to common settings, reduced optimizations, adds debug symbols)
DEBUG_FLAGS = -O0 -g -DDEBUG

# Source files
SOURCES = main.cpp \
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
          triangle.cpp \q
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
		  benchmark.cpp

# Object files for release build
OBJECTS = $(SOURCES:.cpp=.o)

# Object files for debug build (separate set)
OBJECTS_DEBUG = $(SOURCES:.cpp=.debug.o)

# Headers (if you have headers that affect all sources, list them here)
HEADERS = constants.h

# Final executables
TARGET = ../MC-DC_Simulator
DEBUG_TARGET = ../MC-DC_Simulator_debug

# Default rule: build the release version
all: $(TARGET)

# Release link step
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJECTS)

# Release compilation step
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Debug target
debug: $(DEBUG_TARGET)

# Debug link step (using debug objects)
$(DEBUG_TARGET): $(OBJECTS_DEBUG)
	$(CXX) $(CXXFLAGS) $(DEBUG_FLAGS) -o $@ $(OBJECTS_DEBUG)

# Debug compilation step (note: uses .debug.o)
%.debug.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(DEBUG_FLAGS) -c $< -o $@

# Clean up intermediate and final files
clean:
	rm -f $(OBJECTS) $(OBJECTS_DEBUG) $(TARGET) $(DEBUG_TARGET)
