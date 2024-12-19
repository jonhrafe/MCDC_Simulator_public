#!/bin/bash -l

# Move into the source directory
cd src

# Common source files
COMMON_FILES="simulablesequence.cpp vertex.cpp obstacle.cpp collision.cpp scheme.cpp voxel.cpp cylinder.cpp walker.cpp \
mcsimulation.cpp parallelmcsimulation.cpp trajectory.cpp triangle.cpp parameters.cpp plyobstacle.cpp pgsesequence.cpp \
dynamicsSimulation.cpp simerrno.cpp collisionsphere.cpp cylindergammadistribution.cpp sentinel.cpp subdivision.cpp \
gradientwaveform.cpp propagator.cpp sphere.cpp sphere.h spheregammadistribution.cpp benchmark.cpp"

# Build MC-DC_Simulator
echo "Compiling MC-DC_Simulator..."
g++ -O3 -std=c++17 -march=native -flto -lpthread -pthread -w -I. main.cpp $COMMON_FILES -o ../MC-DC_Simulator
if [ $? -ne 0 ]; then
    echo "Error compiling MC-DC_Simulator"
    exit 1
fi

# Build dataSynth
echo "Compiling dataSynth..."
g++ -O3 -std=c++17 -march=native -flto -lpthread -pthread -w -I. dataSynth.cpp $COMMON_FILES -o ../dataSynth
if [ $? -ne 0 ]; then
    echo "Error compiling dataSynth"
    exit 1
fi

echo "Compilation finished successfully."