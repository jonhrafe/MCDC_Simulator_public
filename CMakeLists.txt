CMAKE_MINIMUM_REQUIRED(VERSION 3.10)

project(MC-DC_Simulator)
set(CMAKE_CXX_STANDARD 17)

# Find all the source files for both projects
file(GLOB SOURCES "src/*.cpp")
list(REMOVE_ITEM SOURCES "${CMAKE_SOURCE_DIR}/src/dataSynth.cpp" "${CMAKE_SOURCE_DIR}/src/main.cpp")

#message(STATUS ${SOURCES})

# Add the executable target
add_executable(MC-DC_Simulator ${CMAKE_SOURCE_DIR}/src/main.cpp ${SOURCES})
add_executable(DataSynth ${CMAKE_SOURCE_DIR}/src/dataSynth.cpp ${SOURCES})
target_include_directories(MC-DC_Simulator PRIVATE src)
target_include_directories(DataSynth PRIVATE src)