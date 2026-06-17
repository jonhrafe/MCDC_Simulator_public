cmake_minimum_required(VERSION 3.16)
project(MCDC_Simulator LANGUAGES CXX VERSION 2.0.3)

# ---------------------------------------------------------------------------
# Mirrors compile.sh:
#   g++ -O3 -std=c++17 -march=native -flto -lpthread -pthread -w -I src
#   binaries: MC-DC_Simulator (main.cpp) and dataSynth (dataSynth.cpp),
#             both over the shared source list in src/.
# Adds compile_commands.json export so the clangd LSP plugin works.
# ---------------------------------------------------------------------------

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default to a Release build (matches the -O3 in compile.sh) if not specified.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Always export compile_commands.json (consumed by clangd / clangd-lsp plugin).
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# --- Options ----------------------------------------------------------------
# -march=native: great for local runs, but breaks portable/cluster binaries.
# Kept ON by default to match compile.sh; turn OFF for portable builds.
option(MCDC_NATIVE_ARCH "Compile with -march=native (non-portable)" ON)
# compile.sh used -w (silence all warnings). For the CPU-correctness cleanup
# we want diagnostics, so warnings are ON by default here. Turn OFF to match
# the legacy silent build.
option(MCDC_ENABLE_WARNINGS "Enable -Wall -Wextra" ON)

# --- Shared core library ----------------------------------------------------
# Glob the sources so newly added files (e.g. the upcoming RNG module) are
# picked up automatically. CONFIGURE_DEPENDS re-globs when the build runs.
file(GLOB MCDC_CORE_SOURCES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.cpp")
# The two files below carry main(); they must not be in the shared library.
list(REMOVE_ITEM MCDC_CORE_SOURCES
    "${CMAKE_SOURCE_DIR}/src/main.cpp"
    "${CMAKE_SOURCE_DIR}/src/dataSynth.cpp")

add_library(mcdc_core STATIC ${MCDC_CORE_SOURCES})

# Eigen is bundled at src/Eigen and included both as "Eigen/..." and <Eigen/...>,
# so src/ must be on the include path (PUBLIC -> propagates to the executables).
target_include_directories(mcdc_core PUBLIC "${CMAKE_SOURCE_DIR}/src")

find_package(Threads REQUIRED)
target_link_libraries(mcdc_core PUBLIC Threads::Threads)

if(MCDC_NATIVE_ARCH)
    target_compile_options(mcdc_core PUBLIC $<$<CONFIG:Release>:-march=native>)
endif()

if(MCDC_ENABLE_WARNINGS)
    target_compile_options(mcdc_core PRIVATE -Wall -Wextra)
endif()

# Link-time optimization (the -flto from compile.sh), only if the toolchain
# supports it. Applied to the library and both executables below.
include(CheckIPOSupported)
check_ipo_supported(RESULT MCDC_IPO_OK OUTPUT MCDC_IPO_MSG)
if(MCDC_IPO_OK)
    set_property(TARGET mcdc_core PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
else()
    message(STATUS "IPO/LTO not supported: ${MCDC_IPO_MSG}")
endif()

# --- Executables ------------------------------------------------------------
add_executable(MC-DC_Simulator "${CMAKE_SOURCE_DIR}/src/main.cpp")
target_link_libraries(MC-DC_Simulator PRIVATE mcdc_core)

add_executable(dataSynth "${CMAKE_SOURCE_DIR}/src/dataSynth.cpp")
target_link_libraries(dataSynth PRIVATE mcdc_core)

if(MCDC_IPO_OK)
    set_property(TARGET MC-DC_Simulator PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
    set_property(TARGET dataSynth        PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
endif()

# Place the binaries at the repo root, where compile.sh put them and where the
# existing debug/run scripts expect to find them.
set_target_properties(MC-DC_Simulator dataSynth PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")

# --- Tests ------------------------------------------------------------------
# Run with:  ctest --test-dir build --output-on-failure
enable_testing()

find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    # RNG reproducibility: same seed -> identical trajectories (P0.1/P0.2).
    # EXPECTED TO FAIL until the RNG unification lands; turns green after.
    add_test(
        NAME rng_reproducibility
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_SOURCE_DIR}/tests/test_rng_reproducibility.py"
                --bin $<TARGET_FILE:MC-DC_Simulator>
                --conf "${CMAKE_SOURCE_DIR}/debug/debug.conf"
                --workdir "${CMAKE_SOURCE_DIR}")
    set_tests_properties(rng_reproducibility PROPERTIES
        LABELS "reproducibility" TIMEOUT 300)
else()
    message(WARNING "Python3 not found: skipping the rng_reproducibility test.")
endif()
