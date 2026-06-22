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

# --- Optional CUDA (GPU path) ----------------------------------------------
# Detected, not required: CPU-only machines still build the full simulator.
# When present, src/*.cu is compiled into a separate mcdc_gpu library and the
# MCDC_HAVE_CUDA macro enables the `gpu` config flag.
include(CheckLanguage)
check_language(CUDA)
set(MCDC_HAVE_CUDA OFF)
if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    set(MCDC_HAVE_CUDA ON)
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES 86)   # RTX 3080 (Ampere); override with -DCMAKE_CUDA_ARCHITECTURES=...
    endif()
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    message(STATUS "CUDA found (${CMAKE_CUDA_COMPILER}); GPU path ENABLED, arch ${CMAKE_CUDA_ARCHITECTURES}.")
else()
    message(STATUS "CUDA not found; GPU path DISABLED (CPU-only build).")
endif()

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

# --- GPU library (CUDA .cu sources) -----------------------------------------
# Kept separate from mcdc_core so host LTO and nvcc don't interfere. mcdc_core
# links it and gains MCDC_HAVE_CUDA so the `gpu` path compiles in.
if(MCDC_HAVE_CUDA)
    file(GLOB MCDC_CUDA_SOURCES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.cu")
    add_library(mcdc_gpu STATIC ${MCDC_CUDA_SOURCES})
    target_include_directories(mcdc_gpu PUBLIC "${CMAKE_SOURCE_DIR}/src")
    set_target_properties(mcdc_gpu PROPERTIES
        CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES}
        CUDA_STANDARD 17
        POSITION_INDEPENDENT_CODE ON)
    find_package(CUDAToolkit REQUIRED)
    target_link_libraries(mcdc_gpu PUBLIC CUDA::cudart)
    target_link_libraries(mcdc_core PUBLIC mcdc_gpu)
    target_compile_definitions(mcdc_core PUBLIC MCDC_HAVE_CUDA)
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
    # Green since the RNG unification (P0.1) landed.
    add_test(
        NAME rng_reproducibility
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_SOURCE_DIR}/tests/test_rng_reproducibility.py"
                --bin $<TARGET_FILE:MC-DC_Simulator>
                --conf "${CMAKE_SOURCE_DIR}/tests/accuracy/repro.conf"
                --workdir "${CMAKE_SOURCE_DIR}")
    set_tests_properties(rng_reproducibility PROPERTIES
        LABELS "reproducibility" TIMEOUT 300)

    # Accuracy golden-master: a fixed impermeable sphere-mesh experiment whose
    # diffusion signal (DWI real/intra/extra) must stay within tolerance of the
    # committed baseline. Guards core physics against regressions from refactors,
    # parallelization and the upcoming permeability work; also the reference a
    # GPU port is checked against (loosen --rtol there). Regenerate the baseline
    # only for an intentional, reviewed change:
    #   python3 tests/test_accuracy_golden.py --update
    add_test(
        NAME accuracy_golden
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_SOURCE_DIR}/tests/test_accuracy_golden.py"
                --bin $<TARGET_FILE:MC-DC_Simulator>
                --conf "${CMAKE_SOURCE_DIR}/tests/accuracy/sphere_impermeable.conf"
                --golden-dir "${CMAKE_SOURCE_DIR}/tests/accuracy/golden"
                --workdir "${CMAKE_SOURCE_DIR}")
    set_tests_properties(accuracy_golden PROPERTIES
        LABELS "accuracy" TIMEOUT 600)

    # Second golden: two impermeable meshes (original sphere + a displaced
    # deformed sphere) loaded via a ply_file_list_scale_permeability list, walkers
    # initialised intra, finite T2. Exercises multi-mesh restricted diffusion; the
    # extra signal must stay ~0 (impermeable leak detector). Same harness/baseline
    # workflow as accuracy_golden.
    add_test(
        NAME accuracy_golden_two_meshes
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_SOURCE_DIR}/tests/test_accuracy_golden.py"
                --bin $<TARGET_FILE:MC-DC_Simulator>
                --conf "${CMAKE_SOURCE_DIR}/tests/accuracy/two_meshes_impermeable.conf"
                --golden-dir "${CMAKE_SOURCE_DIR}/tests/accuracy/golden"
                --workdir "${CMAKE_SOURCE_DIR}")
    set_tests_properties(accuracy_golden_two_meshes PROPERTIES
        LABELS "accuracy" TIMEOUT 600)

    # Native analytic obstacles (NOT meshes): impermeable hexagonal packings of
    # the simulator's own Sphere / Cylinder objects, walkers init intra. These
    # exercise Sphere::checkCollision / Cylinder::checkCollision -- a code path
    # entirely separate from the PLY mesh goldens above.
    add_test(
        NAME accuracy_golden_sphere_hex
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_SOURCE_DIR}/tests/test_accuracy_golden.py"
                --bin $<TARGET_FILE:MC-DC_Simulator>
                --conf "${CMAKE_SOURCE_DIR}/tests/accuracy/sphere_hex_impermeable.conf"
                --golden-dir "${CMAKE_SOURCE_DIR}/tests/accuracy/golden"
                --workdir "${CMAKE_SOURCE_DIR}")
    set_tests_properties(accuracy_golden_sphere_hex PROPERTIES
        LABELS "accuracy" TIMEOUT 600)

    add_test(
        NAME accuracy_golden_cylinder_hex
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_SOURCE_DIR}/tests/test_accuracy_golden.py"
                --bin $<TARGET_FILE:MC-DC_Simulator>
                --conf "${CMAKE_SOURCE_DIR}/tests/accuracy/cylinder_hex_impermeable.conf"
                --golden-dir "${CMAKE_SOURCE_DIR}/tests/accuracy/golden"
                --workdir "${CMAKE_SOURCE_DIR}")
    set_tests_properties(accuracy_golden_cylinder_hex PROPERTIES
        LABELS "accuracy" TIMEOUT 600)

    # Powles permeability model: a single permeable sphere with distinct Di/De so
    # the prescribed prob_i_e != prob_e_i (directional). Asserts the empirical
    # per-encounter crossing probability (crossings/hits) matches the prescribed
    # value within tolerance -- i.e. the Powles probabilities are applied
    # correctly. num_process 1 (exact counters).
    add_test(
        NAME accuracy_permeability
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_SOURCE_DIR}/tests/test_permeability_powles.py"
                --bin $<TARGET_FILE:MC-DC_Simulator>
                --conf "${CMAKE_SOURCE_DIR}/tests/accuracy/permeable_sphere.conf"
                --workdir "${CMAKE_SOURCE_DIR}")
    set_tests_properties(accuracy_permeability PROPERTIES
        LABELS "permeability" TIMEOUT 600)
    # Per-obstacle property tracking: distinct Di/T2/kappa per mesh, permeable, with
    # `debug` on. Asserts each compartment maps to a single correct (Di,T2) across
    # all steps (no cross-contamination), crossings occur, and per-obstacle kappa is
    # correct. Guards the per-obstacle property machinery.
    add_test(
        NAME obstacle_properties
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/test_obstacle_properties.py")
    set_tests_properties(obstacle_properties PROPERTIES LABELS "permeability" TIMEOUT 600)

    # Negative tests: malformed configs must be REJECTED (clean non-zero exit +
    # error), valid/benign ones accepted. Guards the input validators and the P0.4
    # enforcement fix (validation must fail in release builds, not just debug).
    add_test(
        NAME format_checks
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/test_format_checks.py")
    set_tests_properties(format_checks PROPERTIES LABELS "validation" TIMEOUT 300)

    # GPU free-diffusion (Phase-1 PoC): the on-device dynamics+synthesis must
    # reproduce S/S0 = exp(-b*D). Self-skips (pass) if the binary was built
    # without CUDA / no GPU, so CPU-only machines still pass the suite.
    add_test(
        NAME gpu_free_diffusion
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/test_gpu_free_diffusion.py")
    set_tests_properties(gpu_free_diffusion PROPERTIES LABELS "gpu" TIMEOUT 300)

    # GPU Phase 2a/2b: restricted diffusion in a single analytic obstacle; GPU vs
    # CPU signal must agree within MC noise (guards the device collision/reflection
    # kernels). Self-skips on non-CUDA builds.
    add_test(
        NAME gpu_restricted_sphere
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/test_gpu_restricted.py"
                "${CMAKE_SOURCE_DIR}/tests/accuracy/sphere_restricted.conf")
    set_tests_properties(gpu_restricted_sphere PROPERTIES LABELS "gpu" TIMEOUT 300)

    add_test(
        NAME gpu_restricted_cylinder
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/test_gpu_restricted.py"
                "${CMAKE_SOURCE_DIR}/tests/accuracy/cylinder_restricted.conf")
    set_tests_properties(gpu_restricted_cylinder PROPERTIES LABELS "gpu" TIMEOUT 300)

    # GPU Phase 2c: multiple analytic obstacles (per-step nearest-collision over an
    # array). Sphere hex packing, and several cylinders with a large axial voxel so
    # the (un-implemented) periodic BC is negligible. GPU vs CPU within MC noise.
    add_test(
        NAME gpu_hex_spheres
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/test_gpu_restricted.py"
                "${CMAKE_SOURCE_DIR}/tests/accuracy/sphere_hex_gpu.conf")
    set_tests_properties(gpu_hex_spheres PROPERTIES LABELS "gpu" TIMEOUT 300)

    add_test(
        NAME gpu_multi_cylinder
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/test_gpu_restricted.py"
                "${CMAKE_SOURCE_DIR}/tests/accuracy/multi_cylinder.conf")
    set_tests_properties(gpu_multi_cylinder PROPERTIES LABELS "gpu" TIMEOUT 300)

    # GPU Phase 2d: restricted diffusion inside a triangle MESH (unitMesh sphere,
    # brute-force Moller-Trumbore over all triangles). GPU vs CPU within noise.
    add_test(
        NAME gpu_mesh_sphere
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/test_gpu_restricted.py"
                "${CMAKE_SOURCE_DIR}/tests/accuracy/mesh_restricted.conf")
    set_tests_properties(gpu_mesh_sphere PROPERTIES LABELS "gpu" TIMEOUT 300)

    # GPU Phase 2e: voxel PERIODIC boundary conditions (free diffusion, tight box).
    # The walker wraps (pos_v) while the real position (pos_r) stays unwrapped; GPU
    # must reproduce the CPU periodic-medium signal. (Mirror BC is not yet matched
    # to the CPU's unfolding convention.)
    add_test(
        NAME gpu_periodic_box
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/test_gpu_restricted.py"
                "${CMAKE_SOURCE_DIR}/tests/accuracy/free_periodic.conf")
    set_tests_properties(gpu_periodic_box PROPERTIES LABELS "gpu" TIMEOUT 300)
else()
    message(WARNING "Python3 not found: skipping the rng_reproducibility and accuracy_golden tests.")
endif()
