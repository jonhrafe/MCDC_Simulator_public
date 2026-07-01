cmake_minimum_required(VERSION 3.16)
project(MCDC_Simulator LANGUAGES CXX VERSION 2.2.0)

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

# Compiler-specific flags. The GCC/Clang flags below (-march=native, -Wall) are not valid for MSVC.
if(MSVC)
    # M_PI is non-standard and must be enabled before <cmath>; /permissive- makes MSVC accept the
    # standard alternative operator tokens (and / or / not) used in the sources.
    target_compile_definitions(mcdc_core PUBLIC _USE_MATH_DEFINES)
    target_compile_options(mcdc_core PRIVATE /permissive-)
    # Force-include constants.h into every TU so the fixed-name uint/ulong aliases (which MSVC,
    # unlike Linux/macOS, does not provide via system headers) are always in scope. PUBLIC so the
    # executables get it too; constants.h is on the include path (src/) and is include-guarded.
    target_compile_options(mcdc_core PUBLIC /FIconstants.h)
    if(MCDC_ENABLE_WARNINGS)
        target_compile_options(mcdc_core PRIVATE /W3)
    endif()
else()
    if(MCDC_NATIVE_ARCH)
        target_compile_options(mcdc_core PUBLIC $<$<CONFIG:Release>:-march=native>)
    endif()
    if(MCDC_ENABLE_WARNINGS)
        target_compile_options(mcdc_core PRIVATE -Wall -Wextra)
    endif()
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

# --- Tests (development tree only) ------------------------------------------
# The test suite under tests/ ships with the development sources and is excluded
# from release clones. When it is present it is wired up for CTest:
#   ctest --test-dir build --output-on-failure
if(EXISTS "${CMAKE_SOURCE_DIR}/tests")
enable_testing()

find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    # RNG reproducibility: same seed -> identical trajectories.
    # Green since the RNG unification landed.
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

    # Per-obstacle (list-column) permeability: a spheres_list line "x y z r kappa T2" must
    # drive exchange via its per-sphere kappa column (not only the global `permeability`).
    # Regressed silently when the Sphere/Cylinder copy ctors dropped the Obstacle base on
    # push_back. Self-contained (writes its own configs/lists); asserts per-sphere kappa
    # both produces exchange AND matches the same global kappa.
    add_test(
        NAME permeability_list
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/test_permeability_list.py")
    set_tests_properties(permeability_list PROPERTIES LABELS "permeability" TIMEOUT 300)
    # Per-obstacle property tracking: distinct Di/T2/kappa per mesh, permeable, with
    # `debug` on. Asserts each compartment maps to a single correct (Di,T2) across
    # all steps (no cross-contamination), crossings occur, and per-obstacle kappa is
    # correct. Guards the per-obstacle property machinery.
    add_test(
        NAME obstacle_properties
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/test_obstacle_properties.py")
    set_tests_properties(obstacle_properties PROPERTIES LABELS "permeability" TIMEOUT 600)

    # Negative tests: malformed configs must be REJECTED (clean non-zero exit +
    # error), valid/benign ones accepted. Guards the input validators
    # (validation must fail in release builds, not just debug).
    add_test(
        NAME format_checks
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/test_format_checks.py")
    set_tests_properties(format_checks PROPERTIES LABELS "validation" TIMEOUT 300)

    # General WAVEFORM protocol: a PGSE encoded as an arbitrary gradient waveform (bipolar two-lobe)
    # must reproduce both the native PGSE signal and the analytic exp(-bD) on free diffusion. Guards
    # scheme.type==WAVEFORM -> GradientWaveform synthesis.
    add_test(
        NAME waveform_protocol
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/test_waveform.py")
    set_tests_properties(waveform_protocol PROPERTIES LABELS "waveform" TIMEOUT 300)

    # GAMMA-distributed packings (seeded -> reproducible substrate): golden DWI for a gamma sphere
    # and a gamma cylinder packing. Guards the gamma generator + the SI unit homogenization
    # (beta/min_radius in metres).
    add_test(
        NAME accuracy_golden_gamma_sphere
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/test_accuracy_golden.py"
                --bin $<TARGET_FILE:MC-DC_Simulator>
                --conf "${CMAKE_SOURCE_DIR}/tests/accuracy/gamma_sphere.conf"
                --golden-dir "${CMAKE_SOURCE_DIR}/tests/accuracy/golden"
                --workdir "${CMAKE_SOURCE_DIR}")
    set_tests_properties(accuracy_golden_gamma_sphere PROPERTIES LABELS "accuracy" TIMEOUT 600)

    add_test(
        NAME accuracy_golden_gamma_cylinder
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/test_accuracy_golden.py"
                --bin $<TARGET_FILE:MC-DC_Simulator>
                --conf "${CMAKE_SOURCE_DIR}/tests/accuracy/gamma_cylinder.conf"
                --golden-dir "${CMAKE_SOURCE_DIR}/tests/accuracy/golden"
                --workdir "${CMAKE_SOURCE_DIR}")
    set_tests_properties(accuracy_golden_gamma_cylinder PROPERTIES LABELS "accuracy" TIMEOUT 600)
else()
    message(WARNING "Python3 not found: skipping the rng_reproducibility and accuracy_golden tests.")
endif()
endif()  # tests/ present (development tree)
