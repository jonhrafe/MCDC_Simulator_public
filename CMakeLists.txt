cmake_minimum_required(VERSION 3.16)
project(MCDC_Simulator LANGUAGES CXX VERSION 2.3.2)

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

# --- Tests -------------------------------------------------------------------
# The test suite under tests/ (Python scripts + tests/accuracy references) is part of the
# repository and is wired up for CTest:
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

    # Unit check: the T2 column of spheres_list / cylinders_list is in SECONDS (SI .conf).
    add_test(
        NAME list_units
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/test_list_units.py")
    set_tests_properties(list_units PROPERTIES LABELS "validation" TIMEOUT 120)

    # Permeable-membrane equilibrium (volume preservation) with D_intra != D_extra: the intra
    # density of a uniformly seeded permeable sphere must stay at its volume fraction, with no
    # walkers discarded as illegal crossings. `_mesh` repeats it on a triangulated PLY sphere and
    # additionally checks that an all-intra start relaxes like the analytic sphere.
    add_test(
        NAME permeability_equilibrium
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_SOURCE_DIR}/tests/test_permeability_equilibrium.py"
                --bin $<TARGET_FILE:MC-DC_Simulator>
                --scheme "${CMAKE_SOURCE_DIR}/tests/accuracy/PGSE_sample_scheme.scheme")
    set_tests_properties(permeability_equilibrium PROPERTIES LABELS "permeability" TIMEOUT 600)
    add_test(
        NAME permeability_equilibrium_mesh
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_SOURCE_DIR}/tests/test_permeability_equilibrium.py"
                --substrate mesh
                --bin $<TARGET_FILE:MC-DC_Simulator>
                --scheme "${CMAKE_SOURCE_DIR}/tests/accuracy/PGSE_sample_scheme.scheme")
    set_tests_properties(permeability_equilibrium_mesh PROPERTIES LABELS "permeability" TIMEOUT 1500)

    # The same equilibrium with steps comparable to the obstacle (r = 1 um, steps 0.77 / 0.39 um):
    # a walker that crossed the membrane can meet it again in the same step and must again cross
    # or reflect. Sphere, z-cylinder and PLY sphere; the inside fraction must stay at the volume
    # fraction at 10 and 50 ms.
    add_test(
        NAME permeability_coarse_step
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/test_permeability_coarse_step.py"
                --bin $<TARGET_FILE:MC-DC_Simulator>)
    set_tests_properties(permeability_coarse_step PROPERTIES LABELS "validation;permeability" TIMEOUT 900)

    # Population-history features (CPU, self-contained): the four-population dynamical-history split
    # (`crossed_filter` partitions the walkers; crossed channels empty without permeability)
    # and explicit radii injection (`radii_file`: originals == injected list, replicas flagged).
    add_test(
        NAME population_history
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/test_population_history.py")
    set_tests_properties(population_history PROPERTIES LABELS "population" TIMEOUT 300)
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

    # Logging / extra-output options (write_bin, trajectories, <log> positions/propagator,
    # log_phase_shift, log_opp, separate_signals, compute_volume), each checked against a known
    # answer in an SI .conf.
    add_test(
        NAME log_outputs
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/test_log_outputs.py"
                --bin $<TARGET_FILE:MC-DC_Simulator>)
    set_tests_properties(log_outputs PROPERTIES LABELS "validation" TIMEOUT 300)

    # Walker-initialisation options (<delta>, ini_walkers_file, <regular_sampling>, <spawning_area>,
    # intra/extra on spheres and on a concave mesh) and the rejected combinations.
    add_test(
        NAME init_options
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/test_init_options.py"
                --bin $<TARGET_FILE:MC-DC_Simulator>)
    set_tests_properties(init_options PROPERTIES LABELS "validation" TIMEOUT 300)

    # Substrate and run-control options: SI vs use_mm_ms equivalence (lists, PLY lists, positions
    # file), ply / ply_file_list / ply_file_list_scale_permeability, periodic and symmetric
    # boundaries, oriented_cylinders_list, max_sim_time, verbatim, coll_sphere_dist.
    add_test(
        NAME substrate_options
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/test_substrate_options.py"
                --bin $<TARGET_FILE:MC-DC_Simulator>)
    set_tests_properties(substrate_options PROPERTIES LABELS "validation" TIMEOUT 300)

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
    # Analytic validation tier (tests/validation): converged-N runs against analytical oracles
    # (free exp(-bD), van Gelderen / Neuman GPA, Karger exchange, T2 closed form, packings, a
    # real cell mesh, diffusion-time dependence, the bounding box). Part of the default suite.
    # Standalone scripts (numpy only); --gpu 0 is the CPU path. Same list as the GPU line.
    foreach(vcase bounding_box cross_representation cylinder_vangelderen free_exp karger_exchange
                  mesh_cell packing_exchange packing_restricted sphere_neuman t2_decay time_dependence)
        add_test(
            NAME validation_${vcase}
            COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/validation/test_${vcase}.py"
                    --bin $<TARGET_FILE:MC-DC_Simulator> --gpu 0
            WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
        set_tests_properties(validation_${vcase} PROPERTIES LABELS "validation_analytic;tier3_analytic" TIMEOUT 1200)
    endforeach()

    # Gamma packing generator: achieved ICVF, radius distribution, non-overlap, auto voxel, determinism.
    add_test(
        NAME validation_gamma_packing
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/validation/test_gamma_packing.py"
                --bin $<TARGET_FILE:MC-DC_Simulator>
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
    set_tests_properties(validation_gamma_packing PROPERTIES LABELS "validation_analytic;tier3_analytic" TIMEOUT 600)

    # Cross-version tier: the current build must reproduce an untouched historical release
    # (v1.5 / v2.0 built from this repo's own history by tests/tools/build_old_version.py and
    # cached under tests/crossversion/old_binaries/) on the same experiment, within Monte Carlo
    # noise. The first run compiles the old release (a few minutes); later runs reuse the cache.
    foreach(xv v1.5 v2.0)
        add_test(
            NAME crossversion_free_diffusion_${xv}
            COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/tests/test_crossversion.py"
                    --bin $<TARGET_FILE:MC-DC_Simulator> --version ${xv} --case free_diffusion
            WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
        set_tests_properties(crossversion_free_diffusion_${xv} PROPERTIES
            LABELS "crossversion;tier4_crossversion" TIMEOUT 1800 RUN_SERIAL TRUE)
    endforeach()
else()
    message(WARNING "Python3 not found: skipping the rng_reproducibility and accuracy_golden tests.")
endif()
endif()  # tests/ present (development tree)
