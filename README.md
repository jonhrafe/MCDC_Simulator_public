![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/jonhrafe/MCDC_Simulator_public)
![GitHub issues](https://img.shields.io/github/issues/jonhrafe/MCDC_Simulator_public)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/jonhrafe/MCDC_Simulator_public/graphs/commit-activity)
![GitHub last commit](https://img.shields.io/github/last-commit/jonhrafe/MCDC_Simulator_public)
![GitHub top language](https://img.shields.io/github/languages/top/jonhrafe/MCDC_Simulator_public)
![GitHub](https://img.shields.io/github/license/jonhrafe/MCDC_Simulator_public)

<p align="center">
  <a href="https://www.frontiersin.org/articles/10.3389/fninf.2020.00008/">
    <img src="https://user-images.githubusercontent.com/4105920/68854670-d2f40280-06dc-11ea-8b45-9253fb6eec41.png" alt="MC/DC logo" width="150" height="240">
  </a>
</p>
<p align="center">
  User Friendly, Fast and Robust <strong>M</strong>onte <strong>C</strong>arlo<strong> D</strong>iffusion and <strong>C</strong>ollision Simulator
  <br>
  <a href="https://github.com/jonhrafe/MCDC_Simulator_public/tree/master/instructions"><strong>Explore MC/DC docs and examples »</strong></a>
  <br>
  <br>
  <a href="https://github.com/jonhrafe/MCDC_Simulator_public/issues">Report bug</a>
  ·
  <a href="https://github.com/jonhrafe/MCDC_Simulator_public/issues">Request feature</a>
</p>


## Table of contents
- [Introduction](#introduction)
- [Features](#features)
- [Quick start](#quick-start)
- [Units](#units)
- [Documentation](#documentation)
- [Bugs and feature requests](#bugs-and-feature-requests)
- [Versioning](#versioning)
- [Developers](#developers)
- [Copyright and license](#copyright-and-license)

## Introduction
The <strong>M</strong>onte <strong>C</strong>arlo<strong> D</strong>iffusion and <strong>C</strong>ollision simulator (MC/DC) is a C++ open-source **Diffusion-Weighted Magnetic Resonance Imaging** (DW-MRI) Monte-Carlo simulator. For an in-detail explanation of the numerical framework and the basics of DW-MRI we refer the visitor to the following publication: [https://doi.org/10.3389/fninf.2020.00008](https://doi.org/10.3389/fninf.2020.00008)

### Citation:

> AUTHOR=Rafael-Patino Jonathan, Romascano David, Ramirez-Manzanares
> Alonso, Canales-Rodríguez Erick Jorge, Girard Gabriel, Thiran
> Jean-Philippe
> TITLE=Robust Monte-Carlo Simulations in Diffusion-MRI:
> Effect of the Substrate Complexity and Parameter Choice on the
> Reproducibility of Results
> JOURNAL=Frontiers in Neuroinformatics
> VOLUME=14       YEAR=2020 PAGES=8
> URL=https://www.frontiersin.org/article/10.3389/fninf.2020.00008
> DOI=10.3389/fninf.2020.00008     ISSN=1662-5196

## Features

MC/DC diffuses spin packets through a user-defined substrate and synthesizes the resulting
DW-MRI signal. The only external dependency is the C++ **Eigen** template library for linear
algebra, which is bundled in the repository (`src/Eigen`); everything else — from the geometry
and collision handling to the signal synthesis — is built from scratch.

- **Substrates:** free diffusion; analytic **spheres** and **cylinders**; triangulated **PLY meshes**
  (single or multiple); **hexagonal packings**; and **gamma-distributed** spheres and cylinders.
- **Microstructure physics:** per-compartment **diffusivity**, **T2 relaxation**, and membrane
  **permeability** (Powles crossing), with intra/extra compartment tracking. Per-obstacle properties
  can be assigned from the geometry lists.
- **Boundary conditions:** periodic and mirror voxel walls (continuous real position for the
  signal/MSD, wrapped position for collisions).
- **Sequences:** Pulsed-Gradient Spin-Echo (**PGSE** / Stejskal–Tanner) and **general gradient
  waveforms**, read from a scheme file.
- **Outputs:** the complex DW-MRI signal (real / imaginary), per-compartment **separated signals**,
  optional **voxel-subdivision** volumetric signal (for parametric maps / movies), per-spin
  **trajectories**, and the diffusion **propagator**.
- **Reproducible & tested:** a single seeded RNG makes runs reproducible, and a CTest suite guards
  the physics against golden references.

All components are organized as class-oriented modules with abstract base prototypes, so the
substrates, sequences and outputs can be extended by inheritance.

## Quick start

 - **[Build the simulator](instructions/compilation.md)** — CMake (recommended) or a single-command build.
 - **[Getting started](instructions/GettingStarted.md)** — run your first (free-diffusion) simulation and learn the `.conf` parameters.
 - **[Tutorial: gamma-distributed cylinders](instructions/example_intra-axonal_initialization.md)** — intra-axonal initialization on a generated substrate.

Once built, a simulation is launched with a single configuration file:

```bash
./MC-DC_Simulator docs/conf_file_examples/freeDiffusion.conf
```

## Units

By default a `.conf` file (and its scheme file) is written in **standard SI units**: metres,
seconds and Tesla. Values are scaled silently to the internal working units (mm, ms) on load.
This applies to **all** lengths, including substrate-generation parameters (voxel, sampling area,
sphere/cylinder/hex radii, gamma `beta`/`min_radius`, and geometry-file scale factors, which are
"metres per file unit"). The only exception is permeability (a velocity: m/s equals mm/ms
numerically, so it is scale-invariant).

Set `use_mm_ms 1` to declare that a file is already in the internal units (mm, ms) and skip the
scaling. The legacy `scale_from_stu` flag is still accepted (it maps to the inverse of `use_mm_ms`)
but is deprecated.

## Documentation

Per-parameter and per-block reference lives in the [Getting started](instructions/GettingStarted.md)
page and in the commented example configurations under
[`docs/conf_file_examples/`](docs/conf_file_examples/). Source-level API documentation can be
generated with **doxygen** from the in-source comments.

## Bugs and feature requests

Have a bug or a feature request? Please search the [existing and closed issues](https://github.com/jonhrafe/MCDC_Simulator_public/issues) first. If your problem or idea is not addressed yet, [open a new issue](https://github.com/jonhrafe/MCDC_Simulator_public/issues).

## Versioning

The first version, **1.42**, is the one released with the paper and is preserved for full
reproducibility of the published results at
[https://github.com/jonhrafe/Robust-Monte-Carlo-Simulations](https://github.com/jonhrafe/Robust-Monte-Carlo-Simulations).
Subsequent **2.x** releases add features and changed the configuration conventions — most notably,
**`.conf` files are now in SI units by default** (see [Units](#units)). Configuration files written
for 1.42 may need updating; the bundled examples reflect the current conventions.

## Developer(s)

**Jonathan Rafael-Patino**
- [https://people.epfl.ch/jonathan.patinolopez?lang=en](https://people.epfl.ch/jonathan.patinolopez?lang=en)
- [https://www.linkedin.com/in/jonhrafe/](https://www.linkedin.com/in/jonhrafe/)

## Copyright and license

**GNU Lesser General Public License v2.1**

Primarily used for software libraries, the GNU LGPL requires that derived works be licensed under the same license, but works that only link to it do not fall under this restriction.
