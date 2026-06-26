[![Build](https://github.com/jonhrafe/MCDC_Simulator_public/actions/workflows/build.yml/badge.svg?branch=mcdc2_dev)](https://github.com/jonhrafe/MCDC_Simulator_public/actions/workflows/build.yml)
[![Latest release](https://img.shields.io/github/v/release/jonhrafe/MCDC_Simulator_public)](https://github.com/jonhrafe/MCDC_Simulator_public/releases/latest)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/jonhrafe/MCDC_Simulator_public)
![GitHub issues](https://img.shields.io/github/issues/jonhrafe/MCDC_Simulator_public)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/jonhrafe/MCDC_Simulator_public/graphs/commit-activity)
![GitHub last commit](https://img.shields.io/github/last-commit/jonhrafe/MCDC_Simulator_public/mcdc2_dev)
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
  <a href="https://github.com/jonhrafe/MCDC_Simulator_public/tree/mcdc2_dev/instructions"><strong>Explore MC/DC docs and examples »</strong></a>
  <br>
  <br>
  <a href="https://github.com/jonhrafe/MCDC_Simulator_public/issues">Report bug</a>
  ·
  <a href="https://github.com/jonhrafe/MCDC_Simulator_public/issues">Request feature</a>
</p>

> [!IMPORTANT]
> **This is MC/DC `v2.2.0` — the current stable release**, now on the default branch (`mcdc2_dev`): SI units, the modernized PGSE protocol, multi-compartment substrates, four tutorials, `pip install`, and cross-platform binaries on the [Releases page](https://github.com/jonhrafe/MCDC_Simulator_public/releases/latest).
> **Looking for the previous stable version?** It remains available on the [**`master`**](https://github.com/jonhrafe/MCDC_Simulator_public/tree/master) branch.

## Table of contents
- [Introduction](#introduction)
- [Quick start](#quick-start)
- [Status](#status)
- [Simulations](#whats-included)
- [Bugs and feature requests](#bugs-and-feature-requests)
- [Documentation](#documentation)
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

## Quick start
Several quick start options are available:

 - [Download a pre-compiled version if available](https://github.com/jonhrafe/MCDC_Simulator_public/releases)
 - **Install with pip** — `pip install .` builds the simulator and adds the `mcdc` command ([details](instructions/compilation.md#easiest-install-with-pip))
 - [Compile the sources](instructions/compilation.md) with CMake
 - Read the [Getting started page](instructions/GettingStarted.md) for information on the basic parameters needed.
 - [Tutorial: Simulation in free diffusion](instructions/tutorial_free_diffusion.md)
 - [Tutorial: Simulation in gamma distributed cylinders](instructions/tutorial_gamma_cylinders.md)
 - [Tutorial: Simulation in PLY models](instructions/tutorial_ply_meshes.md)
 - [Tutorial: Multi-compartment substrate (multi-D, multi-T2, multi-PLY)](instructions/tutorial_multi_compartment.md)

> **Note on units.** A `.conf` file (and its scheme) is in **standard SI units** by default — metres,
> seconds, Tesla — scaled internally to mm/ms on load. This applies to all lengths, including the
> substrate-generation scales (`ply_scale` and the sphere/cylinder list scale are "metres per file
> unit", and gamma `beta`/`min_radius` are in metres). Set `use_mm_ms 1` to declare a file already in
> the internal units; the legacy `scale_from_stu` flag is still accepted but deprecated.

## Status
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/jonhrafe/MCDC_Simulator_public/graphs/commit-activity)
![GitHub last commit](https://img.shields.io/github/last-commit/jonhrafe/MCDC_Simulator_public)
![GitHub issues](https://img.shields.io/github/issues/jonhrafe/MCDC_Simulator_public)
[![ForTheBadge built-with-science](http://ForTheBadge.com/images/badges/built-with-science.svg)](https://www.frontiersin.org/articles/10.3389/fninf.2020.00008/full)

## What's included

The MC/DC Simulator uses the C++ **Eigen** template library for linear algebra as its only external
dependency, which is bundled in the repository (`src/Eigen`); all other components — from the
3D-mesh handling to the MRI signal synthesis — are built from scratch as class-oriented modules with
abstract base prototypes that can be inherited and re-implemented to extend the simulator.

MC/DC diffuses spin packets through a user-defined substrate and synthesizes the resulting DW-MRI
signal. It supports:

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

## Bugs and feature requests

Have a bug or a feature request? Please search the [existing and closed issues](https://github.com/jonhrafe/MCDC_Simulator_public/issues) first. If your problem or idea is not addressed yet, [open a new issue](https://github.com/jonhrafe/MCDC_Simulator_public/issues).

## Documentation

Per-parameter and per-block reference lives in the [Getting started](instructions/GettingStarted.md)
page and in the commented example configurations under
[`docs/conf_file_examples/`](docs/conf_file_examples/). Source-level API documentation can be
generated with **doxygen** from the in-source comments.

## Versioning

The first version, **1.42**, is the one released with the paper and is preserved for full
reproducibility of the published results at
[https://github.com/jonhrafe/Robust-Monte-Carlo-Simulations](https://github.com/jonhrafe/Robust-Monte-Carlo-Simulations).
Subsequent **2.x** releases add features and changed the configuration conventions — most notably,
**`.conf` files are now in SI units by default** (see the note under [Quick start](#quick-start)).
Configuration files written for 1.42 may need updating; the bundled examples reflect the current
conventions.

## Developer(s)

**Jonathan Rafael-Patino**
- [https://people.epfl.ch/jonathan.patinolopez?lang=en](https://people.epfl.ch/jonathan.patinolopez?lang=en)
- [https://www.linkedin.com/in/jonhrafe/](https://www.linkedin.com/in/jonhrafe/)

## Copyright and license

**GNU Lesser General Public License v2.1**

Primarily used for software libraries, the GNU LGPL requires that derived works be licensed under the same license, but works that only link to it do not fall under this restriction.
