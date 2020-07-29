![GitHub commits since latest release (by date)](https://img.shields.io/github/commits-since/jonhrafe/MCDC_Simulator_public/1.42)
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
  <a href="https://github.com/jonhrafe/MCDC_Simulator_public/projects/1">Request feature</a>
</p>


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
The <strong>M</strong>onte <strong>C</strong>arlo<strong> D</strong>iffusion and <strong>C</strong>ollision simulator (MC/DC), is an C++ open source **Diffusion Weighted Magnetic Resonance Imaging** (DW-MRI) Monte Carlo Simulator. For a in-detail explanation of the numerical framework and basics on DW-MRI imaging we refer the visitor to the following publication:  <img src="https://www.frontiersin.org/Areas/Home/Content/Images/logo-home.svg" alt="Frontiers" width="100" height="20"></a>  [https://doi.org/10.3389/fninf.2020.00008](https://doi.org/10.3389/fninf.2020.00008 ) 

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

 - [Download a pre-compiled version if available](https://github.com/jonhrafe/Robust-Monte-Carlo-Simulations/releases)
 - [Compile the sources](https://github.com/jonhrafe/MCDC_Simulator_public/blob/master/instructions/compilation.md)
 - Read the [Getting started page](https://github.com/jonhrafe/MCDC_Simulator_public/blob/master/instructions/GettingStarted.md) for information on the basic parameters needed.
 - [Tutorial: Simulation in free diffusion](https://github.com/jonhrafe/MCDC_Simulator_public/blob/master/instructions/GettingStarted.md)
 - [Tutorial: Simulation in gamma distributed cylinders](https://github.com/jonhrafe/MCDC_Simulator_public/blob/master/instructions/example_intra-axonal_initialization.md)
 - [Tutorial: Simulation in PLY models](https://github.com/jonhrafe/MCDC_Simulator_public/blob/master/instructions/example_intra-axonal_initialization.md)

## Status
![GitHub commits since latest release (by date)](https://img.shields.io/github/commits-since/jonhrafe/MCDC_Simulator_public/1.42)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/jonhrafe/MCDC_Simulator_public)
![GitHub issues](https://img.shields.io/github/issues/jonhrafe/MCDC_Simulator_public)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/jonhrafe/MCDC_Simulator_public/graphs/commit-activity)
![GitHub last commit](https://img.shields.io/github/last-commit/jonhrafe/MCDC_Simulator_public)
[![ForTheBadge built-with-science](http://ForTheBadge.com/images/badges/built-with-science.svg)](https://www.frontiersin.org/articles/10.3389/fninf.2020.00008/full)

## What's included

The MC/DC Simulator only uses C++ Eigen template library for linear algebra as an external dependency, which is included in the repository; all other basic structures ---from 3d-mesh-models handling to the MRI signal synthesis---are built from scratch in the following hierarchical structure.  

![enter image description here](https://user-images.githubusercontent.com/4105920/68950663-e9718b00-07bc-11ea-8b27-db16733c103c.png)

All the displayed components depicted above are organized in separated class-oriented modules with abstract-based prototypes that can inherited and  re-implemented to augment the functionality and scope of the a simulations. 

## Bugs and feature requests

Have a bug or a feature request? Please first read the [ongoing development chart](https://github.com/jonhrafe/MCDC_Simulator_public/projects/1) and search for existing and closed issues. If your problem or idea is not addressed yet, [please open a new issue](https://github.com/jonhrafe/MCDC_Simulator_public/issues).

## Documentation

The MC/DC technical documentation is automatically generated using **doxygen** and the source-files comments, and can be found [in here](https://github.com/jonhrafe/MCDC_Simulator_public/blob/master/docs/html/index.html) in html format, or [in here](https://github.com/jonhrafe/MCDC_Simulator_public/tree/master/docs/latex) as a .tex source file.

## Versioning

For full reproducibility of the results presented on: [https://www.frontiersin.org/articles/10.3389/fninf.2020.00008/full](https://www.frontiersin.org/articles/10.3389/fninf.2020.00008/full), the first version, 1.42, is the same version released in the paper's GitHub repository:  [https://github.com/jonhrafe/Robust-Monte-Carlo-Simulations](https://github.com/jonhrafe/Robust-Monte-Carlo-Simulations). 

However, **future releases may not follow the conf_file conventions used** and included in the aforementioned repository. 

## Developer(s)

**Jonathan Rafael-Patino**
- [https://people.epfl.ch/jonathan.patinolopez?lang=en](https://people.epfl.ch/jonathan.patinolopez?lang=en)
- [https://www.linkedin.com/in/jonhrafe/](https://www.linkedin.com/in/jonhrafe/)


## Copyright and license

**GNU Lesser General Public License v2.1**

Primarily used for software libraries, the GNU LGPL requires that derived works be licensed under the same license, but works that only link to it do not fall under this restriction. 
