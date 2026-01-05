# THEMOSIA Winter School 2026  
## Coarse-Grained Modeling of Solid–Ionic Liquid Interfaces with Martini 3

This repository contains the **hands-on course material** for the **THEMOSIA Winter School 2026**, focused on coarse-grained (CG) molecular modeling of **solid–ionic liquid interfaces** using the **Martini 3 force field**.

The course combines **lectures** and **hands-on tutorials**, guiding participants through the construction, simulation, and analysis of CG models for interfaces involving **silica, graphene, and imidazolium-based ionic liquids**.

---

## Course Information

- **Event:** THEMOSIA Winter School 2026
- **Location:** ENS de Lyon (Monod site), Lyon, France
- **Dates:** January 19–23, 2026
  
- **Coarse-grained Module:** January 22–23, 2026
     
- **Instructor:**  
  - Dr. Paulo Cesar Telles de Souza (ENS de Lyon, CNRS, France)

The coarse-grained module follows introductory molecular dynamics lectures earlier in the week.

## Course Scope and Objectives

The main objective of this course is to provide participants with **practical experience** in coarse-grained molecular simulations of interfaces.

Topics covered include:

- Fundamental principles of **coarse-graining**
- Overview and philosophy of the **Martini 3 force field**
- Modeling **solid surfaces** (silica and graphene) in CG simulations
- Coarse-grained representations of **imidazolium-based ionic liquids**
  - C2, C4, and C8 alkyl chain lengths
- Construction of **solid–ionic liquid interfaces**
- Running CG molecular dynamics simulations
- Analysis of **interfacial structure and dynamics**
- Interpretation and limitations of CG results

The emphasis is on **model construction, physical insight, and critical interpretation**, rather than force-field development.

---

## Hands-on Tutorials Overview

The hands-on material is organized into a sequence of tutorials located in the `hands_on/` directory.

### Tutorial 1 – Building Solid Surfaces
- Construction of **silica slabs** using a custom CG builder
- Exploration of surface parameters:
  - bead types
  - surface charges
  - rigidity and restraints
- Construction of **graphene surfaces** using Martini 3 models

### Tutorial 2 – Ionic Liquids in Coarse-Grained Representation
- Imidazolium-based ionic liquids with different alkyl chain lengths:
  - C2
  - C4
  - C8
- Assembly of bulk ionic liquid systems
- Discussion of mapping choices and CG limitations

### Tutorial 3 – Building Solid–Ionic Liquid Interfaces
- Combining solid slabs with ionic liquids
- Introduction of vacuum / gas regions
- Boundary conditions and simulation setup

### Tutorial 4 – Running CG Molecular Dynamics
- Running Martini 3 simulations using:
  - **GROMACS**
  - **OpenMM**
- Differences between CG and all-atom MD simulations
- Practical considerations (time step, thermostats, restraints)

### Tutorial 5 – Analysis of Interfacial Structure and Dynamics
- Partial density profiles
- Radial distribution functions (RDFs)
- Diffusion coefficients and mean-square displacement (MSD)
- Surface-induced structural organization

### Bonus Tutorial – From Atomistic to CG Potentials
- Conceptual overview of **iterative Boltzmann inversion (IBI)**
- Tabulated potentials derived from atomistic reference data
- Discussion of transferability and limitations

## Repository Organization

```
themosia2026-cg-interfaces/
├── docs/                 # Course documentation
├── environment/          # Environment setup (conda, dependencies)
├── external/             # External codes and submodules
├── hands_on/
│   ├── 00_templates/     # Models, scripts, and templates (do not edit)
│   ├── 01_build_slabs/
│   ├── 02_build_interface/
│   ├── 03_run_simulations/
│   └── 04_analysis/
└── README.md
```

Participants are expected to **build their own systems** during the tutorials.

The `00_templates/` directory contains reference scripts, models, and templates that should be **copied or linked**, not modified directly.

---

## Key References

Participants are encouraged to consult the following references:

- Souza et al., *Nature Methods* (2021): Martini 3 force field  
- Padua et al., *Green Chemistry* (2020): Ionic liquids and interfaces  
- Relevant JCTC and interface modeling literature

A detailed list of references is provided in `docs/references.md`.

---

## Key Codes and Resources

- **Martini 3 Force Field**
- **GROMACS**
- **OpenMM**
- Martini 3 graphene models
- OpenMM implementation of Martini
- Atomistic IL–surface reference workflows

Installation instructions and links are provided in the documentation.

---

## Notes for Participants

- These materials are intended **for educational purposes**
- Coarse-grained models have **known limitations**
- Results should be interpreted with care
- If you use these workflows or models in research, please **cite the original authors and force-field developers**

---

## Acknowledgements

This course material was prepared for the **THEMOSIA Winter School 2026** and builds upon contributions from the Martini community and the molecular simulation community at large.
