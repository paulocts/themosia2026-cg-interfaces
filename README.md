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


## [Hands-on Sessions Overview](hands_on/)

The hands-on material is organized into a sequence of tutorials located in the `hands_on/` directory.

### - [Hands-on 01 – Building solid surfaces]

### - [Hands-on 02 – Ionic Liquids in Coarse-Grained Representation]

### - [Hands-on 03 – Building Solid–Ionic Liquid Interfaces]

### - [Hands-on 04 – Running CG Molecular Dynamics]

### - [Hands-on 05 –  Analysis of Interfacial Structure and Dynamics]

## Repository Organization

```
themosia2026-cg-interfaces/
├── docs/                  # Course documentation
├── environment/           # Environment setup (conda, dependencies)
├── external/              # External codes and submodules
├── hands_on/
│   ├── 00_templates/      # Models, scripts, and templates (do not edit)
│   ├── 01_build_slabs/    # Building solid surfaces (silica, graphene)
│   ├── 02_ionic_liquids/  # Bulk ionic liquids (C2, C4, C8)
│   ├── 03_build_interface/# Solid–ionic liquid interfaces
│   ├── 04_run_simulations/# Running CG MD (GROMACS / OpenMM)
│   └── 05_analysis/       # Analysis of interfacial properties
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

## Notes for Participants

- These materials are intended **for educational purposes**
- Coarse-grained models have **known limitations**
- Results should be interpreted with care
- If you use these workflows or models in research, please **cite the original authors and force-field developers**

---

## Acknowledgements

This course material was prepared for the **THEMOSIA Winter School 2026** and builds upon contributions from the Martini community and the molecular simulation community at large.
