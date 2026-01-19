# Hands-on Sessions Overview

## Before you start

### Where are the input files?
All reference models, scripts, and template files are available in:

- [`hands_on/00_templates/`](00_templates/)

This folder contains force-field files, molecule structures, MDP templates, and helper scripts.  
In general, you should **not edit files inside `00_templates/`**.

### Recommended workflow (one folder per system)
To keep things organized, create a separate working directory for each interface you build. For instance

```bash
mkdir -p work
cd work
mkdir silica_C2mimi-BF4

# Then along the hands-on, copy template files you need into your working folder. For instance, to copy XXX.pdb
cp -r ../00_templates/pdb/XXX.pdb silica_C2mimi-BF4
cd silica_C2mimi-BF4

> If you are not familiar with Linux commands (`cd`, `ls`, `cp`, `mkdir`), please ask the instructors for help.

Now, checkout the overview of each Hands-on software enviroment.

## [-Hands-on 01 – Building solid surfaces](01_build_slabs/)
- Construction of Martini 3 models of **silica slabs** using a custom CG builder
- Exploration of surface parameters:
  - chesmitry and bead types
  - surface density
  - surface charge
- Construction of **graphite slabs** using Martini 3 graphene models

## [-Hands-on 02 – Building Solid–Ionic Liquid Interfaces](02_build_interface/)
- Imidazolium-based ionic liquids with different alkyl chain lengths.
- Combining solid slabs with ionic liquids
- Introduction of vacuum / gas region
- Simulation Boxes and topology files

## [-Hands-on 03 – Running CG Molecular Dynamics](03_run_md_and_analysis/)
- Running Martini 3 simulations using GROMACS
- Differences between CG and all-atom MD simulations
- Practical considerations (time step, thermostats, restraints)
- Analysis
  - Partial density profiles

---

# Software environment
The hands-on sessions were tested on Debian Linux systems with the following
software available system-wide:

- [GROMACS](https://www.gromacs.org/) (2023.3 or newer)
- [VMD](https://www.ks.uiuc.edu/Research/vmd/])
- [Packmol](https://m3g.github.io/packmol/)
- [Python](https://www.python.org/)
- [MDAnalysis](https://www.mdanalysis.org/)

> **Note:** All tools are already installed in the computer, except for Packmol and MDAnalysis, You may use them via a dedicated Python virtual environment.
> You can activate it with:
>
> ```bash
> source /projects/themosia/CG/virtual_env/bin/activate
> ```

---



