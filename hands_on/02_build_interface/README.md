# Hands-on 02 — Building Solid–Ionic Liquid Interfaces

## Objective

In this hands-on, we build **solid–ionic liquid interfaces** using coarse-grained (CG) models within the Martini 3 framework. Starting from the solid surfaces constructed in Hands-on 01 (silica or graphite), we add Martini 3 models of imidazolium-based ionic liquids and prepare systems suitable for molecular dynamics simulations.

The goal is to investigate **how solid surfaces influence the organization of ionic liquids**, focusing on interfacial layering, density profiles, and structural ordering.

---

## Background and modeling choices


The coarse-grained models of ionic liquids used in this hands-on are based on the Martini 3 models developed and validated by [Vazquez-Salazar et al. Green Chem., 2020](https://doi.org/10.1039/D0GC01823F). This model was later updated, with a new version released with the major publication of Martini 3 force-field, [P.C.T. et al.Nature Methods 2021](https://doi.org/10.1038/s41592-021-01098-3).These models describe a family of **imidazolium-based ionic liquids** of the form **[C<sub>n</sub>mim][BF<sub>4</sub>]**, where **C<sub>n</sub>mim<sup>+</sup>** denotes the *1-alkyl-3-methylimidazolium* cation with an alkyl chain of length *n* (with *n* = 2, 4, or 8 in this course), and **BF<sub>4</sub><sup>−</sup>** is the *tetrafluoroborate* anion.

These coarse-grained models reproduce key bulk properties of imidazolium-based ionic liquids, including density, nanostructural heterogeneity, and chain-length-dependent ordering (Figure 1). In this course, we exploit these validated bulk properties to investigate how solid surfaces further perturb ionic-liquid organization at solid–ionic liquid interfaces.

  <p align="center">
  <img src="figures/structure_IL.png" width="700">
  </p>

<p align="center">
<em><strong>Figure 1</strong>: Structural characterization  of Martini 3 CG models of ionic liquids. The model reasinably reproduce experimental macroscopic phases(A)  and some features of local structural observed in atomics models (B). Figure is adapted from Vazquez-Salazar et al. Green Chem., 2020. </em>
</p>

### Ionic liquids considered

We focus on **imidazolium-based ionic liquids** paired with the **BF<sub>4</sub><sup>−</sup>** anion, using two representative alkyl chain lengths on the cation:

- [C<sub>2</sub>mim][BF<sub>4</sub>] or **C2**: 1-ethyl-3-methylimidazolium tetrafluoroborate  
- [C<sub>8</sub>mim][BF<sub>4</sub>] or **C8**: 1-octyl-3-methylimidazolium tetrafluoroborate  

These two systems span increasing degrees of **nanostructural heterogeneity** in the bulk ionic liquid, while **not forming lamellar phases under the conditions considered**. Lamellar ordering typically emerges only for longer alkyl chains (e.g. C12).

This makes C2 and C8 ideal model systems to explore whether **solid surfaces can induce additional ordering or layering** in ionic liquids that is absent in the bulk.

As a bonus (if time allows), the `00_templates/` directory also includes coarse-grained models for **C4 and C12 imidazolium-based ionic liquids**, as well as a **phosphonium-based cations**(specifically trihexyltetradecylphosphonium), which can be explored for comparison or extended studies.


### Coarse-grained ionic liquid model

The coarse-grained ionic liquid models used here were built based on the following key features:

- Careful **mapping and bonded geometry choices** to reproduce molecular volume, solvent-accessible surface area (SASA), and packing behavior.
- Coupled with the previous point, careful parametrization of the bonded terms was also considered in order to reproduce distributions of bond distances, angles, and dihedrals obtained from atomistic simulations, such as those discussed in the previous days of this school.
- Compatibility with the Martini 3 interaction matrix, respecting the expected bead assignments according to the chemistry being represented.
- Use of **partial charges (q labels)** to represent charge delocalization on the imidazolium ring. The values of the charges in the imidazolium ring (+0.5 on each nitrogen-based bead) were computed according to a variant of the Dipole Preserving Charge method, based on quantum chemical calculations similar to those discussed on the first day of this school. This allowed for an improved electrostatic description, including possible quadrupolar effects. In this sense, simple dipole and quadrupole estimates can be used to help define an appropriate set of partial charges for the beads representing the imidazolium ions.

<p align="center">
  <img src="figures/model_IL.png" width="900">
</p>

<p align="center">
<em><strong>Figure 2</strong>: Martini 3 CG models of ionic liquids. (A) CG models of trihexyltetradecylphosphonium and 1,3-dialkylimidazolium
cations, as well as tetrafluoroborate anions. The Martini bead types and sizes are indicated. Blue indicates positively charged groups, while red and gray indicate negatively charged and nonpolar groups, respectively. (B) Molecular surfaces (also called Connolly surfaces) of atomistic and CG structures of the C2 cation and the [BF<sub>4</sub>]<sup>−</sup> anion. Figure adapted from Vazquez-Salazar et al., Green Chem. (2020) and Souza et al., Nat. Methods (2021).</em>
</p>

As the focus of this hands-on is **not coarse-grained parametrization**, but rather **interfacial organization**, we will not extend this discussion further. Extensive examples of parametrization strategies are available in the tutorials provided by the Martini Force Field Initiative website.

---

## Interface construction strategy

To keep the setup robust and reproducible, we adopt the following strategy:

- The **solid slab** (silica or graphite) is kept fixed during construction
- Ionic liquids are inserted **above the surface** using **PACKMOL**
- **Ion pairs (cation + anion)** are used during packing to guarantee charge neutrality
- The simulation geometry is kept **identical for all systems**, enabling direct comparison

### Fixed geometry for all systems

- Solid slab: 10 × 10 × 2 nm
- Ionic liquid region: 15 nm
- Vacuum region: 10 nm
- Total box height: 20 nm
- Periodic boundary conditions in x and y

---

## Step 1 — Extract surface dimensions

Before building the interface, extract the lateral box dimensions from the surface file generated in Hands-on 01:

```bash
tail -n 1 silica.gro
```

The x and y values will be reused to define the interface box, ensuring perfect lateral matching between surface and ionic liquid.

---

## Step 2 — PACKMOL input: slab + ionic liquid + vacuum

PACKMOL is used to place ionic liquid **ion pairs** above the surface while leaving a vacuum region empty.

Because Martini beads are larger and softer than atomistic particles, we use a **larger tolerance** than typical atomistic setups.

### Example PACKMOL input (C4–BF₄)

```text
tolerance 2.8
filetype pdb
output silica_IL_vacuum.pdb
seed 12345

structure silica.pdb
  fixed 0.0 0.0 0.0   0.0 0.0 0.0
end structure

structure C4-BF4.pdb
  number 4000
  inside box 0.0 0.0 25.0   100.1700 99.1426 170.0
end structure
```

Notes:

- The slab occupies approximately z = 0–20 Å
- Ionic liquid starts slightly above the surface (z = 25 Å) to avoid overlaps
- The region above z = 170 Å is left empty and acts as vacuum
- Charge neutrality is guaranteed by packing **ion pairs**

Run PACKMOL:

```bash
packmol < input.inp
```

---

## Step 3 — Define the final simulation box

Convert the PACKMOL output to a GROMACS `.gro` file and define the final box dimensions:

```bash
gmx editconf -f silica_IL_vacuum.pdb -o box.gro -box 10.01700 9.91426 20.00000 -noc
```

The x and y values must match those extracted from the surface file.

---

## Step 4 — Reorder atoms for GROMACS topology

PACKMOL does not guarantee molecule ordering. For GROMACS, we reorder the `.gro` file so that atoms appear in the following order:

1. Silica slab  
2. Cations  
3. Anions  

Residue names used:

- Silica: 1SIS
- Cation: 1BIM
- Anion: 2BF4

### Reordering procedure

```bash
title=$(head -n 1 box.gro)
box=$(tail -n 1 box.gro)

tail -n +3 box.gro | head -n -1 > atoms_only.tmp

{
  grep "^[[:space:]]*1SIS" atoms_only.tmp
  grep "^[[:space:]]*1BIM" atoms_only.tmp
  grep "^[[:space:]]*2BF4" atoms_only.tmp
} > atoms_reordered.tmp

natoms=$(wc -l < atoms_reordered.tmp)

{
  echo "$title"
  printf "%5d\n" "$natoms"
  cat atoms_reordered.tmp
  echo "$box"
} > start.gro

rm atoms_only.tmp atoms_reordered.tmp
```

---

## Step 5 — Topology file

Example topology file:

```text
#include "martini_v3.0.itp"
#include "cation_C4C1imin.itp"
#include "anion.itp"
#include "silica.itp"

[ system ]
Interface silica C4-BF4

[ molecules ]
SILICA_SLAB 1
BIM   4000
BF4   4000
```

---

## Final remarks and discussion points

- The number of ion pairs is chosen to give a **reasonable filling** of the ionic liquid region. At the coarse-grained level, this is a modeling choice rather than a strict physical constraint.
- The same geometry and protocol are used for C2, C4, and C8, enabling direct comparison.
- In production studies, the number of ion pairs could be adjusted to reproduce bulk densities; this is not required for this hands-on.

### Questions for discussion

- How would this protocol need to be adapted for **charged silica surfaces**?
- What changes (if any) are required when replacing silica with **graphite**?
- How might surface chemistry influence ionic liquid layering and orientation?
