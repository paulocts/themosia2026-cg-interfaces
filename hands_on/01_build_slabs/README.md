# Hands-on 01 — Building Solid Surfaces

## 1. Objective

In this hands-on session you will build a **Martini 3 coarse-grained (CG) silica slab** that will later be used to assemble **solid–ionic liquid interfaces**. The material is based on the work of  [Cambiaso et al 2025](https://doi.org/10.1016/j.surfin.2025.106997).

By the end of this hands-on you should be able to:

- Generate a CG silica slab with a defined size and thickness.
- Create a surface model with **different bead types in the core and at the surface**.
- Control basic “surface properties” that influence interfacial behavior:
  - **bead type** (surface chemistry)
  - **bead spacing / surface density** (structure)
- Produce the files needed for the next hands-ons (coordinates + topology).

---

## 2. Background: CG surfaces in Martini

In coarse-grained simulations, a solid surface can be represented as a **mesoscopic model** designed to reproduce key behaviors (e.g., wetting, layering, adsorption), rather than a one-to-one mapping of every atom. When building the CG silica slabs, we consider two possible simulation box geometries: **orthorhombic** and **hexagonal (triclinic)**. These choices correspond to different modeling intentions.

Hexagonal (triclinic) boxes are the natural choice if one aims to represent a true crystalline α-quartz surface. Such a geometry would be appropriate if the goal were to model a specific quartz crystal face with atomistic fidelity. Orthorhombic boxes, on the other hand, are primarily a **computational convenience**. A triangular lattice can be embedded in an orthorhombic supercell without altering the local bead arrangement or surface bead density. In this case, the box shape does not represent a real crystal symmetry; it simply provides a rectangular container for the coarse-grained surface.

In this context, the term *crystalline* refers to the **local bead arrangement and bead density** used to construct the slab, and **not** to long-range crystallographic symmetry as in atomistic models. Cambiaso et al. (2025) used orthorhombic boxes regardless of whether the surface was described as crystalline or amorphous. The focus of that work—and of the present hands-on—is on reproducing macroscopic interfacial properties (e.g., wetting behavior and contact angles) through the correct surface bead density, rather than on preserving crystallographic symmetry or atom-to-bead mapping.

This means that several choices are **part of the model definition**, not “details”:

- **Surface bead density**, defined by the spacing between beads.
- **Surface bead types**, which define the effective surface chemistry (and possibly surface charge).
- **Surface rigidity**, i.e. whether the surface is treated as rigid, restrained, or flexible.

These choices can strongly affect interfacial structure and dynamics, and they are expected to be tested and justified depending on the target application.

In this course, we will focus on three practical knobs that you can tune and immediately interpret:

1. **Surface bead type (chemistry)**
   - We assign different bead types to the **surface** and the **core** to represent different effective interactions with the liquid (following the modeling logic used for silica surfaces in the reference work). 

2. **Bead spacing / surface density (structure)**
   - The distance between surface beads controls the effective roughness and interaction density seen by the liquid, which can change layering and organization near the interface. 

Once the solid–ionic liquid systems are assembled (later hands-ons), we will quantify surface-induced effects using:

- **partial density profiles** along the surface normal (layering, enrichment/depletion)
- optional: RDFs, diffusion, and orientational order parameters

For now, the goal is to build a clean, reproducible silica slab model with controllable parameters.

---

## 3. Building a silica slab

### 3.1 Overview of the slab model

The slab model used in this tutorial has the following features:

- **Crystalline lattice**  
  Beads are arranged in stacked triangular layers along the surface normal (*z*).

- **Core vs surface beads**  
  - **Core beads** represent the bulk of the solid and provide mechanical stability.
  - **Surface beads** (top and bottom layers) define the effective surface chemistry experienced by the liquid.

- **Bonded solid network**  
  Beads within a cutoff distance are connected by harmonic bonds, forming a mechanically stable solid.

- **Slab-safe periodicity**  
  Bonds are constructed considering periodic boundary conditions **only in the x–y directions**, avoiding unphysical bonded connections across the vacuum gap along *z*.

---

### 3.2 Generating the slab using the builder script

In this hands-on, we generate a **reference silica slab** that closely follows the original model used in the literature. Copy build_silica_slab.py from 00_templates/. You can first run

```bash
python3 build_silica_slab.py -h
```
To see the options of the code avaiable. But to start with an example, let's run:

```bash
python3 build_silica_slab.py \
  -type crystalline \
  -x 10.0 -y 10.0 -z 2.0 \
  -core N2 \
  -surface N5 \
  -bonds -ef 10000 -rcut 0.65 \
  -pbc_bonds xy \
  -pr 10000 \
  -o silica
```

**Key options explained:**

- `-x 10.0 -y 10.0 -z 2.0 `  
  Gives the dimensions of the box in nm, with 2520 beads distributed in 5 planes. As each bead can be considered more or less centered in the Si positions, this would be a crystal equivalent to 7560. So equivalent to a mapping of 3 non-hydrogen atoms represented by 1 bead (roughly like SiO₂) 

- `-core N2` / `-surface N5`  
  Assign different bead types to the bulk and surface to represent effective surface chemistry.

- `-bonds -ef 10000 -rcut 0.65`  
  Create a stiff bonded network between neighboring beads (within 6.5 Å), ensuring the slab behaves as a solid.

- `-pbc_bonds xy`  
  Apply periodic bonding only in the lateral directions, which is safe for slab geometries.

- `-pr 10000`  
  Add **position restraints** to all slab beads (10 000 kJ mol⁻¹ nm⁻²).  
  These restraints are intended **only for energy minimization and equilibration** and are activated by defining `POSRES` in the GROMACS `.mdp` file. Alternatively, you could use position restraints along the production, but without the adding a bonded network.

---

### 3.3 Output files

The builder produces three files:

- **silica.gro**  
  GROMACS coordinate file (units: nm)

- **silica.pdb**  
  Coordinate file in Å, useful for visualization in VMD.

- **silica.itp**  
  Topology file defining bead types, bonds, masses, charges, and optional position restraints

These files will be reused in later hands-ons to assemble solid–ionic liquid systems.

## 4. Modifying surface parameters

Now that you have learned how the builder works, you can generate **controlled variations of surface properties**.  
The examples below illustrate simple and physically motivated modifications that can be explored systematically.

---

### 4.1 Surface bead type (effective chemistry)

A first and intuitive modification is to change the **surface bead type**, which controls the effective chemistry experienced by the liquid.

#### Example: hydrophobized silica surface

Experimentally, silica surfaces are often hydrophobized by trimethylsilylation, replacing surface silanol groups with –Si(CH₃)₃ moieties.

At the CG level, this can be mimicked by assigning a **nonpolar bead type** to the surface:

```bash
-surface C1
```

This represents a hydrophobic surface and is expected to reduce wetting and interfacial layering.

---

### 4.2 Bead spacing and surface density

Bead spacing controls the **effective surface density and roughness** experienced by the liquid at the interface.

By default, the bead spacing used in the silica slab builder corresponds to the values employed by Cambiaso et al. (2025), which were shown to reproduce macroscopic interfacial properties such as wetting behavior.

An alternative modeling choice is to construct a surface whose bead spacing reflects typical **Si–Si distances** in crystalline silica. This can be compared directly with the crystal structures used in the atomistic MD hands-on from the previous day. A representative estimate for this spacing is approximately **0.292 nm**.

The bead spacing can be modified using the `-a` flag in the slab builder:

```bash
-a 0.292
```

When reducing the bead spacing, it is recommended to also adjust the bonding cutoff to avoid creating an excessively connected network. For example, the cutoff distance can be reduced to -rcut 0.53

---

### 4.3 Surface charge and functionalization

Surface charge is introduced in an **effective, coarse-grained manner**, representing the average charge state of surface functional groups.

#### Example: amino-functionalized silica surface

Assuming pH ≈ pKₐ and partial protonation of amine groups, a reasonable CG approximation is to assign a **partial positive charge** to surface beads.

To represent roughly **50% protonation**, use:

```bash
-top N4dq -qtop +0.5
```

This creates a **partially charged surface**, sufficient to induce electrostatic ordering of the ionic liquid without assuming full protonation.

---

### 4.4 Partial and heterogeneous surface modification

Real silica surfaces are rarely perfectly homogeneous.  
To mimic **partial functionalization or mixed surface chemistry**, the builder allows random modification of **only a fraction of surface beads**.

You can define:
- a second bead type on the surface
- the fraction of beads to modify
- optionally, a different charge for the modified beads

If the fraction is not specified, the code defaults to an **equal (50/50) distribution**.

#### Example: 30% amino-functionalized patches on the top surface

```bash
-top2 Q4p -qtop2 +1.0 -ftop2 0.30
```

This produces a surface where 30% of top-layer beads are randomly functionalized, while the remainder retain the original surface chemistry.

This option is useful for exploring:
- heterogeneous functionalization
- imperfect or patchy surface treatments
- sensitivity of interfacial structure to surface disorder

---

## 4.5 Notes and modeling perspective

The default silica slab and the surface modifications introduced here are **model definitions**, not fine-tuning parameters.

They should be chosen based on:
- the physical question being addressed
- consistency with the chosen level of coarse-graining
- systematic comparison across variants

The silica models presented here are best viewed as **toy models**: they sacrifice chemical detail in favor of simplicity, transparency, and tunability.

The strength of the CG approach lies in its ability to **identify qualitative trends**, not to reproduce atomic-scale structure.

## 5. A note on graphite surfaces (based on the Martini 3 graphene model)

In addition to silica, we will also use **graphitic surfaces** as an alternative solid substrate. In this course, graphite is modeled using the **Martini 3 graphene model**, which represents an infinite, periodic graphene sheet and can be used as an effective model for graphitic surfaces.

The graphene model employed here follows the Martini 3 parametrization developed by [Shrestha et al. (2025)](https://doi.org/10.1016/j.surfin.2025.106997), which can reproduces resonable well some structural, elastic, and adsorption trends of graphene within the Martini framework. The model uses a hexagonal lattice of tiny (T) beads and is designed to be used as a **rigid, periodic surface**, making it well suited for interfacial simulations. The model has well-defined mapping rules and should even be able to reproduce overall packing of graphite and even stacking distances.

In contrast to silica, **no parameter exploration is performed for graphite** in this course. The goal is simply to generate a surface of comparable lateral size to the silica slab, which can then be used in solid–ionic liquid interface simulations.

---

### 5.1 Building a periodic graphene sheet

Periodic graphene sheets are generated using `martini3-graphene-periodic.py` (available at  
https://github.com/MoMS-MMSB/Martini3-Graphene). A copy of the script is also provided in `00_templates/`.

A typical command to generate a periodic graphene sheet is:

```bash
python3 martini3-graphene-periodic.py \
  -x 10 \
  -y 10 \
  -o graphene
```

After generation, inspect the box dimensions written in the last line of the generated `.gro` file:

```bash
tail -n 1 graphene.gro
```

The reported *x* and *y* values will be reused when redefining the box height.

---

### 5.2 Defining the interlayer spacing

The graphene builder ensures correct lateral periodicity (*x* and *y*), but the box height (*z*) must be defined manually before stacking layers.

In the Martini 3 graphene model, carbon atoms are represented by **Tiny beads** with a Lennard–Jones size parameter of approximately  
\(\sigma \approx 0.34\) nm. For a standard Lennard–Jones 12–6 potential, the minimum of the interaction potential occurs at:

\[
r_\text{min} = 2^{1/6}\sigma \approx 0.382\ \text{nm}
\]

We use this value as an effective **interlayer spacing** when constructing a graphitic slab.

Redefine the box height using the *x* and *y* values obtained from `tail -n 1`:

```bash
gmx editconf -f graphene.gro -box X Y 0.382 -o graphene_1layer.gro
```

where `X` and `Y` are taken directly from the last line of `graphene.gro`.

---

### 5.3 Stacking graphene layers to obtain graphite

A graphitic slab is constructed by stacking multiple graphene layers along the *z* direction. In this course, we use **five layers**, which is sufficient to represent a graphite substrate at the coarse-grained level.

```bash
gmx genconf -f graphene_1layer.gro -o graphite.gro -nbox 1 1 5
```

This produces a five-layer graphite slab with a total thickness of  
\(5 \times 0.382 \approx 1.9\) nm.


