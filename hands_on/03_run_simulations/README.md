# Hands-on 03 — Running CG Molecular Dynamics and Analysis

## Objective

In this hands-on, we run short **coarse-grained (Martini 3)** molecular dynamics simulations of the solid–ionic liquid interfaces built in Hands-on 02 and perform the **minimum analysis** needed to assess interfacial organization.

The central question is:

> Do silica or graphite surfaces induce additional ordering in imidazolium-based ionic liquids (C2, C4, C8)?

---

## Simulation protocol (GROMACS)

We use a simple and consistent protocol for all systems:

1. Energy minimization
2. Short equilibration (NVT, 10 fs time step)
3. Production run (NVT, 20 fs time step)

Interfaces include vacuum along *z*, therefore **NVT simulations are used throughout**.

Position restraints are applied to the solid during minimization and equilibration (silica; and also graphite when applicable).

---

## Running the simulations

### Energy minimization

```bash
gmx grompp -f min.mdp -p topol.top -c start.gro -o min.tpr -maxwarn 1 -r start.gro
gmx mdrun  -nt 8 -v -deffnm min
```

### Equilibration (NVT, 10 fs, Berendsen thermostat)

```bash
gmx grompp -f eq.mdp -p topol.top -c min.gro -o eq.tpr -maxwarn 1 -r start.gro
gmx mdrun  -nt 8 -v -deffnm eq
```

### Production (NVT, 20 fs, v-rescale thermostat)

```bash
gmx grompp -f md.mdp -p topol.top -c eq.gro -o md.tpr -maxwarn 1
gmx mdrun  -nt 8 -v -deffnm md
```

---

## Notes on CG simulations at interfaces

- Martini models allow **larger time steps** than atomistic MD.
- Berendsen thermostat is used only for equilibration.
- v-rescale thermostat is used for production to ensure correct ensemble sampling.
- Periodic boundary conditions are applied in **x and y only**.
- The vacuum region along *z* prevents artificial interactions between periodic images.

---

## Optional: Martini simulations with OpenMM

Martini can also be run efficiently using OpenMM. An example tutorial is available at:

```
https://github.com/maccallumlab/martini_openmm/tree/master/tutorial
```

This is optional for the course; the reference workflow uses GROMACS.

---

## Analysis: interfacial structure

### Index groups

We define index groups to separate solid, cation, anion, and specific bead types within the ionic liquid.

```bash
gmx make_ndx -f min.gro -o index.ndx
```

Typical groups include:

- Solid beads (e.g., SIS, SIC)
- Cation (BIM)
- Anion (BF4)
- Imidazolium head beads (e.g., SI1, SI2)
- Alkyl tail beads (e.g., SI3, SI4, … depending on chain length)

---

### Density profiles along z

Number density profiles along the **z direction** are used to quantify layering near the surface.

Example command (adjust group numbers as needed):

```bash
echo 3 2 6 7 5 | gmx density   -s md.tpr   -f md.xtc   -n index.ndx   -ng 5   -center   -dens number   -d Z   -sl 500
```

Notes:

- `-center` aligns the profiles using the reference solid group.
- `-sl 500` provides sufficient resolution for long boxes.
- Negative values in plots usually correspond to the **z-axis**, not negative densities.

---

## Minimum deliverables

By the end of this hands-on, each group should have:

- A short production trajectory (`md.xtc`, `md.gro`)
- Density profiles along z for:
  - Solid reference
  - Cation head beads
  - Cation tail beads
  - Anion
- A short interpretation (2–3 sentences) addressing:
  - Presence or absence of layering
  - Effect of alkyl chain length (C2 vs C4 vs C8)
  - Differences between silica and graphite

---

## Discussion prompts

- Does the surface induce ordering not present in bulk?
- How does alkyl chain length affect interfacial structure?
- Are silica and graphite interfaces qualitatively different?
- What changes would be required for charged surfaces?
