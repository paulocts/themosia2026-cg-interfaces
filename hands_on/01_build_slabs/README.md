# Hands-on 01 — Building Solid Surfaces

## 1. Objective

In this hands-on session you will build a **Martini 3 coarse-grained (CG) silica slab** that will later be used to assemble **solid–ionic liquid interfaces**. The material is based on the work of  [Cambiaso et al 2025](https://doi.org/10.1016/j.surfin.2025.106997).

By the end of this hands-on you should be able to:

- Generate a CG silica slab with a defined size and thickness.
- Create a surface model with **different bead types in the core and at the surface**.
- Control basic “surface properties” that influence interfacial behavior:
  - **bead type** (surface chemistry)
  - **bead spacing / surface density** (structure)
  - **rigidity** (bonding / restraints)
- Produce the files needed for the next hands-ons (coordinates + topology).

---

## 2. Background: CG surfaces in Martini

In coarse-grained simulations, a solid surface can be represented as a **mesoscopic model** designed to reproduce key behaviors (e.g., wetting, layering, adsorption), rather than a one-to-one mapping of every atom. This means that several choices are **part of the model definition**, not “details”:

the **surface bead density**, which can be defined based on the spacing between beads.
- the **surface bead types **, which can define the preferential (effective chemistry)even some charge.
  whether the surface is treated as **rigid, restrained, or flexible**

These choices can strongly affect interfacial structure and dynamics, and they are expected to be tested and justified depending on the target application.


In this course, we will focus on three practical knobs that you can tune and immediately interpret:

1. **Surface bead type (chemistry)**
   - We assign different bead types to the **surface** and the **core** to represent different effective interactions with the liquid (following the modeling logic used for silica surfaces in the reference work). 

2. **Bead spacing / surface density (structure)**
   - The distance between surface beads controls the effective roughness and interaction density seen by the liquid, which can change layering and organization near the interface. :contentReference[oaicite:2]{index=2}

Once the solid–ionic liquid systems are assembled (later hands-ons), we will quantify surface-induced effects using:

- **partial density profiles** along the surface normal (layering, enrichment/depletion)
- optional: RDFs, diffusion, and orientational order parameters

For now, the goal is to build a clean, reproducible silica slab model with controllable parameters.

---


