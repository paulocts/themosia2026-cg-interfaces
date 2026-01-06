# Hands-on 01 — Building Solid Surfaces (Silica first)

## 1. Objective

In this hands-on session you will build a **Martini 3 coarse-grained (CG) silica slab** that will later be used to assemble **solid–ionic liquid interfaces**.

By the end of this hands-on you should be able to:

- Generate a CG silica slab with a defined size and thickness.
- Create a surface model with **different bead types in the core and at the surface**.
- Control basic “surface knobs” that influence interfacial behavior:
  - **bead type** (surface chemistry)
  - **bead spacing / surface density** (structure)
  - **rigidity** (bonding / restraints)
- Produce the files needed for the next hands-ons (coordinates + topology).

> Note: In this session we focus on **silica only** and explore parameters.  
> Graphene will be introduced later as an additional surface; it will be a quick build to obtain the required files (less exploration).

---

## 2. Background: CG surfaces in Martini

### 2.1 CG surfaces are *models*, not atomistic “mappings”

In coarse-grained simulations, a solid surface is typically represented as a **mesoscopic model** designed to reproduce key behaviors (e.g., wetting, layering, adsorption), rather than a one-to-one mapping of every atom. This means that several choices are **part of the model definition**, not “details”:

- the **surface bead density** (spacing between beads)
- the **surface bead types / interactions** (effective chemistry)
- whether the surface is treated as **rigid, restrained, or flexible**

These choices can strongly affect interfacial structure and dynamics, and they are expected to be tested and justified depending on the target application. :contentReference[oaicite:0]{index=0}

### 2.2 Key “surface knobs” we will use in this course

In this course, we will focus on three practical knobs that you can tune and immediately interpret:

1. **Surface bead type (chemistry)**
   - We assign different bead types to the **surface** and the **core** to represent different effective interactions with the liquid (following the modeling logic used for silica surfaces in the reference work). :contentReference[oaicite:1]{index=1}

2. **Bead spacing / surface density (structure)**
   - The distance between surface beads controls the effective roughness and interaction density seen by the liquid, which can change layering and organization near the interface. :contentReference[oaicite:2]{index=2}

3. **Rigidity (mechanics)**
   - A surface can be made rigid (frozen / strongly restrained) or mechanically coherent (bonded network). This affects how the surface responds to the liquid and can influence interfacial ordering. :contentReference[oaicite:3]{index=3}

### 2.3 What we will measure later

Once the solid–ionic liquid systems are assembled (later hands-ons), we will quantify surface-induced effects using:

- **partial density profiles** along the surface normal (layering, enrichment/depletion)
- optional: RDFs, diffusion, and orientational order parameters

For now, the goal is to build a clean, reproducible silica slab model with controllable parameters.

