#!/usr/bin/env python3
import argparse
import math
import random
import sys


def _choose_fraction(indices, frac, rng):
    """Return a set of selected indices (subset) of length round(frac*n)."""
    n = len(indices)
    if n == 0:
        return set()
    if frac <= 0.0:
        return set()
    if frac >= 1.0:
        return set(indices)
    k = int(round(frac * n))
    k = max(0, min(n, k))
    return set(rng.sample(indices, k))


def _name_triplet(prefix: str):
    """
    Return a 3-name cycle suitable for PDB/GRO (<=4 chars).
    Examples: prefix='S' -> S1,S2,S3 ; prefix='TM' -> TM1,TM2,TM3
    """
    prefix = (prefix or "X").strip()
    if len(prefix) > 3:
        prefix = prefix[:3]
    return [f"{prefix}1", f"{prefix}2", f"{prefix}3"]


def build_slab(
    slab_type,
    target_x,
    target_y,
    target_z,
    core_type_cli,
    surface_type_cli,
    top_type_cli,
    bottom_type_cli,
    basename,
    pr=None,
    a_user=None,
    dens=0.9,
    use_bonds=False,
    ef=10000.0,
    rcut=0.65,
    pbc_bonds="none",  # none | xy | xyz
    core_charge=0.0,
    surface_charge=0.0,
    top_charge=None,
    bottom_charge=None,
    # Random surface modification (two-state surface)
    top_mod_type=None,
    bottom_mod_type=None,
    top_mod_charge=None,
    bottom_mod_charge=None,
    top_mod_frac=None,
    bottom_mod_frac=None,
    seed=1,
    # Naming for visualization
    surf_name_prefix="S",      # base surface names: S1,S2,S3
    top2_name_prefix="T",      # modified TOP names: T1,T2,T3
    bottom2_name_prefix="B",   # modified BOTTOM names: B1,B2,B3
):
    """
    Two-state surface option (base + modified subset) with distinct atom names:
      - Base surface beads use names like S1/S2/S3 (configurable via -surf_name_prefix).
      - Modified TOP subset uses T1/T2/T3 (configurable via -top2_name_prefix).
      - Modified BOTTOM subset uses B1/B2/B3 (configurable via -bottom2_name_prefix).

    Fractions:
      - If a mod type is provided but fraction is omitted, defaults to 0.5 (equal split).
      - If no mod type is provided, fraction is ignored.
    """
    rng = random.Random(seed)

    # --- Lattice spacing ---
    a_ref = 0.53
    a0 = a_user if a_user is not None else a_ref * dens

    # Default bead types
    if slab_type == "crystalline":
        default_core_type = "N2"
        default_surface_type = "N2"
    elif slab_type == "amorphous":
        default_core_type = "N5"
        default_surface_type = "N5"
    else:
        raise ValueError("Unknown slab type")

    core_type = core_type_cli or default_core_type
    surface_type = surface_type_cli or default_surface_type

    # Top/bottom overrides
    top_type = top_type_cli or surface_type
    bottom_type = bottom_type_cli or surface_type

    # Charges
    top_q = surface_charge if top_charge is None else top_charge
    bottom_q = surface_charge if bottom_charge is None else bottom_charge

    # If mod type provided but fraction omitted -> equal split
    if top_mod_type is not None and top_mod_frac is None:
        top_mod_frac = 0.5
    if bottom_mod_type is not None and bottom_mod_frac is None:
        bottom_mod_frac = 0.5

    # If mod charge omitted, fall back to base top/bottom charges
    if top_mod_type is not None and top_mod_charge is None:
        top_mod_charge = top_q
    if bottom_mod_type is not None and bottom_mod_charge is None:
        bottom_mod_charge = bottom_q

    core_resname = "SIC"
    surf_resname = "SIS"

    # --- Lattice ---
    dx = a0
    dy = a0 * math.sqrt(3.0) / 2.0
    dz = a0 * math.sqrt(2.0 / 3.0)

    nx = max(1, int(round(target_x / dx)))
    ny = max(1, int(round(target_y / dy)))
    n_layers = max(1, int(round(target_z / dz)))

    box_x_nm = nx * dx
    box_y_nm = ny * dy
    slab_thickness_nm = n_layers * dz
    box_z_nm = max(target_z, slab_thickness_nm)

    atoms = []
    atom_id = 1
    for layer in range(n_layers):
        z = layer * dz
        for j in range(ny):
            for i in range(nx):
                x = i * dx + (0.5 * a0 if (j % 2 == 1) else 0.0)
                y = j * dy
                atoms.append({"id": atom_id, "x": x, "y": y, "z": z, "layer": layer})
                atom_id += 1

    # Tag regions
    top_ids = []
    bottom_ids = []
    for a in atoms:
        if a["layer"] == 0:
            a["region"] = "bottom"
            bottom_ids.append(a["id"])
        elif a["layer"] == n_layers - 1:
            a["region"] = "top"
            top_ids.append(a["id"])
        else:
            a["region"] = "core"

    # Select modified subsets
    top_mod_set = set()
    bottom_mod_set = set()
    if top_mod_type is not None:
        top_mod_set = _choose_fraction(top_ids, float(top_mod_frac), rng)
    if bottom_mod_type is not None:
        bottom_mod_set = _choose_fraction(bottom_ids, float(bottom_mod_frac), rng)

    # Center slab
    zs = [a["z"] for a in atoms]
    shift_z = 0.5 * box_z_nm - 0.5 * (min(zs) + max(zs))
    for a in atoms:
        a["z"] += shift_z

    # Names + types + charges
    core_names = ["C1", "C2", "C3"]
    surf_names_base = _name_triplet(surf_name_prefix)
    surf_names_top2 = _name_triplet(top2_name_prefix)
    surf_names_bottom2 = _name_triplet(bottom2_name_prefix)

    ci = 0
    ti_base = ti_mod = 0
    bi_base = bi_mod = 0

    for a in atoms:
        a["resid"] = 1
        if a["region"] == "core":
            a["resname"] = core_resname
            a["name"] = core_names[ci % 3]
            a["bead_type"] = core_type
            a["charge"] = float(core_charge)
            ci += 1

        elif a["region"] == "top":
            a["resname"] = surf_resname
            if a["id"] in top_mod_set:
                # Modified TOP bead
                a["name"] = surf_names_top2[ti_mod % 3]
                a["bead_type"] = top_mod_type
                a["charge"] = float(top_mod_charge)
                ti_mod += 1
            else:
                # Base TOP bead
                a["name"] = surf_names_base[ti_base % 3]
                a["bead_type"] = top_type
                a["charge"] = float(top_q)
                ti_base += 1

        else:  # bottom
            a["resname"] = surf_resname
            if a["id"] in bottom_mod_set:
                # Modified BOTTOM bead
                a["name"] = surf_names_bottom2[bi_mod % 3]
                a["bead_type"] = bottom_mod_type
                a["charge"] = float(bottom_mod_charge)
                bi_mod += 1
            else:
                # Base BOTTOM bead
                a["name"] = surf_names_base[bi_base % 3]
                a["bead_type"] = bottom_type
                a["charge"] = float(bottom_q)
                bi_base += 1

    # --- Neighbour detection (distance based) ---
    positions = [(a["x"], a["y"], a["z"]) for a in atoms]
    Lx, Ly, Lz = box_x_nm, box_y_nm, box_z_nm
    rcut2 = rcut * rcut

    def min_image(dxv, dyv, dzv):
        if pbc_bonds in ("xy", "xyz"):
            dxv -= Lx * round(dxv / Lx)
            dyv -= Ly * round(dyv / Ly)
        if pbc_bonds == "xyz":
            dzv -= Lz * round(dzv / Lz)
        return dxv, dyv, dzv

    neighbour_pairs = []
    for i in range(len(atoms)):
        xi, yi, zi = positions[i]
        for j in range(i + 1, len(atoms)):
            xj, yj, zj = positions[j]
            dxv = xj - xi
            dyv = yj - yi
            dzv = zj - zi
            dxv, dyv, dzv = min_image(dxv, dyv, dzv)
            r2 = dxv * dxv + dyv * dyv + dzv * dzv
            if r2 <= rcut2:
                neighbour_pairs.append((i + 1, j + 1, math.sqrt(r2)))

    # --- Write PDB (Å) ---
    with open(f"{basename}.pdb", "w") as f:
        f.write(
            "CRYST1{:9.3f}{:9.3f}{:9.3f} 90.00 90.00 90.00 P 1           1\n".format(
                box_x_nm * 10, box_y_nm * 10, box_z_nm * 10
            )
        )
        for a in atoms:
            f.write(
                "ATOM  {:5d} {:>4s} {:>3s} A{:4d}    {:8.3f}{:8.3f}{:8.3f}\n".format(
                    a["id"], a["name"], a["resname"], a["resid"],
                    a["x"] * 10, a["y"] * 10, a["z"] * 10
                )
            )
        if use_bonds:
            for i, j, _ in neighbour_pairs:
                f.write(f"CONECT{i:5d}{j:5d}\n")
        f.write("END\n")

    # --- Write GRO (nm) ---
    with open(f"{basename}.gro", "w") as g:
        g.write("Silica slab\n")
        g.write(f"{len(atoms)}\n")
        for a in atoms:
            g.write(
                "{:5d}{:<5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(
                    a["resid"], a["resname"], a["name"], a["id"], a["x"], a["y"], a["z"]
                )
            )
        g.write("{:10.5f}{:10.5f}{:10.5f}\n".format(box_x_nm, box_y_nm, box_z_nm))

    # --- Write ITP ---
    with open(f"{basename}.itp", "w") as f:
        f.write("; Orthorhombic silica slab\n")
        f.write("[ moleculetype ]\n")
        f.write("SILICA_SLAB  1\n\n")

        f.write("[ atoms ]\n")
        f.write("; nr  type  resnr  resid  atom  cgnr  charge  mass\n")
        for a in atoms:
            f.write(
                "{:5d} {:>6s} {:4d} {:>4s} {:>4s} {:4d} {:10.6f} {:7.1f}\n".format(
                    a["id"], a["bead_type"], a["resid"], a["resname"], a["name"],
                    a["id"], float(a["charge"]), 72.0
                )
            )

        if use_bonds:
            f.write("\n[ bonds ]\n")
            f.write("; i j funct b0(nm) kb(kJ/mol/nm^2)\n")
            for i, j, r in neighbour_pairs:
                f.write("{:5d}{:6d}{:6d}{:11.5f}{:11.1f}\n".format(i, j, 1, r, ef))

        if pr is not None:
            f.write("\n#ifdef POSRES\n")
            f.write("[ position_restraints ]\n")
            f.write("; ai  funct  fcx    fcy    fcz   (kJ/mol/nm^2)\n")
            for a in atoms:
                f.write("{:5d} 1 {:7.1f} {:7.1f} {:7.1f}\n".format(a["id"], pr, pr, pr))
            f.write("#endif\n")

    # Print a short summary
    if top_mod_type is not None:
        print(f"Top surface: {len(top_mod_set)}/{len(top_ids)} beads modified "
              f"({(len(top_mod_set)/max(1,len(top_ids))):.3f}) -> {top_mod_type} names {surf_names_top2[0]}/{surf_names_top2[1]}/{surf_names_top2[2]}")
    if bottom_mod_type is not None:
        print(f"Bottom surface: {len(bottom_mod_set)}/{len(bottom_ids)} beads modified "
              f"({(len(bottom_mod_set)/max(1,len(bottom_ids))):.3f}) -> {bottom_mod_type} names {surf_names_bottom2[0]}/{surf_names_bottom2[1]}/{surf_names_bottom2[2]}")


def main():
    p = argparse.ArgumentParser(
        description="Build a silica slab (PDB+GRO+ITP) with optional bonds, optional POSRES, and optional random surface modification with distinct atom names for visualization."
    )
    p.add_argument("-type", required=True, choices=["crystalline", "amorphous"])
    p.add_argument("-x", type=float, required=True, help="Target size along x (nm).")
    p.add_argument("-y", type=float, required=True, help="Target size along y (nm).")
    p.add_argument("-z", type=float, required=True, help="Target thickness along z (nm).")

    p.add_argument("-core", help="Core bead type override (default depends on -type).")
    p.add_argument("-surface", help="Surface bead type override (both top+bottom unless overridden).")
    p.add_argument("-top", dest="top_type", help="Top surface bead type override.")
    p.add_argument("-bottom", dest="bottom_type", help="Bottom surface bead type override.")

    p.add_argument("-qcore", type=float, default=0.0, help="Core bead partial charge (default 0.0).")
    p.add_argument("-qsurf", type=float, default=0.0, help="Surface partial charge for both top+bottom (default 0.0).")
    p.add_argument("-qtop", type=float, default=None, help="Top surface partial charge override.")
    p.add_argument("-qbottom", type=float, default=None, help="Bottom surface partial charge override.")

    # Random modifications (two-state surface)
    p.add_argument("-top2", dest="top_mod_type", default=None,
                   help="Second bead type on TOP surface (randomly assigned to a fraction).")
    p.add_argument("-bottom2", dest="bottom_mod_type", default=None,
                   help="Second bead type on BOTTOM surface (randomly assigned to a fraction).")
    p.add_argument("-qtop2", dest="top_mod_charge", type=float, default=None,
                   help="Charge for TOP second bead type (default: same as top base charge).")
    p.add_argument("-qbottom2", dest="bottom_mod_charge", type=float, default=None,
                   help="Charge for BOTTOM second bead type (default: same as bottom base charge).")
    p.add_argument("-ftop2", dest="top_mod_frac", type=float, default=None,
                   help="Fraction of TOP surface beads assigned to -top2 (0..1). If omitted, defaults to 0.5 when -top2 is set.")
    p.add_argument("-fbottom2", dest="bottom_mod_frac", type=float, default=None,
                   help="Fraction of BOTTOM surface beads assigned to -bottom2 (0..1). If omitted, defaults to 0.5 when -bottom2 is set.")
    p.add_argument("-seed", type=int, default=1, help="Random seed for surface modification (default 1).")

    # Naming options (visualization)
    p.add_argument("-surf_name_prefix", default="S",
                   help="Atom name prefix for base surface beads (default 'S' -> S1/S2/S3).")
    p.add_argument("-top2_name_prefix", default="T",
                   help="Atom name prefix for modified TOP beads (default 'T' -> T1/T2/T3).")
    p.add_argument("-bottom2_name_prefix", default="B",
                   help="Atom name prefix for modified BOTTOM beads (default 'B' -> B1/B2/B3).")

    p.add_argument("-pr", type=float, help="Position restraint k (kJ/mol/nm^2). Emits POSRES-guarded block.")
    p.add_argument("-a", type=float, help="Absolute lattice spacing a0 (nm). Overrides -dens if given.")
    p.add_argument("-dens", type=float, default=0.9, help="Density multiplier on reference a_ref=0.53 nm.")
    p.add_argument("-bonds", action="store_true", help="Write [ bonds ] for neighbours within -rcut.")
    p.add_argument("-ef", type=float, default=10000.0, help="Bond force constant (kJ/mol/nm^2).")
    p.add_argument("-rcut", type=float, default=0.65, help="Neighbour cutoff for bonds (nm). Default 0.65 (=6.5 Å).")
    p.add_argument("-pbc_bonds", choices=["none", "xy", "xyz"], default="none",
                   help="Apply minimum-image when building bonds: none|xy|xyz (default none).")
    p.add_argument("-o", default="silica_slab", help="Basename for output files.")

    args = p.parse_args()

    if args.x <= 0 or args.y <= 0 or args.z <= 0:
        print("Error: -x, -y, -z must be > 0.", file=sys.stderr)
        sys.exit(1)

    if args.top_mod_frac is not None and not (0.0 <= args.top_mod_frac <= 1.0):
        print("Error: -ftop2 must be in [0, 1].", file=sys.stderr)
        sys.exit(1)
    if args.bottom_mod_frac is not None and not (0.0 <= args.bottom_mod_frac <= 1.0):
        print("Error: -fbottom2 must be in [0, 1].", file=sys.stderr)
        sys.exit(1)

    build_slab(
        slab_type=args.type,
        target_x=args.x,
        target_y=args.y,
        target_z=args.z,
        core_type_cli=args.core,
        surface_type_cli=args.surface,
        top_type_cli=args.top_type,
        bottom_type_cli=args.bottom_type,
        basename=args.o,
        pr=args.pr,
        a_user=args.a,
        dens=args.dens,
        use_bonds=args.bonds,
        ef=args.ef,
        rcut=args.rcut,
        pbc_bonds=args.pbc_bonds,
        core_charge=args.qcore,
        surface_charge=args.qsurf,
        top_charge=args.qtop,
        bottom_charge=args.qbottom,
        top_mod_type=args.top_mod_type,
        bottom_mod_type=args.bottom_mod_type,
        top_mod_charge=args.top_mod_charge,
        bottom_mod_charge=args.bottom_mod_charge,
        top_mod_frac=args.top_mod_frac,
        bottom_mod_frac=args.bottom_mod_frac,
        seed=args.seed,
        surf_name_prefix=args.surf_name_prefix,
        top2_name_prefix=args.top2_name_prefix,
        bottom2_name_prefix=args.bottom2_name_prefix,
    )


if __name__ == "__main__":
    main()
