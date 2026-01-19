#!/usr/bin/env bash
# Reorder atoms in a .gro file by residue-name patterns (extended regex).
#
# Usage:
#   ./reordering.sh input.gro "SOLID_REGEX" "CATION_REGEX" "ANION_REGEX" [output.gro]
#
# Examples:
#   # Silica
#   ./reordering.sh box.gro "1SI" "1BIM" "2BF4" start.gro
#
#   # Graphite (5 layers: 1GRA..5GRA)
#   ./reordering.sh box.gro "[1-5]GRA" "1EIM" "2BF4" start.gro

set -euo pipefail

if [[ $# -lt 4 ]]; then
  echo "Usage: $0 <input.gro> <solid_regex> <cation_regex> <anion_regex> [output.gro]" >&2
  exit 1
fi

in_gro="$1"
solid="$2"
cation="$3"
anion="$4"
out_gro="${5:-start.gro}"

title=$(head -n 1 "$in_gro")
box=$(tail -n 1 "$in_gro")

# Keep only atom lines (skip title + natoms; drop last box line)
tail -n +3 "$in_gro" | head -n -1 > atoms_only.tmp

# Reorder using regex patterns anchored at start of line (after optional spaces)
{
  grep -E "^[[:space:]]*${solid}" atoms_only.tmp || true
  grep -E "^[[:space:]]*${cation}" atoms_only.tmp || true
  grep -E "^[[:space:]]*${anion}" atoms_only.tmp || true
} > atoms_reordered.tmp

natoms=$(wc -l < atoms_reordered.tmp | tr -d '[:space:]')

{
  echo "$title"
  printf "%5d\n" "$natoms"
  cat atoms_reordered.tmp
  echo "$box"
} > "$out_gro"

rm -f atoms_only.tmp atoms_reordered.tmp

echo "Wrote: $out_gro (natoms=$natoms)"
