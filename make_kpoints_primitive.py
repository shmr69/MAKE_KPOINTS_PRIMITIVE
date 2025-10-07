#!/usr/bin/env python3
"""
make_kpoints_primitive.py

Work around the issue of compatible k-grid sampling between orthorhombic body-centred conventional cells and associated triclinic primitive cells.
Generate explicit KPOINTS for the primitive cell from a conventional-cell Gamma-centred grid.
Supports:
  - Optional fractional grid shifts (--shift) e.g., 0.5 0.5 0.5 for Monkhorst–Pack
  - Optional weight normalization (--normalize-weights): Normalize k-point weights so they sum to 1.0
  - Two rescaling options:
        * --match-conventional-spacing: Rescale primitive grid to match conventional directional spacing
        * --match-min-distance: Rescale primitive grid to match min inter-k-point distance of conventional grid
  - Writes POSCAR_prim for the primitive cell
"""

import argparse
import numpy as np
from collections import defaultdict
from itertools import combinations

def read_lattice_from_poscar(filename):
    """
    Read the conventional lattice parameters from a POSCAR file.
    Inputs:
        - filename of the POSCAR file of the conventional body-centred cell
    Outputs:
        - lattice matrix as 3x3 np.ndarray
    """
    with open(filename, "r") as f:
        lines = [l.strip() for l in f.readlines()]
    scale = float(lines[1].split()[0])
    a1 = np.array([float(x) for x in lines[2].split()])
    a2 = np.array([float(x) for x in lines[3].split()])
    a3 = np.array([float(x) for x in lines[4].split()])
    A = np.vstack([a1, a2, a3]) * scale
    return A

def read_poscar_atoms(filename):
    """
    Read the atomic positions from the use-supplied POSCAR file of the conventional body-centred unit cell.
    Inputs:
        - filename of the POSCAR file of the conventional body-centred cell
    Outputs:
        - scale factor : float
        - atom types : list[str]
        - atom multiplicities : list[int]
        - positions : Nx3 np.ndarray
    """
    with open(filename, "r") as f:
        lines = [l.rstrip() for l in f.readlines()]
    scale = float(lines[1])
    atom_types = lines[5].split()
    atom_counts = [int(x) for x in lines[6].split()]
    coord_start = 8
    positions = []
    for count in atom_counts:
        for _ in range(count):
            pos = np.array([float(x) for x in lines[coord_start].split()[:3]])
            positions.append(pos)
            coord_start += 1
    positions = np.array(positions)
    return scale, atom_types, atom_counts, positions

def remove_duplicates_fractional(positions, tol=1e-5):
    """
    Remove duplicate atoms in fractional coordinates considering periodic boundary conditions.
    Inputs:
        positions: Nx3 array of fractional coordinates
        tol: tolerance for minimum-image distance
    Outputs:
        unique_positions: Mx3 array of unique positions
    """
    positions = positions % 1.0  # wrap into [0,1)
    unique_positions = []
    for pos in positions:
        # Compute minimum-image distance to all previously accepted atoms
        is_duplicate = False
        for u in unique_positions:
            delta = np.abs(pos - u)
            delta = np.minimum(delta, 1.0 - delta)  # minimum image
            if np.all(delta < tol):
                is_duplicate = True
                break
        if not is_duplicate:
            unique_positions.append(pos)
    return np.array(unique_positions)


def write_poscar(filename, scale, lattice, atom_types, atom_counts, positions):
    """
    Write the transformed primitive cell to a POSCAR file.
    Inputs:
        - filename : str
        - scale : float
        - lattice matrix : 3x3 array-like
        - atom types : list[str]
        - atom multiplicities : list[int]
        - positions : Nx3 np.ndarray
    Outputs: None
    """
    tol = 1e-5  # tolerance for duplicate detection in fractional coordinates
    # split positions by atom type
    pos_by_type = []
    start = 0
    for count in atom_counts:
        pos_by_type.append(positions[start:start+count])
        start += count

    unique_positions_by_type = []
    new_counts = []
    for pos_type in pos_by_type:
        unique_pos = remove_duplicates_fractional(pos_type, tol)
        unique_positions_by_type.append(unique_pos)
        new_counts.append(len(unique_pos))
    
    # concatenate unique positions in order of atom types
    final_positions = np.vstack(unique_positions_by_type)
    # write POSCAR
    with open(filename, "w") as f:
        f.write(f"Primitive cell generated from {poscar_file} (duplicates removed)\n")
        f.write(f"{scale:.16f}\n")
        for vec in lattice:
            f.write("  " + "  ".join(f"{x:.16f}" for x in vec) + "\n")
        f.write("  " + "  ".join(atom_types) + "\n")
        f.write("  " + "  ".join(str(x) for x in new_counts) + "\n")
        f.write("Direct\n")
        for pos in final_positions:
            f.write("  " + "  ".join(f"{x:.16f}" for x in pos) + "\n")

def reciprocal_lattice(A):
    """
    Calculate the recprocal lattice given a real-space lattice.
    Inputs:
        - real-space lattice matrix : 3x3 array-like
    Outputs:
        - reciprocal lattice matrix : 3x3 array-like
    """
    vol = np.dot(A[0], np.cross(A[1], A[2]))
    B = 2*np.pi*np.array([
        np.cross(A[1], A[2]),
        np.cross(A[2], A[0]),
        np.cross(A[0], A[1])
    ]).T / vol
    return B

def min_k_distance(kcart):
    """
    Calculates the minimum distance between k-points.
    Inputs:
        - list of k-point coordinates in cartesian coordinates : Nx3 array-like 
    Oututs:
        - minumum distance for the given set : float
    """
    min_dist = np.inf
    for p, q in combinations(kcart, 2):
        d = np.linalg.norm(p - q)
        if d < min_dist:
            min_dist = d
    return min_dist

def generate_kgrid(n1, n2, n3, shift):
    """
    Generates a regular k-grid in fractional coordinates given number of spacings along each reciprocal lattice vector.
    Inputs:
        - number of subdivisions n1, n2, n3 along each reciprocal lattice vector : float
        - shift (for grid offsets, e.g. (1/2,1/2,1/2) for MP grids) : float
    Outputs:
        - regular k-point grid : np.ndarray
    """
    kgrid = []
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                frac = ((np.array([i, j, k]) + shift) / np.array([n1, n2, n3])) % 1.0
                kgrid.append(frac)
    return np.array(kgrid)



# CLI argument parsing
parser = argparse.ArgumentParser(
    description="Generate explicit KPOINTS for primitive cell from conventional Gamma-centred grid."
)
parser.add_argument("poscar", type=str, help="POSCAR file of the conventional cell")
parser.add_argument("n1", type=int, help="Number of k-points along a")
parser.add_argument("n2", type=int, help="Number of k-points along b")
parser.add_argument("n3", type=int, help="Number of k-points along c")
parser.add_argument("--shift", nargs=3, type=float, default=[0.0, 0.0, 0.0],
                    help="Fractional shift (e.g., 0.5 0.5 0.5 for Monkhorst–Pack)")
parser.add_argument("--normalize-weights", action="store_true",
                    help="Normalize k-point weights so they sum to 1.0")
group = parser.add_mutually_exclusive_group()
group.add_argument("--match-conventional-spacing", action="store_true",
                   help="Rescale primitive grid to match conventional directional spacing")
group.add_argument("--match-min-distance", action="store_true",
                   help="Rescale primitive grid to match min inter-k-point distance of conventional grid")

parser.add_argument("--init-prim-grid", nargs=3, type=int, default=[2,2,2],
                    help="Initial primitive grid (used only with --match-min-distance), default: 2 2 2")

args = parser.parse_args()

n1, n2, n3 = args.n1, args.n2, args.n3
shift = np.array(args.shift)
poscar_file = args.poscar
normalize = args.normalize_weights
match_spacing = args.match_conventional_spacing
match_min_distance = args.match_min_distance
init_prim_grid = np.array(args.init_prim_grid)


# read lattice and atomic positions
A_conv = read_lattice_from_poscar(poscar_file)
scale, atom_types, atom_counts, positions = read_poscar_atoms(poscar_file)

# transformation matrix from body-centred conventional to triclinic primitive cell
T = np.array([
    [-0.5,  0.5,  0.5],
    [ 0.5, -0.5,  0.5],
    [ 0.5,  0.5, -0.5]
])

# reciprocal lattices
B_conv = reciprocal_lattice(A_conv)
B_prim = B_conv @ T

# generate conventional k-grid (if rescaling options not requested)
k_conv = generate_kgrid(n1, n2, n3, shift)
k_cart_conv = (B_conv @ k_conv.T).T
min_dist_conv = min_k_distance(k_cart_conv)


# transform to primitive fractional coordinates
k_prim = (T.T @ k_conv.T).T % 1.0
k_cart_prim = (B_prim @ k_prim.T).T
min_dist_prim = min_k_distance(k_cart_prim)


# collapse duplicates and assign weights
tol_round = 12
weights = defaultdict(int)
for kp in k_prim:
    key = tuple(np.round(kp, tol_round))
    weights[key] += 1

unique_k = np.array(list(weights.keys()))
multiplicities = np.array(list(weights.values()), dtype=float)
if normalize:
    multiplicities /= multiplicities.sum()

# calculate directional spacings along each reciprocal lattice vector
spacing_conv = np.linalg.norm(B_conv, axis=0) / np.array([n1, n2, n3])
spacing_prim = np.linalg.norm(B_prim, axis=0) / np.array([n1, n2, n3])

print("\nReciprocal-lattice spacing comparison (Å⁻¹):")
print("  Directional (conv):", "  ".join(f"{x:.4f}" for x in spacing_conv))
print("  Directional (prim, original):", "  ".join(f"{x:.4f}" for x in spacing_prim))
print("  Ratio (prim/conv): ", "  ".join(f"{(p/c):.3f}" for p,c in zip(spacing_prim, spacing_conv)))


# rescaling options
if match_spacing:
    # --match-conventional-spacing
    lens_prim = np.linalg.norm(B_prim, axis=0)
    factor = min_dist_conv / min_dist_prim
    n_prim_new = np.round(np.array([n1, n2, n3]) * factor * (lens_prim / lens_prim.mean())).astype(int)
    n_prim_new = np.maximum(n_prim_new, 1)
    print(f"Rescaling primitive grid (match-conventional-spacing): old n={n1,n2,n3}, new n={tuple(map(int,n_prim_new))}")
    k_prim = generate_kgrid(n_prim_new[0], n_prim_new[1], n_prim_new[2], shift)
    k_cart_prim = (B_prim @ k_prim.T).T
    min_dist_prim = min_k_distance(k_cart_prim)
    # directional spacings for the rescaled grid
    spacing_prim_rescaled = np.linalg.norm(B_prim, axis=0) / np.array([n_prim_new[0], n_prim_new[1], n_prim_new[2]])
    print("  Directional reciprocal spacing of rescaled primitive grid (Å⁻¹):", 
          "  ".join(f"{x:.4f}" for x in spacing_prim_rescaled))
    # anisotropy check
    max_s = spacing_prim_rescaled.max()
    min_s = spacing_prim_rescaled.min()
    ratio = max_s / min_s if min_s > 0 else np.inf
    if ratio > 2.0:
        print(f"\n WARNING: Primitive reciprocal grid is anisotropic (max/min spacing ratio = {ratio:.2f}).")
        print("   Consider adjusting the grid or transformation if this is undesirable.\n")
elif match_min_distance:
    # --match-min-distance (density-matched)
    # start with small initial grid along primitive axes
    print(f"Using initial primitive grid: {tuple(init_prim_grid)}")
    k_test = generate_kgrid(init_prim_grid[0], init_prim_grid[1], init_prim_grid[2], shift)
    k_cart_test = (B_prim @ k_test.T).T
    d_prim_init = min_k_distance(k_cart_test)
    factor = d_prim_init / min_dist_conv
    n_prim_new = np.round(init_prim_grid / factor).astype(int)
    n_prim_new = np.maximum(n_prim_new, 1)
    print(f"Rescaling primitive grid (match-min-distance): old n={n1,n2,n3}, new n={tuple(map(int,n_prim_new))}")
    k_prim = generate_kgrid(n_prim_new[0], n_prim_new[1], n_prim_new[2], shift)
    k_cart_prim = (B_prim @ k_prim.T).T
    min_dist_prim = min_k_distance(k_cart_prim)
    # directional spacings for the rescaled grid
    spacing_prim_rescaled = np.linalg.norm(B_prim, axis=0) / np.array([n_prim_new[0], n_prim_new[1], n_prim_new[2]])
    print("  Directional reciprocal spacing of rescaled primitive grid (Å⁻¹):", 
          "  ".join(f"{x:.4f}" for x in spacing_prim_rescaled))
    # anisotropy check
    max_s = spacing_prim_rescaled.max()
    min_s = spacing_prim_rescaled.min()
    ratio = max_s / min_s if min_s > 0 else np.inf
    if ratio > 2.0:
        print(f"\n WARNING: Primitive reciprocal grid is anisotropic (max/min spacing ratio = {ratio:.2f}).")
        print("   Consider adjusting the grid or transformation if this is undesirable.\n")

print(f"  Exact minimum inter-k-point distance (conventional grid): {min_dist_conv:.4f} Å⁻¹")
print(f"  Exact minimum inter-k-point distance (primitive grid):   {min_dist_prim:.4f} Å⁻¹\n")





# write KPOINTS_explicit
outfile = "KPOINTS_explicit"
lines = []
lines.append(f"Explicit primitive grid (shift {shift.tolist()})\n")
lines.append(f"{len(unique_k)}\n")
lines.append("Reciprocal\n")
for kvec, wt in zip(unique_k, multiplicities):
    if normalize:
        lines.append("{: .12f} {: .12f} {: .12f} {: .6f}\n".format(kvec[0], kvec[1], kvec[2], wt))
    else:
        lines.append("{: .12f} {: .12f} {: .12f} {:d}\n".format(kvec[0], kvec[1], kvec[2], int(wt)))

with open(outfile, "w") as f:
    f.writelines(lines)
print(f"Wrote {outfile} with {len(unique_k)} unique k-points (weights sum = {multiplicities.sum():.4f})")


# compare total number of k-points
total_conv_kpoints = len(k_conv)
total_prim_kpoints = len(k_prim)
print(f"\nTotal number of k-points:")
print(f"  Conventional grid: {total_conv_kpoints}")
print(f"  Primitive grid:    {total_prim_kpoints}")



# write primitive POSCAR
A_prim = T @ A_conv
positions_cart = (A_conv.T @ positions.T).T
positions_prim = np.linalg.solve(A_prim.T, positions_cart.T).T % 1.0
write_poscar("POSCAR_prim", 1.0, A_prim, atom_types, atom_counts, positions_prim)
print("\nWrote primitive cell structure to POSCAR_prim")
