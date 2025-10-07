# MAKE_KPOINTS_PRIMITIVE

## Overview
Many crystal structures have body-centered orthorhombic (I-centered) cells, while the primitive cell is smaller and may be triclinic. If one naively applies a k-point grid in the primitive cell using VASP’s automatic Monkhorst–Pack/Gamma-centered grids, the reciprocal lattice symmetry can be broken. This can result in warnings such as:
```
Your reciprocal lattice and k-lattice belong to different lattice classes:
The reciprocal lattice is face-centered orthorhombic, whereas your k-lattice is triclinic.
```
The goal of this script is to transform a conventional-cell k-grid into the primitive cell, preserving symmetry as much as possible.\
The script applies the primitive transformation matrix to the conventional k-points. The k-point positions are then explicitly calculated in fractional coordinates of the primitive cell. The script can rescale the grid to approximately match either the conventional-cell directional spacing or the minimum inter-k-point distance, keeping the mesh balanced and maintaining symmetry. Duplicate k-points are collapsed, and their weights are correctly summed.

## Usage
Thte script is designed to be executed directly from command-line and all arguments are provided via CLI. For the most basic functionality, one needs to specify:\
- `poscar`: conventional (body-centered) cell POSCAR file 
- `n1 n2 n3`: grid dimensions of the Gamma-centered/Monkhorst-Pack grid used in the conventional cell
- `--shift` (*optional*): fractional shift of the grid (e.g. `0.5 0.5 0.5` for MP-grids)
- `--normalize-weights` (*optional*): by default, the k-point weights are all equal to 1, One can optionally normalise them to make the sum equal to 1.

The results are written to a file `KPOINTS_explicit` and can directly be used in a VASP calculation. The script also transforms the conventional unit cell to a primitive one and removes any duplicate atoms, and writes this to a file `POSCAR_prim`. It is recommended to use this POSCAR in the calculation.
### Resaling options
The script offers three distinct ways to generate primitive-cell k-grids:
1. **Direct conversion**: Simply transforms the conventional grid to the primitive cell using the conventional to primitive transformation matrix. However, this may result in non-uniform spacings and does not rescale the grid.
2. **Match Conventional Directional Spacing** (`--match-conventional-spacing`): Rescales the primitive k-grid so that the directional spacing along each axis is approximately the same as in the conventional cell.\
example usage:
```
python make_kpoints_primitive.py POSCAR 2 5 5 --match-conventional-spacing
```
3. **Match Minimum Inter-K-Point Distance** (`--match-min-distance`): Rescales the primitive k-grid so that the minimum distance between k-points matches the conventional grid. This should generate a k-grid with a density that roughly matches the one of the conventional cell, using the minimal number of k-points. The directional spacings may differ slightly, but the minimum spacing in 3D is matched. This works by supplying an initial primitive grid (by default this is set to $2\times 2\times 2$, but one can adjust this using the `--init-prim-grid` flag) as a starting point, which is then used to scale the primitive k-point numbers uniformly $n_{\mathrm{prim,new}}=\mathrm{round}(\frac{d_{\mathrm{conv}}}{d_{\mathrm{prim,initial}}})$. This avoids over-dense grids sometimes obtained with the previous method.\
example usage:
```
python make_kpoints_primitive.py POSCAR 2 5 5 --match-min-distance --init-prim-grid 3 3 3
```

### Spacing and convergence checks
The script prints Directional reciprocal spacings and exact minimum inter-k-point distances in both conventional and primitive cells. One can use this information to decide on the optimum grid for the primitive cell (the minimum distance is often more relevant for total energy and force convergece).\
An anisotropy check is also performed, warning if the ratio between the largest and smallest directional spacing exceeds 2.\
***Important note:***\
Even after transformation, reciprocal spacings in the primitive cell may differ from the conventional cell. The user must perform k-point convergence tests independently for each primitive grid before production calculations. One should also check that all symmetry oprations found in the conventional cell are preserved in the primitive cell with the given grid (this information is usually printed in the OUTCAR file, depending on the VASP version one might have to enable more verbose output).