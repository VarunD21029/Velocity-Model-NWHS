# Fault damage zones in the NW Himalaya velocity model

This folder adds low-velocity fault damage zones to the tomography model in
`../velocity_model_file.zip`. The result is
`../velocity_model_file_faults.zip` (`velocity_model_file_faults/tomography_model.xyz`).
It is a drop-in replacement for the original file: copy it to `DATA/tomo_files/` and keep
`MODEL = tomo` in `Par_file`.

## What was changed

* **Faults included** (traces and geometry in `faults_nwh.json`):

  | Fault | Trace source | Dip / down-dip geometry |
  |---|---|---|
  | MFT → MHT | GEM Global Active Faults / HimaTibetMap; extended NW by hand towards Jammu | 30° ramp to 5 km, 6° décollement to 12 km, 15° mid-crustal ramp to 20 km |
  | Jwalamukhi Thrust (JMT) | approximate, hand-placed | 35° → 25° (listric), soles into MHT |
  | Main Boundary Thrust (MBT) | approximate, hand-placed | 45° → 35° → 25° (listric), soles into MHT |
  | Main Central Thrust (MCT) | approximate, hand-placed | 35° → 25°, soles into MHT |
  | South Tibetan Detachment (STD) | approximate, hand-placed | 25° N, to 15 km |
  | Kaurik-Chango Fault | GEM / HimaTibetMap | 60° W normal, to 15 km |
  | Karakoram Fault | GEM / HimaTibetMap | 80° NE, to 25 km |

  The JMT, MBT, MCT and STD traces were placed by hand from published regional maps,
  so they are only accurate to a few km. Replace them with your own digitised traces
  if you have them.
* **Vs reduction**: 10 % at the surface trace. The reduction follows the fault plane
  down-dip and is 10 % down to 5 km, easing linearly to 5 % at 20 km
  (`vs_reduction_vs_depth`). For a uniform 10 % along the whole plane, set
  `"vs_reduction_vs_depth": [[0, 10]]`.
* **Vp and density** are unchanged by default (`vp_reduction_ratio`,
  `rho_reduction_ratio` = 0), so Vp/Vs rises inside the damage zones.
* **Zone width**: physical half-width 1.5 km normal to the fault, plus a cosine taper of
  0.75 km. The tomography grid is coarse (≈4.1 × 4.7 km horizontally, 0.5 km vertically).
  So steep faults are widened to at least half a horizontal grid step, so they show up
  as continuous zones instead of missing grid nodes.
* Where zones overlap, the **largest** reduction is used. They never stack, so no node is
  reduced by more than 10 %.
* Only the Vs column and the Vs min/max in the header change. Coordinates, record order
  and number format are kept exactly as in the original file.

About 6 % of the grid nodes are affected. QC plots are in `figures/`: map slices at
0/5/10/15 km and three S–N cross-sections.

## Re-running / changing the model

```bash
python3 add_fault_zones.py \
    --input ../velocity_model_file.zip \
    --faults faults_nwh.json \
    --output tomography_model.xyz \
    --figures figures          # optional, needs matplotlib
```

Requires Python 3 + numpy. To change the result, edit `faults_nwh.json`: add or move
traces, change dips, depth limits, zone width, the depth profile of the reduction, or
the Vp/density coupling. Then re-run the script.

## Notes on the original model

These points are not changed by this script, but you should check them:

* The records in `tomography_model.xyz` are ordered with **northing varying fastest**,
  then easting, then depth. SPECFEM3D's `tomo` reader assumes x (easting) varies fastest.
* The tomography grid starts at easting 537 km (about 75.4°E). The mesh in `Mesh_Par_file`
  starts at 74.5°E (easting about 452 km). SPECFEM clamps points outside the tomography
  grid to its nearest edge values.
* The tomography top is z = 0 (sea level). SPECFEM uses the z = 0 values for topography
  above sea level, so the surface reduction carries up to the free surface.
