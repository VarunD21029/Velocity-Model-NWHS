# Fault damage zones in the NW Himalaya velocity model

This folder adds low-velocity fault damage zones to the tomography model in
`../velocity_model_file.zip`. The result is
`../velocity_model_file_faults.zip` (`velocity_model_file_faults/tomography_model.xyz`).
It is a drop-in replacement for the original file: copy it to `DATA/tomo_files/` and keep
`MODEL = tomo` in `Par_file`.

## What was changed

* **One connected thrust system.** The thrusts are modelled as a single system rooted
  in the Main Himalayan Thrust (MHT), the basal décollement:

  ```
  S                                                                        N
  MFT     JMT  MBT                    MCT                STD
  0 km \     \    \                     \                  \
  5 km  \_____\____\____                 \                  \ (ends 15 km)
                       ‾‾‾‾‾‾‾‾‾‾‾‾\       \
  ~11 km      upper flat (3.5°)      \      \
                         mid-crustal  \     \
                         ramp (15°)    \     \
  22 km                                 \_____\________________
                                          lower flat (6°) → 30 km
  ```

  | Element | Geometry |
  |---|---|
  | MFT frontal ramp | 30° from the surface to 5 km |
  | Upper flat (décollement) | ~3.5°, from 5 km down to ~10–12 km at the top of the ramp |
  | **Mid-crustal ramp** | 15°, from ~10–12 km down to 22 km. Its top lies 10 km south of the MCT surface trace, so it follows the front of the Higher Himalaya along strike rather than sitting a fixed distance from the MFT. That matters because the Sub- and Lesser Himalaya are much wider in the Kangra re-entrant than in Garhwal–Kumaon |
  | Lower flat | 6°, from 22 to 30 km beneath the Higher and Tethyan Himalaya |
  | JMT, MBT (splays) | listric (35–45° → 25°), rooting into the upper flat at ~6–10 km |
  | MCT (splay) | 35° → 25°, rooting into the MHT just below the ramp at ~23 km |

  The splays have no depth limit of their own. Each one continues until it meets the MHT
  and stops there, merging into the MHT damage zone. The script prints the depths where
  they join. `figures/mht_depth_map.png` shows the MHT depth contours and the
  surface projection of the ramp.

  The ramp position, dips and depths are all in the `mht_geometry` block of
  `faults_nwh.json`. This geometry follows the flat–ramp–flat MHT seen in balanced
  sections (Powers et al., 1998; Srivastava & Mitra, 1994) and in receiver-function
  images (Caldwell et al., 2013). The mid-crustal ramp is also where the Himalayan
  microseismicity belt lies, with events such as the 1991 Uttarkashi and 1999 Chamoli
  earthquakes.

* **Faults outside the thrust system** (each kept separate):

  | Fault | Geometry |
  |---|---|
  | South Tibetan Detachment (STD) | 25° N normal-sense detachment, to 15 km |
  | Kaurik-Chango Fault | 60° W normal fault, to 15 km |
  | Karakoram Fault | 80° NE strike-slip fault, to 25 km |

* **Trace sources.** The MFT, Kaurik-Chango and Karakoram traces come from the GEM
  Global Active Faults database (HimaTibetMap). I extended the MFT by hand toward Jammu.
  The JMT, MBT, MCT and STD traces were **placed by hand from published regional maps**,
  so they are only accurate to a few km. The MCT trace also sets where the ramp is,
  so better MCT coordinates move the ramp too. Replace any of these with your own
  digitised traces if you have them.
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

About 8 % of the grid nodes are affected. QC plots are in `figures/`:
* map slices at 0–25 km depth
* the MHT depth map
* three S–N cross-sections with the fault planes drawn on top

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
