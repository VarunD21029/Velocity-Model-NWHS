# Velocity_Model_NWHS
This repository contains the SPECFEM files for the 1905 Kangra Simulations
Current files be used with SPECFEM to generate results from the work: _"A Regional 3D Crustal Velocity Model for the Northwestern Himalayas and Scenario Simulations of the 1905 Kangra Earthquake"_

## Fault damage zones
`velocity_model_file_faults.zip` is the same tomography model with low-velocity fault damage zones added. The thrusts form one connected system rooted in the Main Himalayan Thrust: an MFT frontal ramp, an upper flat, a mid-crustal ramp beneath the Higher Himalaya and a lower flat. The JMT, MBT and MCT root into it. The STD, Kaurik-Chango and Karakoram faults are modelled separately. Vs is reduced by 10 % at the surface traces, and the reduction continues down-dip along each fault plane. See `fault_zones/README.md` for the geometry, the parameters and how to regenerate the file.
