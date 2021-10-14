# Multicycle_dynamic_SSAF_NSJF

## The repository hosts software, fault geometric structure, model, and result for the manuscript 'Observation-constrained multicycle dynamic models of the southern San Andreas and the northern San Jacinto Faults: addressing complexity in paleoearthquake extents and recurrences with realistic fault geometry.' by Dunyu Liu, Benchun Duan, Katerine Scharer, and Doug Yule, which is submitted to JGR Solid Earth.

## ./exe/ 
The direcotry contains the executable of EQdyna version 2.0.2 that is used to produce the Model A, B, and C. With input files located in the same directory, eqdyna2d_2.0.2 should run on linux systems.

## ./fault_geometry/
The directory contains the raw CFM5.2 data used in the project. The MATLAB script calls functions under the subdirectory function/ to convert the UTM coordinates of the fault to the x-y-z coordinates in EQdyna.

## ./input/ 
The directory contains input files required by EQdyna.

## ./result/ 
The directory contains the MATLAB post-processing scripts to show the majority of results in the manuscript. 
Mesh files and results (mesh_simulation_result.zip) are too large for GitHub and they are stored via Pangaea https://issues.pangaea.de/browse/PDI-29895. After unzipping them under the ./result/ directory, the scripts should run and produce plots. 
