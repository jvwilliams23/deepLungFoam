# deepLungFoam

[![DOI](https://zenodo.org/badge/432129853.svg)](https://zenodo.org/badge/latestdoi/432129853)


Here, we have implemented a repository for coupling a 3D respiratory CFD model to a lumped-parameter model for downstream flow. This model accounts for pleural-pressure change with time and varying capacitance and resistance across the lungs during the period of respiration. Such models are widely used in literature, so there are many variations. Here we use the model from Oakes et al. [1].

Credits for code foundations are given to [OpenFOAM-phys-flow](https://github.com/KeepFloyding/OpenFOAM-phys-flow) by KeepFloyding (Andris Piebalgs).

## Installation and usage
There are two solvers. One for single phase flow (deepLungFoam) and one for particle-laden flow (deepLungMPPICFoam). Each one comes with compilation scripts in the directory. 
The single-phase flow solver is less developed than the MPPIC solver, so it may lack some features.

An example case is provided in tutorials/ which shows the boundary conditions and additional dictionaries required to run the simulation.

## References
[1] Oakes, J.M., Roth, S.C. and Shadden, S.C., 2018. Airflow simulations in infant, child, and adult pulmonary conducting airways. Annals of biomedical engineering, 46(3), pp.498-512.
