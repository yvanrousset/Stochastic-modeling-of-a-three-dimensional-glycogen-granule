# Stochastic modeling of a three dimensional glycogen granule

This repository contains jupyter notebook scripts to reproduce all the figures presented in the article "Stochastic modeling of a three-dimensional glycogen granule and impact of the branching enzyme on the structure"

Some functions will be re-defined in some scripts while there could be already defined in the main class from `glycogen_module.py`. It allows some flexibility in the notebooks and the plots.

The parameters file `parameters.json`, is required to initiate a glycogen granule. Although most of the time, parameters will be set inside the different scripts.

By default, the number of simulations as well as the size of the granule is kept low to avoid excessive computational time. We made sure that simulation time, by default, does not exceed 2 or 3 min.

