# Coarse-grained Simulation using OpenMM
## Introduction 

This project contains the source code that reproduces the data of the following two papers: 

- [Wang, Yanming, and Charles L. Brooks III. "Enhanced Sampling Applied to Modeling Allosteric Regulation in Transcription." The journal of physical chemistry letters 10.19 (2019): 5963-5968.](https://doi.org/10.1021/acs.jpclett.9b03618)

- [Wang, Yanming, and Charles L. Brooks. "Electrostatic Forces Control the Negative Allosteric Regulation in a Disordered Protein Switch." The Journal of Physical Chemistry Letters (2020).](https://doi.org/10.1021/acs.jpclett.9b02226)

This project contains C++ code running coarse-grained model simulations based on OpenMM (http://openmm.org/). <br>
Now it supports the Electrostatics Inclusive Karanicolas-Brooks (EIKB) Gō model (Karanicolas, John, and Charles L. Brooks III. Protein Science 11.10 (2002): 2351-2361.) <br>

The main feature of this project is that it supports Hamiltonian replica exchange method with highly flexible Hamiltonian by utilizing the flexible OpenMM API.
