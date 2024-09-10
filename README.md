# QARSEN_lite

Quantum Algorithms for Real Space simulation of Electrons and Nuclei

This repo contains algorithm demos which uses [pyQuEST](https://github.com/rrmeister/pyQuEST) for quantum simulation of particles in real-space grids. It contains the essential elements which went into the numerical simulations reported in [this work](https://arxiv.org/abs/2202.05864).

## Requirements
- pyQuEST
- numpy
- jupyter
- matplotlib

## Installation
More detailed instructions can be found in the first notebook
1. Install [miniconda](https://docs.conda.io/en/latest/miniconda.html) if you don't have it already.
```bash
bash ./install_miniconda.sh
```

2. Create new python environment: 
```bash
# change dir to grid-based-quantum repo
cd grid-based-quantum
# create pyquest conda env
conda create -n -y pyquest python=3.9
# then activate the conda environment with
conda activate pyquest
```

3. Install PyQuEST: 
```bash
pip install ./pyQuEST
```

4. Install necessary dependencies by running:
```bash
pip install numpy==1.26.4 jupyter matplotlib
```

## Examples

[Getting Started](https://github.com/QARSEN-QC/qarsen_lite/blob/main/Getting_started.ipynb) -
Basic methods of defining wavepackets in a grid, using `PhaseFunc` from `pyquest`, extracting phase.

[Custom Potential](https://github.com/QARSEN-QC/qarsen_lite/blob/main/Custom_potentials.ipynb) -
Time propagation of quantum wavepackets under a custom potential, using `pyQuEST` methods.

[Ground state preparation](https://github.com/QARSEN-QC/qarsen_lite/blob/main/Ground_state_preparation.ipynb) -
Driving a loaded wavepacket towards the ground state subject to the interaction potentials, using `pyQuEST` methods to implement PITE.




