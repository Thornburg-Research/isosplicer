# isosplicer
Module to simulate mRNA splicing reaction dynamics for the human genome. Reaction dynamics are simualted using the well-stirred stochastic reaction solvers in Lattice Microbes.

This module includes the capability to simulate mRNA for the vast majority of the human genome.

## Contents
| File-Directory | Description |
|----------|-------------|
| `model/` | Python executables to run the reaction model for mRNA splicing |
| `human_genome/` | Preproccessing instructions and code to download and extract gene, transcript, and exon information from the human genome |
| `publication_simulations/` | Parameter sets and executables used to generate trajectories for publication |
| `isosplicer.yml` | Yaml dependencies file for creating the conda environemnt to run executables in ```model/``` |

## Installation
This is a Python-based model that uses a conda environment to manage libraries.

I recommend working in the main directory for this repository for local installation.

First, create a conda environment from the .yml file provided here. Experienced conda users can change the environment name as desired.

```
conda env create -f isosplicer.yml
```

Once the environment is created, activate your environment.

```
conda activate isosplicer
```

Download/clone Lattice Microbes from the developer: https://github.com/Luthey-Schulten-Lab/Lattice_Microbes

Once you have ```Lattice-Microbes/``` locally, make a build directory and move into it. (Location is not important, but I recommend making it in this directory.)

```
mkdir build_isosplicer
cd build_isosplicer
```

In the build directory, we will install Lattice Microbes.  You will need to use the PATH to where your local Lattice-Microbes repo is located.

```
cmake /PATH/Lattice-Microbes/src/
make
make install
```

If you want to verify installation, you can run the following command and should get a list of ```lm``` input options..

```
lm --help
```

No further dependencies are required to run isosplicer simulations. You are ready to go!
