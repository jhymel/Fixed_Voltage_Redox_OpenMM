# Equilibrate Monolayer
Example of using a Monte Carlo based method for equilibrating electrolyte between the monolayer and electrodes.

Requires cuda Version 10.0, an OpenMM conda environment, and a proper starting pdb file.

In order to run on PACE:
```
module load cuda/10.0
module load anaconda3
conda activate OpenMM
python run_openMM.py example_monolayer.pdb
```

A pbs script is also included so that this can be run using:
```
qsub run.pbs
```
Though modifications to run.pbs will most likely be needed to run in your PBS environment.

Warning: This script includes a relative path pointing to the head of this repository. If run_openMM.py is moved up or down a directory this will break. After cloning this repository it might be a good idea modify run_openMM.py so that the path to the head of the repo is hard coded in.

The input for this run script, example_monolayer.pdb, is a combination of the outputs created in the build_monolayer and build_electrolyte directories. The monolayer_electrolyte_combiner.py is included for combining these, instructions for it's use are included in the file's header.

