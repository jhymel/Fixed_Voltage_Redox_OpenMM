# Randomize monolayer
After the monolayer and electrolyte have been equilibrated using the tools in [equilibrate_monolayer](https://github.com/jhymel/Fixed_Voltage_Redox_OpenMM/tree/main/examples/construct_monolayer_system/equilibrate_monolayer), one further issue needs to be resolved. When the monolayer is created using build_monolayer.py in [build_monolayer](https://github.com/jhymel/Fixed_Voltage_Redox_OpenMM/tree/main/examples/construct_monolayer_system/build_monolayer), the location of ferrocene-terminated chains is highly ordered since the script just loops over lattice sites modulo the Fc11/C11 ratio. This ordering is unphysical and needs to be dealt with prior to further simulations.

The method we've implemented to randomize/disorder the system involves heating to 600 Kelvin(K) for 20 nanoseconds(ns) followed by cooling the system from 600 K to 300 K at a rate of 10 K/ns.

This should only be run after equilibrating the system in [equilibrate_monolayer](https://github.com/jhymel/Fixed_Voltage_Redox_OpenMM/tree/main/examples/construct_monolayer_system/equilibrate_monolayer). After that, the checkpoint file, state.chk, needs to be copied to this directory in order for the simulation to start from the equilibrated geometry.

Two run scripts are included in this directory:
 - simulated_heating.py
   - Heat system to 600 K for 20 ns in order to randomize.
 - simulated_annealing.py
   - Cool system from 600 K to 300 K at a rate of 10 K/ns in order to bring system back to standard working temperature.

Requires cuda Version 10.0, an OpenMM conda environment, starting pdb file, and initial checkpoint file.

In order to run on PACE:
```
module load cuda/10.0
module load anaconda3
conda activate OpenMM
python simulated_heating.py example_monolayer.pdb
python simulated_annealing.py example_monolayer.pdb
```

A pbs scripts are also included so that simulated_heating.py and simulated_annealing.py can be run using:
```
qsub simulated_heating.pbs
qsub simulated_annealing.pbs
```
Though modifications to the pbs scripts will most likely be needed to run in your PBS environment.

Since simulated_annealing.py should only really be run after simulated_heating.py has completed, a better way to run these is as dependent pbs jobs. A very short bash script, run.sh, has been included which does this. It can be run using:
```
./run.sh
```

Warning: This script includes a relative path pointing to the head of this repository. If simulated_heating.py or simulated_annealing.py are moved up or down a directory they will break. After cloning this repository it might be a good idea modify those files so that the path to the head of the repo is hard coded in.
