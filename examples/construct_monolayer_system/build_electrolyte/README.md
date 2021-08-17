# Build Electrolyte
Example of using PackMol to create an electrolyte consisting of tip3p water and NaCl ions.
The pdb file output by this script is intended to be combined with a pdb containing the electrodes and surfactant (created by build_monolayer.py).

In order to run use:
```
./packmol < pack.inp
```

Default parameters are set such that...
 - Box dimensions
   - xmin = 2.0 , ymin = 2.0 , zmin = 25.0
   - xmax = 29.0 , ymax = 29.0 , zmax = 88.0
 - Molecules
   - 2000 H2O
   - 100 Na
   - 100 Cl

These parameters can be changed by modifying pack.inp
