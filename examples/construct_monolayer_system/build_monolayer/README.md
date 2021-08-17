# Build monolayer
Example for building a monolayer of alkanethiols where a certain percentage are terminated with ferrocene.

Requires an active MDAnalyis conda enviroment in order to run.

In order to run on PACE:
```
module load anaconda3
conda activate MDAnalysis
python build_monolayer.py
```

Default parameters are set such that...
 - Target density = 0.046
 - Sufactants used are Fc11 and C11
 - Ferrocene percentage = 0.25

These can be changed by modifying build_monolayer.py.
