# Fixed_Voltage_Redox_OpenMM

General use scripts, libraries, and forcefield files used to calculate free-energies associated with the oxidation/reduction of molecules within a fixed-voltage molecular dynamics simulation. Examples are included for reduction of ferrocene-terminated alkanethiols within an alkanethiol monolayer on gold.

#### Directories of interest
 - examples
   - Contains example run scripts for many different simulations as well as scripts for constructing alkanethiol monolayers
 - lib
   - Contains classes and methods used in the example simulations.
   - Classes inherit from (Fixed_Voltage_OpenMM)[https://github.com/jmcdaniel43/Fixed_Voltage_OpenMM/tree/master] submodule hosted by (Prof. Jesse McDaniel)[https://github.com/jmcdaniel43].
 - ffdir
   - Contains xmls files used in the example simulations.
   - Same ffdir used for systems simulating fixed and flexible surfactants.
