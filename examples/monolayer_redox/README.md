# Monolayer Redox - The meat and potatoes of examples.
Example of how to compute the change in free energy during the transfer of an electron from a redox molecule to the electrode using thermodynamic integration and OpenMM.

#### Background/theory:
We want to compute the change in free energy when a neutral ferrocene (Fc) becomes oxidized (Fc<sup>+</sup>) within our system. Since energy is a state function, the change in free energy is simply the difference in energy between the two states.

![equation](https://latex.codecogs.com/svg.latex?%7B%5Ccolor%7BBlue%7D%20%5CDelta%20A%5E0_%7BFc%20%5Crightarrow%20Fc%5E&plus;%7D%20%3D%20A%5E0_%7BFc%5E&plus;%7D%20-%20A%5E0_%7BFc%7D%7D)

We need some path, thermodynamicially, between these two states. If we introduce a parameter, λ, which couples the two states such that,

![equation](https://latex.codecogs.com/svg.latex?%7B%5Ccolor%7BBlue%7D%20%5Clambda%20%3D%20%5Cleft%5C%7B%5Cbegin%7Bmatrix%7D%201%20%5C%20at%20%5C%20Fc%5E&plus;%20%5C%5C%200%20%5C%20at%20%5C%20Fc%20%5Cend%7Bmatrix%7D%5Cright.)

then the change in free energy can be represented as,

![equation](https://latex.codecogs.com/svg.latex?%5Ccolor%7BBlue%7D%7B%20%5CDelta%20A%5E0_%7BFc%20%5Crightarrow%20Fc%5E&plus;%7D%20%3D%20A%5E0%28%5Clambda%20%3D%201%29%20-%20A%5E0%28%5Clambda%3D0%29%7D)

Now that we are computing a difference between states connected by a coupling parameter, λ, we can employ the machinary of statistical mechanics to represent this as a integral over λ.

![equation](https://latex.codecogs.com/svg.latex?%7B%5Ccolor%7BBlue%7D%20%5CDelta%20A%5E0_%7BFc%20%5Crightarrow%20Fc%5E&plus;%7D%20%3D%20%5Cint_%7B0%7D%5E%7B1%7D%5Cleft%5Clangle%20%5Cfrac%7B%5Cdelta%20A%28%5Cvec%7B%5CGamma%7D%2C%5Clambda%29%7D%7B%5Cdelta%5Clambda%7D%20%5Cright%5Crangle_%5Clambda%20d%5Clambda)

From this, to move along our path between the neutral and oxidized states, all we need to do is linearly scale λ. But what is λ actually changing?

If we assume that the only difference between the neutral and oxidized states is the charges on the atoms, then all λ needs to be is a scaling factor for the charges, q, between the two states.

![equation](https://latex.codecogs.com/svg.latex?%7B%5Ccolor%7BBlue%7D%20%5Cvec%7Bq%7D_%7Bredox%7D%28%5Clambda%29%20%3D%20%281-%5Clambda%29%5C%3A%20%5Cvec%7Bq%7D_%7BFc%7D%20&plus;%20%5Clambda%5C%3A%20%5Cvec%7Bq%7D_%7BFc%5E&plus;%7D)

What we end up doing in practicallity is run molecular dynamics simulations using a series of scaled charges corresponding to λ = [0.0, 0.1, ... , 0.9, 1.0]. Throughout each simulation we compute δA/δλ using a numerical derivative in the form of

![equation](https://latex.codecogs.com/svg.latex?%7B%5Ccolor%7BBlue%7D%20%5Cfrac%7B%5Cdelta%20A%28%5Cvec%7B%5CGamma%7D%2C%5Clambda%29%7D%7B%5Cdelta%5Clambda%7D%20%5Capprox%20%5Cfrac%7BA%28%5Cvec%7B%5CGamma%7D%2C%5Clambda&plus;%5CDelta%5Clambda%29-A%28%5Cvec%7B%5CGamma%7D%2C%5Clambda%29%7D%7B%5CDelta%5Clambda%7D%7D)

By averaging these derivatives from throughout each simulation we get the ensemble average derivative for each value of λ. Then by numerically integrating over these with respect to λ, we can get the total free energy difference.
