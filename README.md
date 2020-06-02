# Optimization Routines involving the 1-D Photochemistry Model Adjoint

This repository holds routines that aim to use the L-BFGS-B optimization routines
to work with Curiosity Rover data and the 1-D photochemistry model and it's 
adjoint, developed by Benjamin M. Taysum.

## Stage One

Create a routine that reads in two atmospheric states, a control and a 
perturbed 1-D model run [nqmx,nlayermx] in shape, to create an atmospheric state
vector in the same shape as in the adjoint/tangent linear model code i.e.
[nlayermx*nqmx]. 

### Begins on 02/06/2020
Stage One now acquires the relevant tlm bin and text and control ncdf file names
from command line prompts in main.f90.
