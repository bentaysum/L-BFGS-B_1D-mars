# Optimization Routines involving the 1-D Photochemistry Model Adjoint

This repository holds routines that aim to use the L-BFGS-B optimization routines
to work with Curiosity Rover data and the 1-D photochemistry model and it's 
adjoint, developed by Benjamin M. Taysum.

## Log 


### 05/06/2020 

netCDF output routine in works. Will output the optimised input state, and it's
forward modelled (via the Tangent Linear Model currently) output state, in the same
format as the 1-D model.

### 04/06/2020 

Routine produces optimized profiles! 

Work on making the routine more interactive. Investigate the new optimised profiles
properly in python, and assess validity with the fully integrated 1-D model.

### 03/06/2020

Optimization routine fails to work. Poorly defined cost function and gradient is probably the
source of the issue. Assess mathematics and ensure we feed the routines the correct variables.

