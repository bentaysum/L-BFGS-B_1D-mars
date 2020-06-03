# Optimization Routines involving the 1-D Photochemistry Model Adjoint

This repository holds routines that aim to use the L-BFGS-B optimization routines
to work with Curiosity Rover data and the 1-D photochemistry model and it's 
adjoint, developed by Benjamin M. Taysum.

## Log 

### 03/06/2020

Optimization routine fails to work. Poorly defined cost function and gradient is probably the
source of the issue. Assess mathematics and ensure we feed the routines the correct variables.

