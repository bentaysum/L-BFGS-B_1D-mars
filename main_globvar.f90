MODULE main_globvar


IMPLICIT NONE 

! TLM/Adjoint Structure
INTEGER nqmx, nlayermx, ndt
REAL lt_init
! Order of tracers in the TLM
CHARACTER(len=15), ALLOCATABLE :: noms(:)
! Mixing Ratio Vectors from the 1D model in TLM format (ndt,nlayermx*nqmx)
REAL, ALLOCATABLE :: pq_c(:,:)
! Mean Molar Mass in 1-D model 
REAL, ALLOCATABLE :: mmean(:,:)
! Tangent Linear/Adjoint Matrix and the Transition matrix P_N
REAL*8, ALLOCATABLE :: TLM(:,:,:), ADJ(:,:,:)
! Sensitivity Vector 
REAL*8, ALLOCATABLE :: hatJ(:,:)
! Time domain limits [forecast time and backtrace time]
INTEGER t_N, t_0
! Index of forecast element 
INTEGER J_idx

END MODULE main_globvar





