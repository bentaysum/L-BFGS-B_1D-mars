MODULE main_globvar


IMPLICIT NONE 

! TLM/Adjoint Structure
INTEGER nqmx, nlayermx, ndt
REAL lt_init
! Order of tracers in the TLM
CHARACTER(len=15), ALLOCATABLE :: noms(:)
! Mixing Ratio Vectors 
REAL, ALLOCATABLE :: pq_c(:,:)
! Tangent Linear/Adjoint Matrix
REAL, ALLOCATABLE :: TLM(:,:,:), ADJ(:,:,:)
! Sensitivity Vector 
REAL, ALLOCATABLE :: hatJ(:,:)

END MODULE main_globvar





