MODULE main_globvar


IMPLICIT NONE 

! TLM/Adjoint Structure
INTEGER, SAVE :: nqmx, nlayermx, ndt
! Order of tracers in the TLM
CHARACTER(len=15), ALLOCATABLE :: noms(:)
! Mixing Ratio Vectors 
REAL, ALLOCATABLE :: pq_c(:,:)

END MODULE main_globvar





