SUBROUTINE adjoint_1D

use main_globvar 

IMPLICIT NONE 

! Local Variables
! ===============
INTEGER NSTEP 
INTEGER STEP0 
INTEGER iostat 
REAL, ALLOCATABLE :: lt_array(:)
INTEGER t, iq 

INTEGER sol_N, sol_0 ! Sol of forecast and backtrace 
REAL lt_N, lt_0 ! Local time of forecast on sol_N and sol_0
INTEGER t_N, t_0 ! Forecast and backtrace timestep 

CHARACTER(len=10) J_TRACER 
INTEGER J_LAYER 
INTEGER J_idx 








! Asks the user which timestep out of NDT they would like to 
! trace the adjoint model backwards towards AND from, what tracer they
! wish to study, and at what model layer.
WRITE(*,"(A20)") "==================================================="
WRITE(*,"(A20)") "Ajoint Model"
WRITE(*,"(A20)") "==================================================="

! First, construct an array of local times corresponding to each timestep 
ALLOCATE(lt_array(ndt)) 

! OPERATES UNDER THE ASSUMPTION THERE ARE 48 TIMESTEPS PER MARTIAN SOL IN 
! THE 1-D MODEL 

! ndt = 0 ; lt_array(t=0) = lt_init
lt_array(1) = lt_init 


DO t = 2, ndt 
    lt_array(t) = lt_array(t-1) + 0.5 
    if ( lt_array(t) >= 24.) lt_array(t) = lt_array(t) - 24. 
ENDDO 

! =====================
! FORECAST TIME : t = N 
! =====================
10 WRITE(*,*) "FORECAST SOL : "
READ(*,*,iostat = iostat) sol_N
IF (iostat .ne. 0) THEN 
    WRITE(*,*) "INTEGER ONLY"
    GOTO 10
ENDIF

t_N = sol_N*48 

IF (t_N > ndt ) THEN 
    WRITE(*,*) "TOO MANY SOLS SELECTED; NOT COVERED BY FORECAST"
    GOTO 10     
ENDIF 

WRITE(*,*) "FORECAST LOCAL TIME : "
READ(*,*,iostat = iostat) lt_N 
IF (iostat .ne. 0) THEN 
    WRITE(*,*) "FLOAT ONLY"
    GOTO 10
ENDIF

t_N = t_N + MINLOC( ABS(lt_array(t_N:t_N + 48) - lt_N) , dim = 1) - 1

! ======================
! BACKTRACE TIME : t = 0
! ======================
20 WRITE(*,*) "BACKTRACE SOL : "
READ(*,*,iostat = iostat) sol_0
IF (iostat .ne. 0) THEN 
    WRITE(*,*) "INTEGER ONLY"
    GOTO 20
ENDIF

t_0 = sol_0*48 

IF (t_0 > t_N ) THEN 
    WRITE(*,*) "BACKTRACE _PROCEEDS_ FORECAST ELEMENT; CHOOSE ANOTHER VALUE"
    GOTO 20     
ENDIF 

WRITE(*,*) "BACKTRACE LOCAL TIME : "
READ(*,*,iostat = iostat) lt_0 
IF (iostat .ne. 0) THEN 
    WRITE(*,*) "FLOAT ONLY"
    GOTO 20
ENDIF

t_0 = t_0 + MINLOC( ABS(lt_array(t_0:t_0 + 48) - lt_0) , dim = 1) - 1

WRITE(*,"(A12,F5.2,I4)") "Forecast: ", lt_array(t_n), t_n 
WRITE(*,"(A12,F5.2,I4)") "Backtrace: ", lt_array(t_0), t_0


! ========================================
! Initialisation of the sensitivity vector
! ======================================== 
!
ALLOCATE(hatJ( t_N-t_0, nqmx*nlayermx))

! Tracer Choice 
30 WRITE(*,*) "Forecast Tracer :"
READ(*,*,iostat=iostat) J_TRACER
IF ( iostat .ne. 0 ) THEN 
    WRITE(*,*) "TRACER NAME AS APPEARS IN 1-D MODEL"
    GOTO 30 
ENDIF 
! Layer Choice 
40 WRITE(*,*) "Forecast Layer :"
READ(*,*,iostat=iostat) J_LAYER
IF ( iostat .ne. 0 ) THEN 
    WRITE(*,*) "TRACER LAYER AS INTEGER"
    GOTO 40 
ENDIF 

! Find index in TLM of forecast element J_l^t 
hatJ(:,:) = 0.E0 
DO iq = 1,nqmx
    IF ( trim(noms(iq)) == trim(J_TRACER) ) THEN 
        
        J_idx = (iq-1)*nlayermx + J_LAYER 
        
        hatJ(t_N-t_0, J_idx) = 1. 
        
        GOTO 50  
    
    ENDIF 

    IF ( iq == nqmx ) THEN 
        WRITE(*,*) "poorly chosen tracer/layer"
        GOTO 30
    ENDIF 

ENDDO 

50 WRITE(*,*) "Sensitivity vector of ", trim(J_TRACER), " at layer ", J_LAYER, " initialised..."

! Engage backtrace of the adjoint 
DO t = (t_N-t_0) - 1, 1, -1 
    hatJ(t,:) = MATMUL(ADJ(t,:,:),hatJ(t+1,:))

    write(*,*) hatJ(t+1, J_idx), MAXVAL(ADJ(t,:,:)), MINVAL(ADJ(t,:,:))

ENDDO 



END 