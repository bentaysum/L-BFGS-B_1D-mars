PROGRAM main 
!
! 02/06/2020 : Ben Taysum
! 
! Routine aims to use a 1-D model control state, 
! a perturbation state, and the TLM/adjoint of the control,
! to minimize a cost function to produce an ideal atmospheric
! state that will best replicate the O2 measurements from
! the Curiosity Rover.

use netcdf 
use main_globvar

IMPLICIT NONE 

!======================================
! Head Directory of all 1-D model files
!======================================
CHARACTER(len=200), PARAMETER :: head_dir = "/exports/csce/datastore/geos/users/s1215319/paper2/curiosity_1Dfiles/"
! ----------------------------------------
! Sub-directories of the controls and tlms
! ----------------------------------------
CHARACTER(len=50), PARAMETER :: control_dir = "control_file/"
CHARACTER(len=50), PARAMETER :: tlm_dir = "tlm_files/"
! --------------------------------------
! Names of files (read at command line)
! --------------------------------------
CHARACTER(len=20), ALLOCATABLE :: NAME_LIST(:)

CHARACTER(len=20) CONTROL_NCDF 
CHARACTER(len=20) TLM_TEXT, TLM_BIN

! Local Variables
! ===============
INTEGER N_files ! Number of files in control directory
INTEGER N_choice ! Choice of file
CHARACTER(len=100) dummy 
INTEGER iostat ! Reading error integer 
INTEGER N_f ,iq, t, i ! loop iterators 
INTEGER a,b ! loop iterators 

REAL*8 Curiosity_O2_mmr ! NEED TO CREATE A ROUTINE TO READ REAL DATA IN

REAL*8, ALLOCATABLE :: dPQ(:) ! Perturbation vector for TLM calculations 

! ==================
! L-BFGS-B Variables
! ==================
INTEGER,PARAMETER :: nmax = 1024 ! Dimension of largest problem to be solved
INTEGER,PARAMETER :: mmax = 17 ! Maximum number of limited memory connections 
CHARACTER(len=60) task ! First entry = START 
                       ! Return of task(1:2) = "FG", user must evaluate function f and gradient g
                       !        at the returned value of x 
                       ! Return of task(1:5) = "NEW_X", an iteration of the algorithm has concluded
                       !        and f and g contain f(x) and g(x) respectively.
                       ! " "    of task(1:4) = "CONV", termination test is satisfied 
                       ! " "    of task(1:4) = "ABNO", termination ended without satisfying termination
                       !        test conditions, x contains best approximation found, and f and g hold 
                       !        f(x) and g(x) respectively. 
                       ! " "    of task(1:5) = "ERROR", routine found errors in input parameters.
CHARACTER(len=60) csave ! character working array 
LOGICAL lsave(4) ! On exit where TASK(1:4) = "NEW_X":
                 ! lsave(1) = .true. - initial x did not satisfy bounds 
                 ! lsave(2) = .true. - problem contains bounds 
                 ! lsave(3) = .true. - each variable has upper and lower bounds
INTEGER n ! Number of variables 
INTEGER m ! Number of corrections used in the limited memory matrix [range 3 <= m <= 20 is recommended]
INTEGER iprint ! Frequency + type of output from the L-BFGS-B routine 
               ! iprint < 0 : none 
               ! iprint = 0 : one line at last iteration 
               ! 0 < iprint < 99 : f and |proj g| at each iteration 
               ! iprint = 99 : details of every iteration bar n-vectors 
               ! iprint = 100 : changes of active set and final x 
               ! iprint > 100 : details of every iteration including x and g 
               ! for iprint > 0, iterate.dat summarises the iteration 
INTEGER,ALLOCATABLE :: nbd(:), iwa(:), isave(:)
REAL*8 f ! value of the function at the point x 
REAL*8 factr ! tolerance for routine to stop Optimization.
             ! routine will stop for :
             !
             ! (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
             !
             ! where epsmch is the machine precision auto set by code. 
             ! factr = 1.D+12 low accuracy 
             !       = 1.D+7 moderate accuracy 
             !       = 1.D+1 high accuracy
REAL*8 pgtol ! Iteration will stop when: 
             ! 
             ! max{|proj g_i | i = 1, ..., n} <= pgtol
             !
             ! and can be suppressed as a stopping condition by settin pgtol = 0

REAL*8,ALLOCATABLE :: x(:),x0(:), l(:), u(:), g(:), dsave(:), wa(:)
! l(:) : lower bounds of all variables (length n)
! u(:) : uppwer bounds of all variables (length n) 
! x(:) : set by user as the initial estimate of the solution vector (length n)
! g(:) : components of the gradient at the point x 

! =========================================================
! Stage 1.0: Selection of the 1-D model files to be studied 
! =========================================================

! Control file, PQ_c
call system('ls ' // trim(head_dir) // trim(control_dir) // ' > controlname_holding.txt')
open(10,FILE="controlname_holding.txt",ACTION="READ")
! Acquire number of files in control Directory
N_files = 0
DO 
    READ(10,'(a)', iostat = iostat) dummy 
    IF (iostat .ne. 0) EXIT 
    N_files = N_files + 1
ENDDO 
CLOSE(10)

ALLOCATE(NAME_LIST(N_files))
! Read the text file again and save the file names in NAME_LIST
open(10,FILE="controlname_holding.txt",ACTION="READ")
write(*,"(a25)")"========================="



10 DO N_f = 1, N_files
    READ(10,'(a)',iostat = iostat) NAME_LIST(N_f)
    ! Display in terminal for user selection 
    WRITE(*,"(I5,A20)") N_f,  TRIM(NAME_LIST(N_f))
   ENDDO

! Choice of file 
20 write(*,"(a25)")"========================="
   write(*,*)"Control file to operate with: "
   read(*,"(I2)",iostat=iostat) N_choice

   IF ( (iostat .ne. 0 ) .or. &
        (N_choice .lt. 1) .or. &
        (N_choice .gt. N_files) ) THEN 
        30 WRITE(*,"(A,I2,A)") "Integer 1 =< i =< ", N_files, " required" 
        GOTO 20 
   ENDIF

CLOSE(10)

CONTROL_NCDF = TRIM(NAME_LIST(N_choice))
TLM_BIN = TRIM(CONTROL_NCDF(1: LEN_TRIM(CONTROL_NCDF) - 3 )) // "_tlm.bin"
TLM_TEXT = TRIM(CONTROL_NCDF(1: LEN_TRIM(CONTROL_NCDF) - 3 )) // "_tlm.txt"

! -------------------------------------------------------
! Stage 1.1 : Reading the TLM text file and allocation of
!             relevant global variables 
! -------------------------------------------------------
OPEN(10,FILE=TRIM(head_dir)//TRIM(tlm_dir)//TRIM(TLM_TEXT),ACTION="READ",iostat=iostat)
IF ( iostat .ne. 0 ) THEN
    WRITE(*,*) "TEXT FILE FOR TLM MISSING"
    WRITE(*,*) TRIM(head_dir)//TRIM(tlm_dir)//TRIM(TLM_TEXT)
    GOTO 10
ELSE 
    WRITE(*,*) "Structure File : ", TRIM(head_dir)//TRIM(tlm_dir)//TRIM(TLM_TEXT)
ENDIF 

! Scratch first two line 
READ(10,"(a15,F6.3)") dummy, lt_init
READ(10,"(a)") dummy
! Total number of model timesteps, tracers, and model layers 
READ(10,"(a15,I5)") dummy, ndt
READ(10,"(a15,I5)") dummy, nqmx
READ(10,"(a15,I5)") dummy, nlayermx
! Allocate size of variables in routine 
ALLOCATE(noms(nqmx))
ALLOCATE(pq_c(ndt,nqmx*nlayermx))
pq_c(:,:) = 0.E0 
ALLOCATE(TLM(nqmx*nlayermx,nqmx*nlayermx,ndt))
ALLOCATE(ADJ(nqmx*nlayermx,nqmx*nlayermx,ndt))
ALLOCATE(mmean(ndt,nlayermx))

! Discard next line
READ(10,"(a)") dummy

! Get order of tracers in the TLM model Structure
DO iq = 1,nqmx
    READ(10,"(a10,a5)") noms(iq), dummy
ENDDO
CLOSE(10)

! -------------------------------------------------------
! Stage 1.2 : Construct the control state vector 
! -------------------------------------------------------
call oneDmgcm_reader(TRIM(head_dir)//TRIM(control_dir)//TRIM(CONTROL_NCDF))


! =======================================================
! Stage 2.0 : Adjoint Model
! =======================================================
WRITE(*,*) "BINARY TLM FILE : ", TRIM(head_dir)//TRIM(tlm_dir)//TRIM(TLM_BIN) 
OPEN(50,FILE = TRIM(head_dir)//TRIM(tlm_dir)//TRIM(TLM_BIN), ACCESS = 'direct', &
        RECL = nqmx*nlayermx*nqmx*nlayermx*8)

DO t = 1, ndt 
    READ(50,rec=t) ( ( TLM(a,b,t), a = 1, nqmx*nlayermx ), b = 1, nqmx*nlayermx )
    ! Adjoint is the TLM transpose 
    write(*,"(F6.2,A1)") 100.*REAL(t)/REAL(ndt), "%"
    ADJ(:,:,t) = TRANSPOSE(TLM(:,:,t))

ENDDO 


call adjoint_1D

! ===================================
! Stage 3.0 : Optimization Procedures
! =================================== 

! Allocation of the variables used in L-BFGS-B 
n = nqmx*nlayermx 
m = 7

ALLOCATE(nbd(nmax))
ALLOCATE(iwa(3*nmax))
ALLOCATE(isave(44))
ALLOCATE(x(nmax))
ALLOCATE(l(nmax))
ALLOCATE(u(nmax))
ALLOCATE(g(nmax))
ALLOCATE(dsave(29))
ALLOCATE(wa(2*mmax*nmax + 5*nmax + 11*mmax*mmax + 8*mmax))

ALLOCATE(dPQ(nqmx*nlayermx))

! Specification of the bounds of the variables.


iprint = 1 
factr = 1.0D+1
pgtol = 0.!1.0D-5

DO i = 1, n 

	! First guess of initial conditions PQ_c in x(:)
	x(i) = PQ_c(t_0,i)*1.D0 
	
	IF ( x(i) < 1.D-30 ) THEN 
		x(i) = 1.D-30
		
		l(i) = 1.D-30 
		u(i) = 1.D-30 
		
		nbd(i) = 2
		CONTINUE 
		
	ENDIF 
	
	l(i) = MAX( 1.D-30, x(i)*5.D-1)
	u(i) = MIN( 9.9D-1, x(i)*1.5D0)
	
	nbd(i) = 2 
	
		
ENDDO 


! NEED REAL ROUTINE TO BE MADE 
Curiosity_O2_mmr = 0.002160D0*(16./mmean(t_N,1)) 

task = 'START'

i = 1

111 CONTINUE 

call setulb(n,m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint, &
            csave,lsave,isave,dsave)
		

IF ( task(1:2) .eq. "FG" ) THEN
	
	IF ( i == 1 ) THEN 
		call costfunction(curiosity_O2_mmr, PQ_c(t_0,:)*1.D0, f) 
		g = hatJ(t_0,:)/Curiosity_O2_mmr 
	
	ELSE 
		call costfunction(curiosity_O2_mmr, x(:n), f) 
		g = hatJ(t_0,:)/Curiosity_O2_mmr 
	ENDIF 
	i = i + 1
	goto 111
	
ENDIF 

IF ( task(1:5) .eq. "NEW_X" ) goto 111 

WRITE(*,*) MAXVAL( x(:n) - PQ_c(t_0,:) ) 

write(*,*) PQ_c(t_N,J_idx), Curiosity_O2_mmr
write(*,*) (PQ_c(t_0,J_idx) + DOT_PRODUCT( hatJ(t_0,:) , x(:n) - PQ_c(t_0,:) ) ) , Curiosity_O2_mmr


write(*,"(16A15)") ( trim(noms(iq)), iq = 1, 16 ) 

DO i = 1, nlayermx
	WRITE(*,"(16F15.7)") ( 100.D0*(x( (iq-1)*nlayermx + i )/PQ_c(t_0, (iq-1)*nlayermx + i ) - 1.D0 ), iq = 1, 16 )
ENDDO 


STOP 

END 


SUBROUTINE costfunction( curiosity_O2, PQ_i, f )

! Calculates the value of the optimizable cost function:
!
! COST_i = [1-D Model Calculated O2] + ( hatJ , PQ_i - PQ_c )
!	   - [Curiosity O2 Measurement]
!
! The value of the gradient at all times is, at the moment, simply
! hatJ, the sensitivity vector calculated by the Adjoint model.
use main_globvar
 
IMPLICIT NONE 

! Input
! ============
REAL*8 curiosity_O2, PQ_i(nqmx*nlayermx)
! Local 
! ============
REAL*8 dPQi(nqmx*nlayermx)

! Output 
! ============
REAL*8 f ! Cost function value


! ===================================================================

dPQi = PQ_i - PQ_c(t_0,:)

! f = PQ_c(t_0,J_idx) + DOT_PRODUCT( hatJ(t_0,:) , dPQi ) &
	! - curiosity_O2
	
f = PQ_c(t_0,J_idx) + DOT_PRODUCT( hatJ(t_0,:) , dPQi )  
	
f = ABS((f/curiosity_O2) - 1.D0)

	
END







