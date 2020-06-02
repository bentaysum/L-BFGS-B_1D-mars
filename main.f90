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
INTEGER N ,iq, t ! loop iterators 
INTEGER a,b ! loop iterators 

REAL*8, ALLOCATABLE :: TLM_double(:,:)










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



10 DO N = 1, N_files
    READ(10,'(a)',iostat = iostat) NAME_LIST(N)
    ! Display in terminal for user selection 
    WRITE(*,"(I5,A20)") N,  TRIM(NAME_LIST(N))
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
! ALLOCATE(TLM_double(nqmx*nlayermx,nqmx*nlayermx))
ALLOCATE(TLM(nqmx*nlayermx,nqmx*nlayermx,ndt))
ALLOCATE(ADJ(nqmx*nlayermx,nqmx*nlayermx,ndt))


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


END