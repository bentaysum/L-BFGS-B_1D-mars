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
INTEGER N ! loop iterator 












! =========================================================
! Stage One: Selection of the 1-D model files to be studied 
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
DO N = 1, N_files
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



END