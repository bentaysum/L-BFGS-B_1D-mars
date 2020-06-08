SUBROUTINE optimised_out( pq_0 ) 

use main_globvar
use netcdf 

IMPLICIT NONE 

! Input ! 
! ===== !
REAL*8 pq_0(nqmx*nlayermx) ! Optimized input mixing ratios 

! Local !
! ===== !
CHARACTER(len=100) FILE_NAME 

INTEGER ncid
INTEGER NY, y_dimid ! Individual tracer mmr's
INTEGER NY_full, y_full_dimid ! Full vector 
INTEGER full_varid, varid(nqmx) 

INTEGER iq, lyr ! Loop iterators 
INTEGER dimids(1)

FILE_NAME = "trial.nc"


! Establish sizes of vectors
NY = nlayermx 
NY_full = nlayermx*nqmx

! Create the file, or overwrite the one already present 
call check( nf90_create(TRIM(FILE_NAME), NF90_CLOBBER, ncid) )

! Define the dimensions. NetCDF will hand back an ID for each. 
call check( nf90_def_dim(ncid, "y", NY, y_dimid) )
call check( nf90_def_dim(ncid, "y_full", NY_full, y_full_dimid))

! Loop over all tracers and define their mixing ratio vectors
dimids = (/y_full_dimid/)
call check( nf90_def_var(ncid, "pq_0", NF90_DOUBLE, dimids , full_varid) )

DO iq = 1, nqmx
	call check( nf90_def_var(ncid, trim(noms(iq)) // "_0", NF90_DOUBLE, y_dimid, &
				varid(iq)) )
ENDDO 

! End define mode. This tells netCDF we are done defining metadata.
call check( nf90_enddef(ncid) )

! Write the data 
! ==============
! Full optimised state vector 
call check( nf90_put_var(ncid,full_varid,pq_0) ) 
! Individual chunks
DO iq = 1,nqmx
	call check( nf90_put_var(ncid,varid(iq), pq_0( (iq-1)*nlayermx + 1 : iq*nlayermx ) ) ) 
ENDDO 

call check( nf90_close(ncid) )

return 







contains

  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  
  
  
END SUBROUTINE 
