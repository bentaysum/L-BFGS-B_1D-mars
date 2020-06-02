SUBROUTINE oneDmgcm_reader(filestring)

use main_globvar

IMPLICIT NONE 

include "netcdf.inc"

! Input Variables
! ===============
CHARACTER(LEN=200) filestring 

! Local Variables 
! ===============
INTEGER ncid
INTEGER retval
INTEGER iq
INTEGER varid 
CHARACTER(len=10) tracer

! Open the netCDF file defined by filestring for read-only access
retval = nf_open(filestring, NF_NOWRITE, ncid)
if (retval .ne. nf_noerr) call handle_err(retval)

! Iterate through the tracers in the order that they appear in 
! the noms global variable to ensure they match with the order
! of the TLM/Adjoint matrix
do iq = 1,nqmx
    tracer = noms(iq)

    ! Get the varid of the tracer mmr variable, based on its name.
    retval = nf_inq_varid(ncid, trim(tracer), varid)
    if (retval .ne. nf_noerr) call handle_err(retval)

    ! Read into the relevant location of the PQ_C[ndt,nqmx*nlayermx] variable
    retval = nf_get_var_REAL(ncid, varid, PQ_c(:,(iq-1)*nlayermx + 1 : iq*nlayermx) )
    if (retval .ne. nf_noerr) call handle_err(retval)

    write(*,*) tracer , MAXVAL(PQ_c(:,(iq-1)*nlayermx + 1 : iq*nlayermx) )


enddo 



STOP

END SUBROUTINE




subroutine handle_err(errcode)
implicit none
include "netcdf.inc"

integer errcode

print *, 'Error: ', nf_strerror(errcode)
stop 2
end



