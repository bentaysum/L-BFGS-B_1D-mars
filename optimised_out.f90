SUBROUTINE optimised_out( pq_0, pq_N ) 

use main_globvar

IMPLICIT NONE 

include 'netcdf.inc'

! Input Variables
! ===============
REAL*8 pq_0(nqmx*nlayermx) ! Optimised input tracer mixing ratios 
REAL*8 pq_N(nqmx*nlayermx) ! pq_0 forward modelled to the forecast time

! Local Variables
! ===============
REAL*8 pq0_iq(nlayermx), pqN_iq(nlayermx) ! Individual tracers 

CHARACTER(len=200) NETCDF_NAME ! NetCDF output name 
INTEGER retval ! Error return

INTEGER iq ! Loop iterator 

INTEGER varid(nqmx) ! Variable ID's for each tracer 
INTEGER t_0id, t_Nid ! Temporal index id's
INTEGER ncid ! file netCDF ID 

INTEGER nx, x_dimid ! Dimensions of the output
INTEGER dimids 

nx = nlayermx

! name of the netcdf file 
NETCDF_NAME = "optimized_state.nc" 


! Construct the file 
retval = nf_create(NETCDF_NAME, NF_CLOBBER, ncid)
if (retval .ne. nf_noerr) call handle_err(retval)

! Define the dimensions. NetCDF will hand back an ID for each. 
retval = nf_def_dim(ncid, "x", NX, x_dimid)
if (retval .ne. nf_noerr) call handle_err(retval)

dimids = x_dimid 
! ===============================================================
! Define the PQ	variables, breaking the pq_0 and pq_N vectors into 
! individual tracer vectors as per the 1-D model output (ease of 
! comparison via Python analysis)
! ===============================================================
DO iq = 1, nqmx 
	
	write(*,*) trim(noms(iq)) // "_0" 
	retval = nf_def_var(ncid, trim(noms(iq)) // "_0", &
					NF_DOUBLE,  dimids, varid(iq))
	if (retval .ne. nf_noerr) call handle_err(retval)
	
	
ENDDO 

! Define the forecast and backtrace temporal indices t_0 and t_N 
retval = nf_def_var(ncid, "t_0", &
			NF_INT, 0, t_0id)
if (retval .ne. nf_noerr) call handle_err(retval)

retval = nf_def_var(ncid, "t_N", &
			NF_INT, 0, t_Nid)
if (retval .ne. nf_noerr) call handle_err(retval)


! ===============
! End define mode 
! ===============
retval = nf_enddef(ncid)
if (retval .ne. nf_noerr) call handle_err(retval)

! ===================================
! Putting the variables into the file
! ===================================
DO iq = 1, nqmx
	
	! Extract input 
	pq0_iq = pq_0( (iq-1)*nlayermx + 1 : iq*nlayermx ) 
	! Extrac forward cast state 
	! pqN_iq = pq_N( (iq-1)*nlayermx + 1 : iq*nlayermx )

	retval = nf_put_var_double(ncid, varid(iq), pq0_iq)
	if (retval .ne. nf_noerr) call handle_err(retval)

	! retval = nf_put_var_double(ncid, varid(iq), pqN_iq)
	! if (retval .ne. nf_noerr) call handle_err(retval)
	
ENDDO 
	

retval = nf_put_var_double(ncid, t_0id, t_0)
if (retval .ne. nf_noerr) call handle_err(retval)

retval = nf_put_var_double(ncid, t_Nid, t_N)
if (retval .ne. nf_noerr) call handle_err(retval)

! Close the file 
retval = nf_close(ncid)
if (retval .ne. nf_noerr) call handle_err(retval)


END SUBROUTINE 
