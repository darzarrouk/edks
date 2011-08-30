	subroutine read_parameters(parameters_filename)

	use parameters
	use wavenumber
	use constants
	implicit none

	character (len = 256)	:: parameters_filename

	! local variables
	integer, parameter	:: UNIT_PARAMETER = 8

	open (status = "old", unit = UNIT_PARAMETER, file = parameters_filename)

	! Model
	read(UNIT_PARAMETER, err = 20, end = 20, fmt = *) model_filename

	! Integration parameters
	read(UNIT_PARAMETER, err = 20, end = 20, fmt = *) facmin
	read(UNIT_PARAMETER, err = 20, end = 20, fmt = *) facmax
	read(UNIT_PARAMETER, err = 20, end = 20, fmt = *) HH

	close(UNIT_PARAMETER)
	return

 20	print *, "Error: reading the parameters file"

	end subroutine read_parameters
	
	subroutine read_model(model_filename)
	
	use model
	implicit none
	
	character (len=256)	:: model_filename
	
	! local variables
	 
	integer, parameter	:: UNIT_MODEL = 9
	integer			:: j
	real			:: factor
	
	open (status = "old", unit = UNIT_MODEL, file = model_filename)
	read (UNIT_MODEL, end = 10, FMT = *)  N, factor

!	if(N > MXMD) then
!		print *, "Error: The model in file: ", model_filename,        &
!                         "       does not fit into the allocated variables;", &
!                         "       you should change the dimensions \n"
!		stop 
!	endif 

!	Format of the model file changed July 1 - 2000
!	read (UNIT_MODEL, end = 10, fmt = *)  (vp0(j), j=1,N)
!	read (UNIT_MODEL, end = 10, fmt = *)  (vs0(j), j=1,N)
!	read (UNIT_MODEL, end = 10, fmt = *)  (rh0(j), j=1,N)
!	if(N > 1) then
!		read (UNIT_MODEL, end = 10, fmt = *)  (th0(j), j=1,N-1)
!		th0 = th0*factor
!	endif

	allocate(vp0(N), vs0(N), rh0(N), th0(N))
	do j=1,N
		read (UNIT_MODEL, end = 10, fmt = *)  rh0(j), vp0(j), vs0(j), th0(j)
	end do
		
	close(UNIT_MODEL)

	vp0 = vp0*factor
	vs0 = vs0*factor
	rh0 = rh0*factor
	th0 = th0*factor

	return
	
 10	print *, "Error reading the model file: ", model_filename
	stop 

	end subroutine read_model
