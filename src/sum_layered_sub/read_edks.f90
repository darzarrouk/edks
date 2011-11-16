!	include		'edks.inc'
	subroutine read_edks(elem_filename, edks)
	use defs_module
	implicit none

!	Mandatory input parameters
	character (len=256)	elem_filename
	type (edks_type)	edks

! 	local variables
	integer jh, jr, jl

	open(unit=10, file='hdr.'//elem_filename, status='old', form='formatted', err=510)
		read(10,*) edks % prefix
		read(10,*) edks % nlayer

		allocate(edks % rho	 (edks % nlayer))
		allocate(edks % alpha	 (edks % nlayer))
		allocate(edks % beta	 (edks % nlayer))
		allocate(edks % thickness(edks % nlayer))

		do jl=1, edks % nlayer - 1
			read(10,*) 	edks % rho(jl), 	&
					edks % alpha(jl), 	&
					edks % beta(jl), 	&
					edks % thickness(jl)
		end do
		read(10,*) 	edks % rho(edks % nlayer), 	&
				edks % alpha(edks % nlayer), 	&
				edks % beta(edks % nlayer) 	
		read(10,*) edks % date
		read(10,*) edks % version
		read(10,*) edks % comment
		read(10,*) edks % depthmin, edks % depthmax,  edks % ndepth
		read(10,*) edks % distamin, edks % distamax,  edks % ndista
	close(10)

	allocate(edks % depths(edks % ndepth))
	allocate(edks % distas(edks % ndista))
	allocate(edks % zrtdsx(edks % ndepth, edks % ndista, 10))

	open(unit=20, file=elem_filename, status='old', form='unformatted')
	read(20) ((edks % depths(jh), edks % distas(jr), edks % zrtdsx(jh,jr,1:10), jr = 1,edks % ndista), jh = 1,edks % ndepth)
	close(20)
	return
 510    write(*,'(a)') 'Error opening: hdr.'//elem_filename
        call exit(1)
 520    write(*,'(a)') 'Error opening: '//elem_filename
        call exit(1)
	end subroutine read_edks
