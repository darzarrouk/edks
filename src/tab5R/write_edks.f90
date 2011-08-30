!	include		'../sum_layered/edks.inc'
	subroutine write_edks(elem_filename, edks)
	use defs_module
	implicit none

!	Mandatory input parameters
	character (len=256)	elem_filename
	type (edks_type)	edks

! 	local variables
	integer jh, jr, jl

	open(unit=10, file='hdr.'//elem_filename, status='new', form='formatted')
		write(10,*) edks % prefix
		write(10,*) edks % nlayer

		do jl=1, edks % nlayer
			write(10,*) 	edks % rho(jl), 	&
					edks % alpha(jl), 	&
					edks % beta(jl), 	&
					edks % thickness(jl)
		end do
		write(10,*) edks % date
		write(10,*) edks % version
		write(10,*) edks % comment
		write(10,*) edks % depthmin, edks % depthmax,  edks % ndepth
		write(10,*) edks % distamin, edks % distamax,  edks % ndista
	close(10)

	open(unit=20, file=elem_filename, status='new', form='unformatted')
	write(20) ((edks % depths(jh), edks % distas(jr), edks % zrtdsx(jh,jr,1:10), jr=1,edks % ndista), jh = 1,edks % ndepth)

	close(20)
	return

	end subroutine write_edks
