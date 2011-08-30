!	include "modules.inc"
	subroutine discotabE(nsource, Pp, Pm, SVp, SVm, SHp, SHm, Ep, Em)

	use model
	implicit none

	integer 		:: nsource
	real, dimension(0:2)	:: Pp,  Pm, SVp, SVm
	real, dimension(0:2)	:: SHp, SHm
	real			:: Ep, Em

	!local variables
	real			:: dd, eps, epd, tmp

	dd  = 	delta(nsource)
	tmp =   (-1.+4.*dd)/2.

	eps =  1.
	epd = eps*dd
	Pp  =  (/ tmp,	-epd, 0.5 /)  ! en accord avec la formule 16 du papier
	Pp  =  Pp/(1.+dd)
	SVp =  (/ -1.5,        	eps,    -0.5 /)
	SVp =  SVp/(1.+dd)
	SHp =  (/  0.0,        	eps,    -1.0 /)
	Ep  =  1.0

	eps = -1.
	epd = eps*dd
	Pm  =  (/ tmp,	-epd, 0.5 /) ! en accord avec la formule 16 du papier
	Pm  =  Pm/(1.+dd)
	SVm =  (/ -1.5,        	eps,    -0.5 /)
	SVm =  SVm/(1.+dd)
	SHm =  (/  0.0,        	eps,    -1.0 /)
	Em  =  1.0

	end subroutine discotabE
