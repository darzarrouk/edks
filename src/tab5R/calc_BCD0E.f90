	subroutine calc_BCD0E(q, w, v, qE, wE, B, C, D, E, F)
	
	use nrtype
	use nr
	use wavenumber
	use parameters
	use constants

	implicit none
	real, dimension(0:2,nk)	:: q, w, v
	real, dimension(nk)	:: qE, wE
	real			:: B(0:2), C(0:2), D(3:4), E, F

	B(0) = sum(w(0,:)*kv)
	B(1) = 0.
	B(2) = 0.
	C(0) = 0.
	C(1) = sum((q(1,:)-v(1,:))*kv)/2.
	C(2) = 0.
	D(3) = C(1)
	D(4) = 0.

	E    = sum(wE(:)*kv)
	F    = 0.

	end subroutine calc_BCD0E
