!	include 'modules.inc'
	subroutine calc_coeffsE(q, w, v, qE, wE, r, B, C, D, E, F)
	
	use nrtype
	use nr
	use wavenumber
	use parameters
	use constants

	implicit none
	real, dimension(0:2,nk)	:: q, w, v
	real, dimension(nk)	:: qE, wE
	real			:: r
	real			:: B(0:2), C(0:2), D(3:4), E, F

	! local variables
	real			:: tmp1, tmp2, tmp3
	real, allocatable	:: j0(:), j1(:), j2(:)

	allocate(j0(nk), j1(nk), j2(nk))
	j0 = bessj0( kv*r)
	j1 = bessj1( kv*r)
	j2 = bessj(2,kv*r)
	B(0) =    sum(   w(0,:)        *j0*kv)
	C(0) =   -sum(   q(0,:)        *j1*kv)

        B(1) =    sum(   w(1,:)        *j1*kv)
        tmp1 =    sum((  q(1,:)+v(1,:))*j1)/r
        tmp2 =    sum(kv*q(1,:)        *j0)
        tmp3 =    sum(kv*       v(1,:) *j0)
        C(1) =   -tmp1 + tmp2
	D(3) =    tmp1 - tmp3

        B(2) =    sum(   w(2,:)        *j2*kv)
        tmp1 =  2*sum((  q(2,:)+v(2,:))*j2)/r
        tmp2 =    sum(kv*q(2,:)        *j1)
        tmp3 =    sum(kv*       v(2,:) *j1)
        C(2) =   -tmp1 + tmp2
	D(4) =    tmp1 - tmp3

	E    =    sum(  wE(:)        *j0*kv)
	F    =   -sum(  qE(:)        *j1*kv)
	deallocate(j0, j1, j2)

	end subroutine calc_coeffsE
