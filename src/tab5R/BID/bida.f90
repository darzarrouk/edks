	program bid
	
	USE nrtype; USE nrutil, ONLY : poly
	implicit none

	integer			:: nk
	real(sp)		:: kval
	real(sp), allocatable	:: kv(:), j2(:)
	REAL(SP), allocatable   :: ax(:),xx(:),z(:)
	REAL(DP), allocatable   :: y(:)
	LOGICAL(LGT), allocatable :: mask(:)
	REAL(DP), DIMENSION(5) :: p = (/1.0_dp,-0.1098628627e-2_dp,&
		0.2734510407e-4_dp,-0.2073370639e-5_dp,0.2093887211e-6_dp/)

	kval = 9.
	nk=89835
	allocate(j2(nk), kv(nk), ax(nk), xx(nk), z(nk), y(nk),mask(nk))
	kv    = kval
	mask = (abs(kv) < 8.0)
	where (mask)
		y=kv**2
		j2=poly(y,p,mask)/poly(y,p,mask)
	elsewhere
		ax=abs(kv)
		z=8.0_sp/ax
		y=z**2
        	j2 = poly(y,p,.not. mask)    ! LINE A
        	j2 = poly(y,p,.not. mask)    ! LINE B
	end where
	deallocate(j2, kv, ax, xx, z, y,mask)
end program bid
