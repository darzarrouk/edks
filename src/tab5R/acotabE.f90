	subroutine acotab(depth,htol)
	use model
	implicit none
	real			:: depth, htol
	
	! local variables
	
	integer, parameter	:: YES= 1
	integer, parameter	:: NON=-1
	integer			:: flag
	real			:: h, dh
	
	h = depth
	N = size(vp0)
	if(N < 1) then
		print *, "Nb. of layers < 0 !! \n"
		stop 
	endif
	
	if(h < 0) then
		print *, "Source depth < 0 !!  \n"
		stop
	endif
	
	if(htol < 0)  then
		print *, "Warning, h-tolerance < 0 \n"
		print *, "         using  100 mts  \n"
		htol = 100
	endif
	
	flag = 0
	if( h < htol) then
	    	M = 0
	    	flag = NON
	elseif(N == 1) then
		M = 1
		flag = YES
	else
		do M=1,N-1
			dh = h - th0(M)
			if(dh > htol) then
				h = dh
			elseif(dh > -htol) then
	!			the source is at one pre-exs. interface
				flag = NON
				exit
			else
				flag = YES
				exit
			endif
		enddo
		if(flag == 0)  then
	!       the source is in the half space
			M = N
			flag = YES
		endif
	endif
	
	if(flag == YES) then
		allocate (vp(N+1), vs(N+1), rh(N+1), th(N))
		vp(1:N)   = vp0
		vs(1:N)   = vs0
		rh(1:N)   = rh0
		if(N > 1) then
			th(1:N-1) = th0(1:N-1)
		endif

		if(M < N) then
			if(M < N-1)  th(M+2:N) = th(M+1:N-1)
			th(M+1)   = th(M) - h
		endif
		th(M)       = h
		vp(M+1:N+1) = vp(M:N)
		vs(M+1:N+1) = vs(M:N)
		rh(M+1:N+1) = rh(M:N)
		N           = N + 1
	else
		allocate (vp(N), vs(N), rh(N), th(N-1))
		vp = vp0
		vs = vs0
		rh = rh0
		th = th0(1:N-1)
	endif
	
	allocate (lambda(N), mu(N), delta(N), gamma(N))
	lambda  = rh*(vp**2 - 2*vs**2)
	mu 	= rh*vs**2
	delta 	= 1 - 2*vs**2/(vp**2+vs**2)
	gamma   = (lambda+mu)/(lambda+2*mu)
	Np = N - 1
	
	end subroutine acotab
