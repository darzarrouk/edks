!	include	"modules.inc"

	recursive subroutine kenpsv(n,m, k, Ruac, Tuac, Rdac, Tdac)

	use constants
	implicit none

	integer 		:: n, m
	real			:: k
	real, dimension(2,2)	:: Ruac, Tuac, Rdac, Tdac

	!local variables
	real, dimension(2,2)	:: dummy, TEMP
	real, dimension(2,2)	:: Ruab, Tuab, Rdab, Tdab
	real, dimension(2,2)	:: Rubc, Tubc, Rdbc, Tdbc
	real, dimension(2,2)	:: Rum,  Tum,  Rdm,  Tdm, Emm1

!	print *, n, m

	if( m < n .or. m < 1 .or. n < 0 ) then
	        print *, "function 'kenpsv' called with n=%d and m=%d\n", n, m
	        print *, "Error: 'n' should be '>=' than '0', \n"
	        print *, "Error: 'm' should be '>=' than '1'  and \n"
	        print *, "Error: 'n' should be '<=' than 'm' \n"
	        stop
	endif

	if ((m == n) .or. (m == 1 .and. n == 0)) then
	        Tdac  = Id2
	        Rdac  = Ze2
	        Tuac  = Id2
	        Ruac  = Ze2
	        if(n == 0) then
	!               when called with n=0; only 'Ruac' is 
	                call RERPSV(1,dummy,Ruac)
	        endif
	else
		call ERTUDPSV(m, k, Emm1,Rum,Tum,Rdm,Tdm)
	
	        Rubc  = Rum
	        Tubc  = matmul(Emm1,Tum)
	        Rdbc  = matmul(Emm1, matmul(Rdm,Emm1))
	        Tdbc  = matmul(Tdm,Emm1)
	
	        call kenpsv(n,m-1,k,Ruab, Tuab, Rdab, Tdab)

                TEMP  = matmul(Ruab,Rdbc)
	        Tdac  =      matmul(matmul(Tdbc,inv22(Id2-TEMP)),Tdab)
	        Rdac  = Rdab+matmul(matmul(matmul(Tuab,Rdbc),inv22(Id2-TEMP)),Tdab)
                TEMP  = matmul(Rdbc,Ruab)
	        Tuac  =      matmul(matmul(Tuab,inv22(Id2-TEMP)),Tubc)
	        Ruac  = Rubc+matmul(matmul(matmul(Tdbc,Ruab),inv22(Id2-TEMP)),Tubc)
	endif

	contains

	function inv22(a) result (ia)

	use constants
	implicit none
	real, dimension(2,2)	:: a, ia

	!local variables
	real			:: det

	det = a(1,1)*a(2,2)-a(1,2)*a(2,1)

	if( det < DET_MIN) then
		print *, "Error: inv22 called with a singular matrix"
		stop
	endif

	ia(1,1) =  a(2,2)
	ia(1,2) = -a(1,2)
	ia(2,1) = -a(2,1)
	ia(2,2) =  a(1,1)

	ia = ia/det

	end function inv22

	end subroutine kenpsv

	recursive subroutine kensh(n,m, k, Ruac, Tuac, Rdac, Tdac) 

	use constants
	implicit none

	integer 		:: n, m
	real			:: k
	real			:: Ruac, Tuac, Rdac, Tdac

	!local variables
	real  			:: dummy
	real  			:: Ruab, Tuab, Rdab, Tdab
	real  			:: Rubc, Tubc, Rdbc, Tdbc
	real  			:: Rum,  Tum,  Rdm,  Tdm, Emm1

	!print *, "function 'kensh(n,m)' called with n=%d and m=%d\n", n, m
	if( m < n .or. m < 1 .or. n < 0 ) then
	        print *, "function 'kensh' called with n=%d and m=%d\n", n, m
	        print *, "Error: 'n' should be '>=' than '0', \n"
	        print *, "Error: 'm' should be '>=' than '1'  and \n"
	        print *, "Error: 'n' should be '<=' than 'm' \n"
	        stop
	endif
	
	if ((m == n) .or. (m == 1 .and. n == 0)) then
	        Tdac  = 1
	        Rdac  = 0
	        Tuac  = 1
	        Ruac  = 0
	        if(n == 0) then
	!               when called with n=0; only 'Ruac' is 
	                call RERSH(dummy, Ruac)
	        endif
	else
	        call ERTUDSH(m, k, Emm1,Rum,Tum,Rdm,Tdm)
	        Rubc  = Rum
	        Tubc  = Emm1*Tum
	        Rdbc  = Emm1*Rdm*Emm1
	        Tdbc  = Tdm*Emm1
	
	        call kensh(n,m-1, k, Ruab, Tuab, Rdab, Tdab)
	        Tdac  = Tdbc/(1.-Ruab*Rdbc)*Tdab
	        Rdac  = Rdab+Tuab*Rdbc/(1-Ruab*Rdbc)*Tdab
	        Tuac  = Tuab/(1.-Rdbc*Ruab)*Tubc
	        Ruac  = Rubc+Tdbc*Ruab/(1-Rdbc*Ruab)*Tubc
	endif

	end subroutine kensh
