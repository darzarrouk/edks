!	include 'modules.inc'

	subroutine calcqwvE(Pp, Pm, SVp, SVm, SHp, SHm, Ep, Em, q, w, v, qE, wE)
	
	use model
	use wavenumber
	use parameters
	use constants
	implicit none
	real, dimension(0:2)	:: Pp, Pm, SVp, SVm, SHp, SHm
	real			:: Ep, Em
	real, dimension(0:2,nk)	:: q, w, v
	real, dimension(nk)	:: qE, wE

	! local variables
	real, dimension(2,2)	:: RUFS,   TURS,   RDRS,   TDRS
	real			:: RUFSsh, TURSsh, RDRSsh, TDRSsh
	real, dimension(2,2)	:: RUSL,   TUSL,   RDSL,   TDSL
	real  			:: RUSLsh, TUSLsh, RDSLsh, TDSLsh
	real, dimension(2,2) 	:: REV,    RTIL
	real 			:: REVsh,  RTILsh
	real, dimension(2,2)	:: MOUT, TEMP1, TEMP2
	real			:: MOUTsh
	integer			:: jm,jk
	real, dimension(2,1)	:: VD,qw
	real, dimension(2,1)	:: VDE,qwE
	real			:: k,VDsh

	!unusefulls
	real, dimension(2,2)	:: TUFS,   RDFS,   TDFS,   RURS
	real			:: TUFSsh, RDFSsh, TDFSsh, RURSsh

	do 	jk=1,nk
		k = kv(jk)

		call kenpsv(0,   M+1,  k, RUFS,   TUFS,   RDFS,   TDFS  ) 
        	call kensh (0,   M+1,  k, RUFSsh, TUFSsh, RDFSsh, TDFSsh) 

		call kenpsv(1,   M+1,  k, RURS,   TURS,   RDRS,   TDRS  ) 
        	call kensh (1,   M+1,  k, RURSsh, TURSsh, RDRSsh, TDRSsh) 

        	call kenpsv(M+1, Np+1, k, RUSL,   TUSL,   RDSL,   TDSL  )
        	call kensh (M+1, Np+1, k, RUSLsh, TUSLsh, RDSLsh, TDSLsh)

        	call RERPSV(1, REV, RTIL)
        	REVsh  = 2.
		RTILsh = 1.
                TEMP1 = matmul(RDRS, RTIL)
                TEMP2 = matmul(RDSL, RUFS)

        	MOUT   =  matmul(matmul(matmul(REV, inv22(Id2 - TEMP1)), 		&
				    	     TURS), inv22(Id2 - TEMP2))
        	MOUTsh =  REVsh / (1. - RDRSsh * RTILsh) * TURSsh / (1. - RDSLsh * RUFSsh)

		do 	jm=0,2
 			VD       = matmul(RDSL,	reshape( (/ Pp(jm),SVp(jm) /), (/ 2, 1 /))) + 	&
						reshape( (/ Pm(jm),SVm(jm) /), (/ 2, 1 /)) 
                	VDsh     = RDSLsh*SHp(jm) + SHm(jm)

                	qw   	 = matmul(MOUT,VD)
                	q(jm,jk) = qw(1,1)
                	w(jm,jk) = qw(2,1)
                	v(jm,jk) = MOUTsh*VDsh

		enddo
 		VDE      = matmul(RDSL,	reshape( (/ Ep, 0. /), (/ 2, 1 /))) + 	&
					reshape( (/ Em, 0. /), (/ 2, 1 /)) 
               	qwE  	 = matmul(MOUT,VDE)
               	qE(jk) = qwE(1,1)
               	wE(jk) = qwE(2,1)
	enddo

!	open(unit = 8, file="kqwv.m", status='unknown')
!	write(8, '(10(e15.9,1x))') (kv(jk), q(:,jk), w(:,jk), v(:,jk), jk=1,nk)

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

	end subroutine calcqwvE
