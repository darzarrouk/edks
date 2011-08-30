!	include "modules.inc"
	
	subroutine ERTUDPSV(jn, k, EM, RU, TU, RD, TD)
	use model
	implicit none
	integer 		:: jn
	real			:: k
	real, dimension(2,2) 	:: EM, RU, TU, RD, TD
	
	! local variables
	real			:: mu_1, Delta_1, d1
	real			:: mu_2, Delta_2
	
	real			:: ex1, Id2(2,2)
	real			:: a, b, c, d, e, f, DR
	
	real			:: rppp, rpsp, rspp, rssp
	real			:: rppm, rpsm, rspm, rssm
	real			:: tppp, tpsp, tspp, tssp
	real			:: tppm, tpsm, tspm, tssm
	data Id2(1,:)/1., 0./, Id2(2,:)/0., 1./
	
	mu_1    = mu(jn-1)
	Delta_1 = delta(jn-1)
	d1      = th(jn-1)
	mu_2    = mu(jn)
	Delta_2 = delta(jn)
	
	
	ex1 	= exp(-k*d1)
	EM   	= ex1*Id2
	
	a 	=                 (mu_1/mu_2+Delta_2)/(1+Delta_2)
	b 	= -2*Delta_1*k*d1*(mu_1/mu_2+Delta_2)/(1+Delta_2)
	c 	=         (mu_1*Delta_1/mu_2-Delta_2)/(1+Delta_2)
	d 	=               (mu_1*Delta_1/mu_2+1)/(1+Delta_2)
	e 	=                       (mu_1/mu_2-1)/(1+Delta_2)
	f 	=        2*Delta_1*k*d1*(mu_1/mu_2-1)/(1+Delta_2)
	
	DR   	=  a*d
	
	rppp 	= b*e/DR
	rpsp 	= -a*e/DR
	rspp 	= -(d*c-b*f)/DR
	rssp 	= -a*f/DR
	RD(1,:)	= (/rppp, rspp/)
	RD(2,:)	= (/rpsp, rssp/)
	
	rppm 	= 0.
	rpsm 	= d*e/DR
	rspm 	= a*c/DR
!	rssm 	= (-a*f-b*e)/DR
	rssm 	= 0.
	RU(1,:)	= (/rppm, rspm/)
	RU(2,:)	= (/rpsm, rssm/)
	
	tppp 	=  a*(a*d-c*e)/DR
!	tpsp 	=  e*(b*e+a*f)/DR
 	tpsp 	=  0.
	tspp 	= -a*(b*d+c*f)/DR
	tssp 	= ((a*d-c*e)*d+(b*e+a*f)*f)/DR
	TD(1,:) = (/tppp, tspp/)
	TD(2,:) = (/tpsp, tssp/)
	
	tppm 	= d/DR
	tpsm 	= 0.
	tspm 	= -b/DR
	tssm 	= a/DR
	TU(1,:)	= (/tppm, tspm/)
	TU(2,:)	= (/tpsm, tssm/)
	END subroutine ERTUDPSV
	
	subroutine ERTUDSH(jn, k, EM, RU, TU, RD, TD)
	use model
	implicit none
	integer 		:: jn
	real			:: k
	real			:: EM, RU, TU, RD, TD
	
	! local variables
	real			:: mu_1, d1
	real			:: mu_2
	real			:: DL
	
	mu_1    = mu(jn-1)
	d1      = th(jn-1)
	mu_2    = mu(jn)
	
	DL      =  mu_1+mu_2
	
	EM      = exp(-k*d1)
	RD      = (mu_1-mu_2)/DL
	TD      = 2*mu_1/DL
	RU      = -RD
	TU      = 2*mu_2/DL
	
	END subroutine ERTUDSH
	
	subroutine RERPSV(jn, REV, R)
	use model
	implicit none
	integer			:: jn
	real, dimension(2,2)	:: REV, R
	
	! local variables
	real			:: idelta
	
	idelta	 = 1./delta(jn)
	REV(1,:) = (/ idelta, -1. /)
	REV(2,:) = (/ idelta,  1. /)
	REV	 = REV*(1.+delta(jn))
	
	R(1,:)	 = (/ 0.,  -delta(jn) /)
	R(2,:)	 = (/-idelta,	 0.  /)
	
	END subroutine RERPSV
	
	subroutine RERSH(REV, R)
	implicit none
	real			:: REV, R
	
	REV	 = 2.
	R	 = 1.
	
	END subroutine RERSH
