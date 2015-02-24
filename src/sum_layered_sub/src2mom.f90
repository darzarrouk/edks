!	Include	'edks.inc'
	Subroutine src2mom 					&
		 ( M, o_Mo, o_MomT, o_nu, o_mu,			&
		 o_Area, o_strike, o_dip, o_rake, o_St, o_Sp, o_deltaV 	)

	Use defs_Module
!	use f90_unix    	!needed on the NAG compiler
	Implicit None
	Optional o_Mo, o_MomT, o_nu, o_mu, 				&
		 o_Area, o_strike, o_dip, o_rake, o_St, o_Sp, o_deltaV

!	Method 1: 	Mo, MomT(1:6),                   nu, mu
!	Method 2: 	Area, strike, dip, rake, St, Sp, nu
!	Method 3: 	Area, strike, dip, rake, St
!	Method 4: 	Area, strike, dip,           Sp, nu
!	Method 5: 	DeltaV,                          nu

!	Mandatory output arguments
	Real, Intent(out)		:: M(6)

!	Optional input arguments
	Real, Intent(in)		:: o_Mo, o_MomT(1:6), o_mu, o_nu
	Real, Intent(in)		:: o_strike, o_dip
	Real, Intent(inout)		:: o_rake
	Real, Intent(inout)		:: o_St, o_Sp
	Real, Intent(in)		:: o_Area, o_DeltaV

!	MECHANISM: 	The  source can be given in several different ways
!
!			in all cases what is given is to be used as the moment
!			for each one of the nsrc sources
!
!	Method 1: 	Mo, MomT(1:6), i                 nu, mu
!	Method 2: 	Area, strike, dip, rake, St, Sp, nu
!	Method 3: 	Area, strike, dip, rake, St
!	Method 4: 	Area, strike, dip,           Sp, nu
!	Method 5: 	DeltaV,                          nu
!
!	Mo:		Coeff to be applied to MomT
!	MomT:		Moment Tensor (vector(1:6)
!	mu, nu:		Rheology used to calculate Mo*MomT
!
!	strike:		azimuth: common for all the point sources
!	dip:		dip:     common for all the point sources
!	rake:		rake:    common for all the point sources
!	Area:		Surface of the Fault
!	St:		Slip tangential to the fault
!	Sp:		Slip perpendicular to the fault
!	DeltaV:		Change in Volume
!
!	OUTPUT PARAMETERS
!	M(6)
!
!	parameters
	Real,	Parameter	:: M_PI 	= 3.14159265358979

!	local variables

!	Source
!	Real			:: M(6)
	Real			:: rad_strike, rad_dip, rad_rake

!	Others
	Integer 		:: srctype
	Real			:: saz, caz, s2az, c2az

!	Local instances of the optional arguments

	Real              	:: Mo, MomT(1:6), mu, nu
	Real              	:: strike, dip
	Real              	:: rake
	Real              	:: St, Sp
	Real              	:: Area, DeltaV

!	Temporary variables for the source geometry
	Real			:: nor(3), sli(3), Maki(3,3), la_mu

        !Counters
        Integer                 :: ix, iy

! Set default value for nu
      nu = 0.25

! Source Parametrization
! srctype
	srctype = 0
	If(         Present(o_Mo)   .And.      Present(o_MomT)   .And.      Present(o_mu)  .And. &
                    Present(o_nu)   .And. .not.Present(o_strike) .And. .not.Present(o_dip) .and. &
               .not.Present(o_rake) .And. .not.Present(o_Area)   .And. .not.Present(o_St)  .And. &
               .not.Present(o_Sp)   .And. .not.Present(o_DeltaV) ) Then
		Mo	= o_Mo
		MomT	= o_MomT
		mu	= o_mu
		nu	= o_nu
		srctype = 1
        Elseif(.not.Present(o_Mo)   .And. .not.Present(o_MomT  ) .And. .not.Present(o_mu ) .And. &
                    Present(o_nu)   .And.      Present(o_strike) .And.      Present(o_dip) .and. &
                    Present(o_rake) .And.      Present(o_Area  ) .And.      Present(o_St ) .And. &
                    Present(o_Sp)   .And. .not.Present(o_DeltaV)  ) Then
                Area 	= o_Area
		strike 	= o_strike
		dip 	= o_dip
		rake 	= o_rake
		St 	= o_St
		Sp	= o_Sp
                nu      = o_nu
		srctype = 2
	Elseif(.not.Present(o_Mo)   .And. .not.Present(o_MomT  ) .And. .not.Present(o_mu ) .And. &
               .not.Present(o_nu)   .And.      Present(o_strike) .And.      Present(o_dip) .and. &
                    Present(o_rake) .And.      Present(o_Area  ) .And.      Present(o_St ) .And. &
               .not.Present(o_Sp)   .And. .not.Present(o_DeltaV)) Then
		Area 	= o_Area
		strike 	= o_strike
		dip 	= o_dip
		rake 	= o_rake
		St 	= o_St
                srctype = 3
        Elseif(.not.Present(o_Mo)   .And. .not.Present(o_MomT  ) .And. .not.Present(o_mu ) .And. &
                    Present(o_nu)   .And.      Present(o_strike) .And.      Present(o_dip) .and. &
               .not.Present(o_rake) .And.      Present(o_Area  ) .And. .not.Present(o_St ) .And. &
                    Present(o_Sp)   .And. .not.Present(o_DeltaV)  ) Then
                Area 	= o_Area
		strike 	= o_strike
		dip 	= o_dip
		Sp	= o_Sp
                nu      = o_nu
                srctype = 4

        Elseif(.not.Present(o_Mo)   .And. .not.Present(o_MomT  ) .And. .not.Present(o_mu)  .And. &
                    Present(o_nu)   .And. .not.Present(o_strike) .And. .not.Present(o_dip) .and. &
               .not.Present(o_rake) .And. .not.Present(o_Area  ) .And. .not.Present(o_St ) .And. &
               .not.Present(o_Sp)   .And.      Present(o_DeltaV)  ) Then
                nu      = o_nu
                DeltaV 	= o_DeltaV
                srctype = 5
        Endif

	If( srctype == 0 ) Then
		Write(*,'(a)') 'The set of furnished parameters to "displa" doesn"t match'
		Write(*,'(a)') 'any of the predefined kind of source method specifications'
		Write(*,'(a)') 'Method 1: 	Mo, MomT(1:6),                   nu, mu'
		Write(*,'(a)') 'Method 2: 	Area, strike, dip, rake, St, Sp, nu'
		Write(*,'(a)') 'Method 3: 	Area, strike, dip, rake, St        '
		Write(*,'(a)') 'Method 4: 	Area, strike, dip,           Sp, nu'
		Write(*,'(a)') 'Method 5: 	DeltaV,                          nu'
		Call Exit(1)
	Endif

!         seismic moment
!         Internally M follows the Aki's convention
!         x: north,  y:  east,    z: down
!         (Aki = Herrman)
	If(srctype .Eq. 1) Then
!               Here     Aki    Harv
!               M(2) =   Mxx =  Mtt
!               M(3) =   Myy =  Mff
!               M(1) =   Mzz =  Mrr
!               M(6) =   Mxy = -Mtf
!               M(4) =   Mxz =  Mrt
!               M(5) =   Myz = -Mrf
                M    =  MomT
		M(6) = -MomT(6)
		M(5) = -MomT(5)
		! The correction for the rheologie
		! should be implemented (Dufumier and Rivera)
	Elseif(srctype .Lt. 5) Then
		If(srctype .Eq. 3) Then
			Sp 	= 0.
		Endif
		If(srctype .Eq. 4) Then
			rake 	= 0.
			St     	= 0.
		Endif
 		rad_strike 	= strike*M_PI/180.
		rad_dip 	= dip   *M_PI/180.
		rad_rake 	= rake  *M_PI/180.
		la_mu		= 2.*nu/(1.-2.*nu)
		! moment tensor of a dislocation
		!lambda * SD*(s.n)Id + mu*SD*s.n'+n.s')

		nor(1)	= -Sin(rad_dip)*Sin(rad_strike)
		nor(2)	=  Sin(rad_dip)*Cos(rad_strike)
		nor(3)	= -Cos(rad_dip)
		sli(1)  =  Cos(rad_rake)*Cos(rad_strike)+Cos(rad_dip)*Sin(rad_rake)*Sin(rad_strike)
		sli(2)  =  Cos(rad_rake)*Sin(rad_strike)-Cos(rad_dip)*Sin(rad_rake)*Cos(rad_strike)
		sli(3)	= -Sin(rad_rake)*Sin(rad_dip)
		Do ix=1,3
			sli(ix) = sli(ix)*St+nor(ix)*Sp
			Do iy=1,3
				Maki(ix,iy) = 0.
			End Do
		End Do

		Do ix=1,3
			Maki(ix,ix) = la_mu*Sp
			Do iy=1,3
				Maki(ix,iy) = Maki(ix,iy) + (nor(ix)*sli(iy)+nor(iy)*sli(ix))
			End Do
		End Do

		M(1) = Maki(3,3)
		M(2) = Maki(1,1)
		M(3) = Maki(2,2)
		M(4) = Maki(1,3)
		M(5) = Maki(2,3)
		M(6) = Maki(1,2)
                M    = Area*M
	Else

	! moment tensor of an explosion
	! positive for explosions
		M(1) 	= DeltaV*2.*(1.+nu)/3./(1.-2.*nu)
		M(2) 	= M(1)
		M(3) 	= M(1)
		M(4) 	= 0.
		M(5) 	= 0.
		M(6) 	= 0.
	Endif

	End Subroutine src2mom
