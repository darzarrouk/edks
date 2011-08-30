! Include	'edks.inc'
Subroutine mom2disp(Ux, Uy, Uz, Xr, Yr, Xs, Ys, Zs, edks, MM)
Use defs_Module 

! use f90_unix    	!needed on the NAG compiler
Implicit None

! Mandatory output arguments
Real, Intent(out)		:: Ux(:,:), Uy(:,:), Uz(:,:)	! nsrc x nrec

! Mandatory input arguments
Real, Intent(in), target       	:: Xr(:), Yr(:)			! nrec
Real, Intent(in), target        :: Xs(:), Ys(:), Zs(:)		! nsrc
Real, Intent(in), target       	:: MM(:,:)			! nsrc x 6
Type (edks_Type), Intent(in)	:: edks

! This subroutine computes the displacements produced by a horizontal
! set of point sources in a given set of recievers.
! It uses a pre-calculated table containing the elementary displacements
! corresponding to the structure and source's depth choosen.
!
! INPUT PARAMETERS:
!
! EDKS: 		name of the array containig the precalculated elementary 
! 		displacements depending on the structure sources's depths
!
! MECHANISM: 	The  source mechanism should be given as a moment tensor M(:,6)
!
! SOURCE GEOMETRY
!
! nsrc:		nb of point sources to be considered: x east, y north, z DOWN !!
! Xs(nsrc), Xs(nsrc), Zs(nsrc): x,y  and z coordinates of the point sources
!
! RECIEVERS
!
! nrec		Number of recievers
! Xr(nrec), Yr(nrec): x (east) and y (north) coordinates of the recievers
!
! OUTPUT PARAMETERS
! Ux(nsrc,nrec), Uy(nsrc,nrec), Uz(nsrc,nrec)   (east, north, UP !!!)
!
! COORDINATES
!
! distances and displacements are given in meters.
! x grows towards the east 	Xs, Xr, Ux
! y grows towards the north	Ys, Yr, Uy
! z grows DOWNWARDS		Zs
! z grows UPWARDS		        Uz
!

Integer			:: nsrc
Integer			:: nrec

! parameters

Real,	Parameter	:: M_PI 	= 3.14159265358979 
Real,	Parameter	:: EPSILON 	= 1e-10

! local variables

Real	 		:: zrtdsx(10),   zrtdsxd1(10),    zrtdsxd2(10)
Real	 		:: zrtdsxh1d1(10), zrtdsxh1d2(10)
Real	 		:: zrtdsxh2d1(10), zrtdsxh2d2(10)
Real			:: ZSS, RSS, TSS, ZDS, RDS, TDS, ZDD, RDD, ZEX, REX
Equivalence		(zrtdsx( 1), ZDD)
Equivalence		(zrtdsx( 2), ZDS)
Equivalence		(zrtdsx( 3), ZSS)
Equivalence		(zrtdsx( 4), RDD)
Equivalence		(zrtdsx( 5), RDS)
Equivalence		(zrtdsx( 6), RSS)
Equivalence		(zrtdsx( 7), TDS)
Equivalence		(zrtdsx( 8), TSS)
Equivalence		(zrtdsx( 9), ZEX)
Equivalence		(zrtdsx(10), REX)

Integer			:: jhbest, jhseco
Integer			:: jdbest, jdseco
Real			:: ratioh, ratiod
Real			:: x, y
Real			:: ws, qr, vt

! Source
Real			:: M(6)        
Real			:: x0, y0, z0
Integer			:: jsrc, jrec, j

! Others
Real			:: r, saz, caz, s2az, c2az

!Counters
Integer                 :: ix, iy

! loop on the sources
nsrc = Size(Xs)
nrec = Size(Xr)

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(MM, Ux, Uy, Uz, Xs, Ys, Zs, Xr, Yr, edks, nrec, nsrc)
!$OMP DO 
Do j = 1,nrec*nsrc
   jrec = int((j-1)/nsrc+1.0)
   jsrc = j-(jrec-1)*nsrc
   
   M  = MM(jsrc,:)
   x0 = Xs(jsrc)
   y0 = Ys(jsrc)
   z0 = Zs(jsrc)

   Call best2(edks % depths, edks % ndepth, z0, jhbest, jhseco)
   ratioh 	= (z0 - edks % depths(jhbest))/(edks % depths(jhseco) - edks % depths(jhbest))
   
   x = Xr(jrec)
   y = Yr(jrec)
   
   ! initialize the displacements
   Ux(jsrc, jrec) = 0.
   Uy(jsrc, jrec) = 0.
   Uz(jsrc, jrec) = 0.
   
   ! geometry source-receiver
   r  = Sqrt((x-x0)**2+(y-y0)**2)
   caz = 1.
   saz = 0.
   If( r .Gt. EPSILON) Then
      caz = (y-y0)/r
      saz = (x-x0)/r
   endif
   c2az = 2.*caz*caz - 1.
   s2az = 2.*saz*caz
   
   ! interpolation in depths and distances
   
   Call best2(edks % distas, edks % ndista, r, jdbest, jdseco)
   ratiod 	= (r - edks % distas(jdbest))/(edks % distas(jdseco) - edks % distas(jdbest))
   
   zrtdsxh1d1 	= edks % zrtdsx(jhbest,jdbest,1:10)
   zrtdsxh1d2 	= edks % zrtdsx(jhbest,jdseco,1:10)
   zrtdsxh2d1 	= edks % zrtdsx(jhseco,jdbest,1:10)
   zrtdsxh2d2 	= edks % zrtdsx(jhseco,jdseco,1:10)
   zrtdsxd1   	= zrtdsxh1d1 + ratioh * ( zrtdsxh2d1 -  zrtdsxh1d1 )
   zrtdsxd2   	= zrtdsxh1d2 + ratioh * ( zrtdsxh2d2 -  zrtdsxh1d2 )
   zrtdsx       = zrtdsxd1   + ratiod * (  zrtdsxd2  -   zrtdsxd1  )
   
   ! reconstruction of the actual displacement Xiaobi vs  Herrman
   ! The coefficients created by tab5 are in Xiaobi's notation
   ! The equations just below follow Herrman's notation; hence

   ZDS=-ZDS
   RDS=-RDS
   TSS=-TSS
   ! Vertical component  (positive down)
   ws =   M(2)*( ZSS*c2az/2. - ZDD/6. + ZEX/3.) &
        + M(3)*(-ZSS*c2az/2. - ZDD/6. + ZEX/3.) &
        + M(1)*( ZDD + ZEX)/3.  &
        + M(6)*  ZSS*s2az &
        + M(4)*  ZDS*caz &
        + M(5)*  ZDS*saz
   
   ! Radial component    (positive away from the source)
   qr = M(2)*( RSS*c2az/2. - RDD/6. + REX/3.) &
        + M(3)*(-RSS*c2az/2. - RDD/6. + REX/3.) &
        + M(1)*( RDD + REX)/3.  &
        + M(6)*  RSS*s2az &
        + M(4)*  RDS*caz &
        + M(5)*  RDS*saz 
   
   ! Tangential component (positive if clockwise from zenithal view)
   vt =   M(2)*TSS*s2az/2. &
        - M(3)*TSS*s2az/2. &
        - M(6)*TSS*c2az &
        + M(4)*TDS*saz &
        - M(5)*TDS*caz 

   ! Cartesian components (x:east, y=north, z=up)
   Ux(jsrc,jrec) =  qr*saz + vt*caz
   Uy(jsrc,jrec) =  qr*caz - vt*saz
   Uz(jsrc,jrec) = -ws
enddo
!$OMP END DO  
!$OMP END PARALLEL
     
Return
end subroutine mom2disp

Subroutine best2(vect, nv, val, jbest, jseco)

Integer		:: nv
Real		:: vect(nv)
Real		:: val

! local variables
Integer 	:: jh, jbest, jseco, jswap
Real		:: dhbest, dhseco, dhswap

jbest = 1;
jseco = 2;
dhbest = Abs(val - vect(jbest))
dhseco = Abs(val - vect(jseco))
Do jh = 3, nv
   If(dhbest > dhseco) Then
      dhswap = dhseco
      dhseco = dhbest
      dhbest = dhswap
      
      jswap  = jseco
      jseco  = jbest
      jbest  = jswap
   endif
   
   dhswap = Abs(val - vect(jh))
   If(dhswap < dhseco) Then
      dhseco = dhswap
      jseco  = jh
   endif
end Do

If(dhbest > dhseco) Then
   dhswap = dhseco
   dhseco = dhbest
   dhbest = dhswap
   
   jswap  = jseco
   jseco  = jbest
   jbest  = jswap
endif

Return
end Subroutine best2
