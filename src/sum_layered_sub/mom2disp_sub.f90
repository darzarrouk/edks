! Include	'edks.inc'
Subroutine mom2disp_sub(Ux, Uy, Uz, Xr, Yr, Xs, Ys, Zs, edks, M, nsrc)
Use defs_Module 

! use f90_unix    	!needed on the NAG compiler
Implicit None

! Mandatory output arguments
Real, Intent(out)		:: Ux(:,:), Uy(:,:), Uz(:,:)	! nsrc x nrec

! Mandatory input arguments
Real, Intent(in), target       	:: Xr(:), Yr(:)			! nrec
Real, Intent(in), target        :: Xs(:), Ys(:), Zs(:)		! nsrc
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
! MECHANISM: 	The  source mechanism is assumed fixed for all sources
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
! z grows DOWNWAzrtdsx(5)		Zs
! z grows UPWAzrtdsx(5)		        Uz
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

! loop on the sources
nrec = Size(Xr)

! initialize the displacements
Ux = 0.0
Uy = 0.0
Uz = 0.0
   
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(M, Ux, Uy, Uz, Xs, Ys, Zs, Xr, Yr, edks, nrec, nsrc)
!$OMP DO 
Do j = 1,nrec*nsrc
   jrec = int((j-1)/nsrc+1.0)
   jsrc = j-(jrec-1)*nsrc
   x0 = Xs(jsrc)
   y0 = Ys(jsrc)
   z0 = Zs(jsrc)

   Call best2(edks % depths, edks % ndepth, z0, jhbest, jhseco)
   ratioh 	= (z0 - edks % depths(jhbest))/(edks % depths(jhseco) - edks % depths(jhbest))
   
   x = Xr(jrec)
   y = Yr(jrec)
   
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
   
   Call best2a(edks % distas, edks % ndista, r, jdbest, jdseco)
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

   zrtdsx(2)=-zrtdsx(2)
   zrtdsx(5)=-zrtdsx(5)
   zrtdsx(8)=-zrtdsx(8)
   ! Vertical component  (positive down)
   ws =   M(2)*( zrtdsx(3)*c2az/2. - zrtdsx(1)/6. + zrtdsx(9)/3.) &
        + M(3)*(-zrtdsx(3)*c2az/2. - zrtdsx(1)/6. + zrtdsx(9)/3.) &
        + M(1)*( zrtdsx(1) + zrtdsx(9))/3.  &
        + M(6)*  zrtdsx(3)*s2az &
        + M(4)*  zrtdsx(2)*caz &
        + M(5)*  zrtdsx(2)*saz
   
   ! Radial component    (positive away from the source)
   qr = M(2)*( zrtdsx(6)*c2az/2. - zrtdsx(4)/6. + zrtdsx(10)/3.) &
        + M(3)*(-zrtdsx(6)*c2az/2. - zrtdsx(4)/6. + zrtdsx(10)/3.) &
        + M(1)*( zrtdsx(4) + zrtdsx(10))/3.  &
        + M(6)*  zrtdsx(6)*s2az &
        + M(4)*  zrtdsx(5)*caz &
        + M(5)*  zrtdsx(5)*saz 
   
   ! Tangential component (positive if clockwise from zenithal view)
   vt =   M(2)*zrtdsx(8)*s2az/2. &
        - M(3)*zrtdsx(8)*s2az/2. &
        - M(6)*zrtdsx(8)*c2az &
        + M(4)*zrtdsx(7)*saz &
        - M(5)*zrtdsx(7)*caz 

   ! Cartesian components (x:east, y=north, z=up)
   Ux(jsrc,jrec) =  qr*saz + vt*caz
   Uy(jsrc,jrec) =  qr*caz - vt*saz
   Uz(jsrc,jrec) = -ws
enddo
!$OMP END DO  
!$OMP END PARALLEL
     
Return
end subroutine mom2disp_sub

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


Subroutine best2a(vect, nv, val, jbest, jseco)

Integer		:: nv
Real		:: vect(nv)
Real		:: val

!	local variables

Integer 	:: jh, jbest, jseco, jswap, jstart, jend
Real		:: dhbest, dhseco, dhswap
Real		:: step

step = vect(2)-vect(1)
jstart = val/step - 3
jstart = max(jstart, 1)
jstart = min(nv-1, jstart)

jend   = jstart + 6
jend   = min(jend, nv)

jbest = jstart;
jseco = jbest + 1;

dhbest = Abs(val - vect(jbest))
dhseco = Abs(val - vect(jseco))

Do jh = jstart+2, jend
   If(dhbest > dhseco) Then
      dhswap = dhseco
      dhseco = dhbest
      dhbest = dhswap
      
      jswap  = jseco
      jseco  = jbest
      jbest  = jswap
   Endif
   
   dhswap = Abs(val - vect(jh))
   
   If(dhswap < dhseco) Then
      dhseco = dhswap
      jseco  = jh
      
   Endif
End Do

If(dhbest > dhseco) Then
   dhswap = dhseco
   dhseco = dhbest
   dhbest = dhswap
   jswap  = jseco
   jseco  = jbest
   jbest  = jswap
Endif

Return

End Subroutine best2a

