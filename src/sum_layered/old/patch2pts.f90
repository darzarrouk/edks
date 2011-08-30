subroutine patch2pts(strike, dip, rake, slip, W, L, &
                     x, y, z, npw, npy, xs, ys, zs, Ms) 

use my_interfaces
use defs_module

implicit none

real*4 :: strike, strike_rad, dip, dip_rad, W, L, rake
real*4 :: x, y, z, slip, Ms(6)
real*4 :: xso, yso, zso
real*4 :: xs(*), ys(*), zs(*)
real*4 :: dw, dy, area
real*4 :: coss, sins, cosd, sind, pi
real*4 :: wax(npw), yax(npy)

integer*4 :: npw, npy, iw, iy, ipt

pi         = 3.14159265358979
strike_rad = strike * pi/180.0 + pi/2 !convert to radians
dip_rad    = dip * pi/180.0           !convert to radians

coss     = cos(strike_rad)
sins     = sin(strike_rad)
cosd     = cos(dip_rad)
sind     = sin(dip_rad) 

Ms = 0
dw   = W/npw
dy   = L/npy
area = W*L/npw/npy !area associated with each point source

! pts will be in center of equally sizes cells

wax = (/ (iw*dw-dw/2, iw=1,npw) /)
yax = (/ (iy*dy-dy/2-L/2, iy=1,npy) /)  !shift origin to center of top 

! define points in local patch coordinates (with z=0 at top of patch)
! then rotate and translate to real coordinates

do iy = 1, npy
   do iw = 1, npw

      yso = yax(iy)
      xso = wax(iw)*cosd
      zso = wax(iw)*sind

      ipt = (iw-1)*npw+iy
      xs(ipt) = xso*sins + yso*coss + x
      ys(ipt) = xso*coss - yso*sins + y 
      zs(ipt) = zso                 + z

   end do 
end do

! convert patch mechanism to moment for use in du_layer
call src2mom(Ms, o_Area=Area, o_strike=strike, o_dip=dip, o_rake=rake, o_St=slip)

!print *, "Ms in patch", Ms

end subroutine patch2pts
