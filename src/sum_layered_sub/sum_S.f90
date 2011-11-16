include 'edks.inc'

program sum_S

use my_interfaces
use defs_module

implicit none

real*4, allocatable :: uxs(:,:), uys(:,:), uzs(:,:)	! nsrc x nrec
real*4, allocatable :: ux(:,:), uy(:,:), uz(:,:)	! npatch x nrec
real*4, allocatable :: xr(:), yr(:)	                ! nrec
real*4, allocatable :: xsp(:), ysp(:), zsp(:)           ! nspp
real*4, allocatable :: x(:), y(:), z(:)                 ! npatch 
type    (edks_type) :: edks

real*4, allocatable :: strike(:), dip(:), strike_r(:), dip_r(:) 
real*4, allocatable :: rake(:), W(:), L(:), slip(:)
real*4 :: dw, dy, M(6)
real*4 :: wg, yg, zg

integer*4 :: npw, npy, nrec, np, nspp, irec, i, ip

character (len=256) :: elem_filename, prefix
character (len=32)  :: argum
integer*4           :: iargc, narg
       
! read in command line parameters
narg = iargc()
if(narg .ne. 6) then
   write(*,'(a)') ' '
   write(*,'(a)') 'Usage: sum_layered edks_name geom_prefix #receivers #patches NPW NPY '
   write(*,'(a)') ' '
   call exit(1)
endif
call getarg(1,elem_filename)
call getarg(2,prefix)
!call getarg(3,argum)
!read(argum,*) nrec
!call getarg(4,argum)
!read(argum,*) np
!call getarg(5,argum)
!read(argum,*) npw
!call getarg(6,argum)
!read(argum,*) npy

np  = 1
npy = 1
npw = 1
nspp = npw*npy ! number of  sources per patch

! read in the edks (the greens functions)
call read_edks(elem_filename, edks)

! read in the receiver locations
call read_receivers(prefix, nrec, xr, yr)
latmin = -20.
latmax =  20.
nlat   = 100
lonmin = -20.
lonmax =  20.
nlon   = 100
dlon   = (lonmax-lonmin)/(nlon-1)
dlat   = (latmax-latmin)/(nlat-1)
nr     = nlon*nlat

allocate(xr(nr), yr(nr))
do jlon = 1,nlon
do jlat = 1,nlat
    jr     = (jlon-1)*nlon+jlat
    xr(jr) = lonmin + dlon*(jlon-1)
    yr(jr) = latmin + dlat*(jlat-1)
end do
end do

allocate( strike(np)    , dip(np)       , strike_r(np)  , dip_r(np) )
allocate( rake(np)      , W(np)         , L(np)         , slip(np)  )
allocate( x(np)         , y(np)         , z(np)                     )
allocate( xsp(nspp)     , ysp(nspp)     , zsp(nspp)                 )
allocate( ux(np,nrec)   , uy(np,nrec)   , uz(np,nrec)               )
allocate( uxs(nspp,nrec), uys(nspp,nrec), uzs(nspp,nrec)            )
strike(1) = 45.
dip(1)    = 45.
rake(1)   = 45.
slip(1)   = 1.
W(1)      = 1.
L(1)      = 1.
x(1)      = 0.
y(1)      = 0.
z(1)      = 0.

! read in the patch information
!call read_patch(prefix, np, x, y, z, strike, dip, rake, W, L, slip)

! need to loop over all patches, get local coordinates and fill single arrays

do ip = 1, np

   call patch2pts(strike(ip), dip(ip), rake(ip), slip(ip), W(ip), L(ip), &
                  x(ip), y(ip), z(ip), npw, npy, xsp, ysp, zsp, M)  

   call mom2disp(uxs, uys, uzs, xr, yr, xsp, ysp, zsp, edks, M)

   ! need to sum pts back to patch contributions where arrays are 
   ! indexed as ux(patch,receiver), uxs(pt_source,receiver)

   ux(ip, :) = sum(uxs, dim=1)
   uy(ip, :) = sum(uys, dim=1)
   uz(ip, :) = sum(uzs, dim=1)

end do

! write out displacement arrays 

call write_disp(prefix, np, nrec, ux, uy, uz)
   
end program sum_S
