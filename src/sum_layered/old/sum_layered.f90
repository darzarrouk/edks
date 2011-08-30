include 'edks.inc'

program sum_layered

use my_interfaces
use defs_module

implicit none

! Mandatory output arguments of mom2disp
real*4, allocatable :: uxs(:,:), uys(:,:), uzs(:,:)	! nsrc x nrec

! Mandatory output arguments of sum_layred
real*4, allocatable :: ux(:,:), uy(:,:), uz(:,:)	! npatch x nrec

! Mandatory input arguments of mom2disp
real*4, allocatable :: xr(:), yr(:)	        ! nrec
real*4, allocatable :: xs(:), ys(:), zs(:)       ! nsrc
real*4, allocatable :: xsp(:), ysp(:), zsp(:)   ! nspp
real*4, allocatable :: x(:), y(:), z(:)         ! npatch 
real*4, allocatable :: M(:,:)			! nsrc x 6
type    (edks_type) :: edks

! Optional input arguments of src3mom
real*4, allocatable :: strike(:), dip(:), strike_r(:), dip_r(:)
real*4, allocatable :: rake(:), W(:), L(:), slip(:)
real*4 :: dw, dy, Ms(6)
real*4 :: wg, yg, zg

integer*4, allocatable :: is(:)
integer*4 :: npw, npy, nsrc, nrec, np, nspp, irec, i, ip

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
call getarg(3,argum)
read(argum,*) nrec
call getarg(4,argum)
read(argum,*) np
call getarg(5,argum)
read(argum,*) npw
call getarg(6,argum)
read(argum,*) npy

nspp = npw*npy ! number of  sources per patch
nsrc = nspp*np ! total number of sources

allocate( is(nspp)      , M(nsrc,6) )
allocate( strike(np)    , dip(np)       , strike_r(np)  , dip_r(np) )
allocate( rake(np)      , W(np)         , L(np)         , slip(np)  )
allocate( x(np)         , y(np)         , z(np)                     )
allocate( xs(nsrc)      , ys(nsrc)      , zs(nsrc)                  )
allocate( xr(nrec)      , yr(nrec)                                  )
allocate( xsp(nspp)     , ysp(nspp)     , zsp(nspp)                 )
allocate( ux(np,nrec)   , uy(np,nrec)   , uz(np,nrec)               )
allocate( uxs(nsrc,nrec), uys(nsrc,nrec), uzs(nsrc,nrec)            )

! read in the edks (the greens functions)
call read_edks(elem_filename, edks)

! read in the receiver locations
call read_receivers(prefix, nrec, xr, yr)

! read in the patch information
call read_patch(prefix, np, x, y, z, strike, dip, rake, W, L, slip)

! need to loop over all patches, get local coordinates and fill single arrays

do ip = 1, np
   is = (/ ((ip-1)*nspp + i, i=1, nspp) /)
   call patch2pts(strike(ip), dip(ip), rake(ip), slip(ip), W(ip), L(ip), &
                  x(ip), y(ip), z(ip), npw, npy, xsp, ysp, zsp, Ms)  
!print *, "Ms in main", Ms
   xs(is)  = xsp
   ys(is)  = ysp
   zs(is)  = zsp
   M(is,:) = spread(Ms,1,nspp)
end do
print *, "Ms out of loop", Ms


call mom2disp(uxs, uys, uzs, xr, yr, xs, ys, zs, edks, M)

! need to sum pts back to patch contributions where arrays are indexed as
! ux(patch,receiver), uxs(source,receiver)
! could be done in terms of reshape and sum, but its too messy for me.

do ip = 1, np
   is = (/ ((ip-1)*nspp + i, i=1, nspp) /)
   do irec = 1, nrec
      ux(ip, irec) = sum(uxs(is,irec))
      uy(ip, irec) = sum(uys(is,irec))
      uz(ip, irec) = sum(uzs(is,irec))
   end do
end do


! write out displacement arrays 

call write_disp(prefix, np, nrec, ux, uy, uz)
   
end program sum_layered
