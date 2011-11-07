include 'edks.inc'

program sum_layered_sub

use my_interfaces
use defs_module

implicit none

real*4, allocatable :: uxs(:,:), uys(:,:), uzs(:,:)	! nsrc x nrec
real*4, allocatable :: ux(:,:), uy(:,:), uz(:,:)	! npatch x nrec
real*4, allocatable :: xr(:), yr(:)	        ! nrec
real*4, allocatable :: xsp(:), ysp(:), zsp(:)   ! nspp
real*4, allocatable :: x(:), y(:), z(:)         ! ntsp 

integer*4, allocatable :: ident(:)         ! ntsp 
type    (edks_type) :: edks

real*4, allocatable :: strike(:), dip(:), area(:) 
real*4, allocatable :: rake(:), slip(:)
real*4 :: dw, dy, M(6)
real*4 :: wg, yg, zg

integer*4 :: npw, npy, nrec, np, nspp, irec, i, ip, ir, nsrc, ntsp

character (len=256) :: elem_filename, prefix
character (len=32)  :: argum
integer*4           :: iargc, narg

! ntsp is the number of total subsources (numpatch x numsubpatch, where numsubpatch is variable by patch)
! nspp is the maximum number of subsources per patch (this number is variable among the different patches, so we need the largest one).
       
! read in command line parameters
narg = iargc()
if(narg .ne. 6) then
   write(*,'(a)') ' '
   write(*,'(a)') 'Usage: sum_layered_sub edks_name geom_prefix #receivers #patches #ntsp #nspp'
   write(*,'(a)') ' '
   write(*,'(a)') '   - #receivers is the number of receivers'
   write(*,'(a)') '   - #patches is the number of (finite) sources to model'
   write(*,'(a)') '   - #ntsp is the number of total subsources (numpatch x numsubpatch,'
   write(*,'(a)') '     where numsubpatch is variable by patch)'
   write(*,'(a)') '   - #nspp is the maximum number of subsources per patch (this number'
   write(*,'(a)') '     is variable among the different patches, so we need the largest one).'
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
read(argum,*) ntsp
call getarg(6,argum)
read(argum,*) nspp

allocate( strike(ntsp)  , dip(ntsp)      )
allocate( rake(ntsp)    , slip(ntsp)     , area(ntsp)  )
allocate( x(ntsp)       , y(ntsp)        , z(ntsp)                     )
allocate( ident(ntsp)                                  )
allocate( xsp(nspp)     , ysp(nspp)      , zsp(nspp)                 )
allocate( ux(np,nrec)   , uy(np,nrec)    , uz(np,nrec)               )
allocate( uxs(nspp,nrec), uys(nspp,nrec) , uzs(nspp,nrec)            )
allocate( xr(nrec)      , yr(nrec)                                  )

! read in the edks (the greens functions)
write(*,'(a)') ' reading EDKS kernels'
call read_edks(elem_filename, edks)

! read in the receiver locations
write(*,'(a)') ' reading the receiver locations'
call read_receivers(prefix, nrec, xr, yr)
! read in the patch information
write(*,'(a)') ' reading the patch information'
call read_patch_sub(prefix, ntsp, x, y, z, strike, dip, rake, area, slip, ident)
! need to loop over all patches, get local coordinates and fill single arrays
do ip = 1, np
   print *, ' selecting sub sources for patch', ip
   call patch2pts_sub(ntsp, nspp, ip, x, y, z, strike, dip, rake, area, slip, ident, xsp, ysp, zsp, M, nsrc)

   print *, ' Calculating ',nsrc, ' sub sources displacements for patch', ip
   call mom2disp_sub(uxs, uys, uzs, xr, yr, xsp, ysp, zsp, edks, M, nsrc)

   ! need to sum pts back to patch contributions where arrays are 
   ! indexed as ux(patch,receiver), uxs(pt_source,receiver)
   ux(ip, :) = sum(uxs, dim=1)
   uy(ip, :) = sum(uys, dim=1)
   uz(ip, :) = sum(uzs, dim=1)

end do

! write out displacement arrays 

call write_disp(prefix, np, nrec, ux, uy, uz)
   
end program sum_layered_sub
