include 'edks.inc'

program sum_layered_volumetric

use my_interfaces
use defs_module

implicit none

real*4, allocatable :: uxs(:,:), uys(:,:), uzs(:,:) ! ntsp x nrec
real*4, allocatable :: ux(:,:), uy(:,:), uz(:,:)    ! npatch x nrec
real*4, allocatable :: xr(:), yr(:)             ! nrec
real*4, allocatable :: xsp(:), ysp(:), zsp(:)   ! nspp
real*4, allocatable :: x(:), y(:), z(:)         ! ntsp 

type (edks_type) :: edks

real*4, allocatable :: dV(:) 
real*4 :: dw, dy, M(6)
real*4 :: wg, yg, zg
real*4 :: nu, bulk, mu

integer*4 :: npw, npy, nrec, np, nspp, irec, i, ip, ir, ntsp

character (len=256) :: elem_filename, prefix
character (len=32)  :: argum
integer*4           :: iargc, narg

! ntsp is the number of total subsources (numpoints x numsubpatch, where numsubpatch is variable by patch)
! nspp is the maximum number of subsources per patch (this number is variable among the different patches, so we need the largest one).
       
! read in command line parameters
narg = iargc()
if(narg .ne. 6) then
   write(*,'(a)') ' '
   write(*,'(a)') 'Usage: sum_layered_volumetric edks_name geom_prefix #receivers #patches #ntsp #nspp'
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

allocate(dV(ntsp))
allocate(x(ntsp), y(ntsp), z(ntsp))
allocate(xsp(nspp), ysp(nspp), zsp(nspp))
allocate(ux(np,nrec), uy(np,nrec), uz(np,nrec))
allocate(uxs(nspp,nrec), uys(nspp,nrec), uzs(nspp,nrec))
allocate(xr(nrec), yr(nrec))

! read in the edks (the greens functions)
write(*,'(2a)') ' reading EDKS kernels from ', elem_filename
call read_edks(elem_filename, edks)

! read in the receiver locations
write(*,'(a)') ' reading the receiver locations'
call read_receivers(prefix, nrec, xr, yr)

! read in the point information
write(*,'(a)') ' reading the point information'
call read_point_sub(prefix, ntsp, x, y, z, dV)

! poisson's ratio
nu = 0.25
!bulk = 53.1d9
!mu = 33.075d9

! need to loop over all patches, get local coordinates and fill single arrays
write(*,'(a,i5,1x,a)') ' looping over the', np, 'points'
do ip = 1,np

   ! Update the moment tensor for this point
   M(1) = dV(ip) * 2.0 * (1.0 + nu) / (3.0 * (1.0 - 2.0 * nu))
   M(2) = M(1)
   M(3) = M(1)
   M(4) = 0.0
   M(5) = 0.0
   M(6) = 0.0   

   ! Calculating sub sources displacements for current point
   xsp(1) = x(ip)
   ysp(1) = y(ip)
   zsp(1) = z(ip)
   call mom2disp_sub(uxs, uys, uzs, xr, yr, xsp, ysp, zsp, edks, M, 1)

   ! need to sum pts back to patch contributions where arrays are 
   ! indexed as ux(patch,receiver), uxs(pt_source,receiver)
   ux(ip,:) = sum(uxs, dim=1)
   uy(ip,:) = sum(uys, dim=1)
   uz(ip,:) = sum(uzs, dim=1)

end do

! write out displacement arrays 
call write_disp(prefix, np, nrec, ux, uy, uz)

deallocate(dV)
deallocate(x, y, z)
deallocate(xsp, ysp, zsp)
deallocate(ux, uy, uz)
deallocate(uxs, uys, uzs)
deallocate(xr, yr)
   
end program sum_layered_volumetric
