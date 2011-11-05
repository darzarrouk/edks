subroutine read_patch_sub(prefix, ntsp, x, y, z, strike, dip, &
                      rake, area, slip, ident)

! prefix  prefix of data files
! ntsp      number of total sub patches to read
! x       m, middle of each patch 
! y       m, middle of each patch
! z       m, middle of each patch, positive down 
! strike  deg
! dip     deg
! rake    deg
! area    m**2
! slip    m, amount of slip (usually set to 1.0)
! ident   integer index to master triangle

implicit none

real*4              :: x(ntsp), y(ntsp), z(ntsp), strike(ntsp), dip(ntsp), rake(ntsp)
real*4              :: slip(ntsp), area(ntsp)
integer*4           :: ntsp, ip, numvar, nbyt, funit, ident(ntsp)
character (len=256) :: prefix, fname

funit = 20
numvar = 9
nbyt = numvar * ntsp * 4

fname = trim(prefix)//'.pat'

open(funit, file=fname, form='unformatted', access='direct', recl=nbyt, status='old')

read(funit,rec=1) (x(ip), y(ip), z(ip), strike(ip), dip(ip), &
                   rake(ip), area(ip), slip(ip), ident(ip), ip=1, ntsp)

close(funit)

end subroutine read_patch_sub
