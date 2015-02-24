subroutine read_point_sub(prefix, ntsp, x, y, z, dV)

! prefix  prefix of data files
! ntsp    number of total sub points to read
! x       m, x-coordinate of point
! y       m, y-coordinate of point
! z       m, z-coordinate of point
! dV      m^3, amount of volume change (usually set to 1)

implicit none

real*4              :: x(ntsp), y(ntsp), z(ntsp), dV(ntsp)
integer*4           :: ntsp, ip, numvar, nbyt, funit
character (len=256) :: prefix, fname

funit = 20
numvar = 4
nbyt = numvar * ntsp * 4

fname = trim(prefix)//'.pat'

open(funit, file=fname, form='unformatted', access='direct', recl=nbyt, status='old')
read(funit,rec=1) (x(ip), y(ip), z(ip), dV(ip), ip=1,ntsp)
close(funit)

end subroutine read_point_sub
