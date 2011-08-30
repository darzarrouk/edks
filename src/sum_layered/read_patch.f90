subroutine read_patch(prefix, np, x, y, z, strike, dip, &
                      rake, width, length, slip)

! prefix  prefix of data files
! np      number of fault patches to read
! x       m, middle of top edge of fault
! y       m, middle of top edge of fault
! z       m, depth of middle of top edge of fault, positive down 
! strike  deg
! dip     deg
! rake    deg
! width   m, down dip extent
! length  m, along strike extent
! slip    m, amount of slip (usually set to 1.0)

implicit none

real*4              :: x(np), y(np), z(np), strike(np), dip(np), rake(np)
real*4              :: slip(np), width(np), length(np)
integer*4           :: np, ip, numvar, nbyt, funit
character (len=256) :: prefix, fname

funit = 20
numvar = 9
nbyt = numvar * np * 4

fname = trim(prefix)//'.pat'

open(funit, file=fname, form='unformatted', access='direct', recl=nbyt, status='old')

read(funit,rec=1) (x(ip), y(ip), z(ip), strike(ip), dip(ip), &
                   rake(ip), width(ip), length(ip), slip(ip), ip=1, np)

close(funit)

end subroutine read_patch
