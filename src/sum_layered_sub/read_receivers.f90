subroutine read_receivers(prefix, nrec, xr, yr)

! prefix   prefix of data files
! nrec     number of recievers to read
! xr       m
! yr       m

implicit none

real*4              :: xr(nrec), yr(nrec)
integer*4           :: nrec, ir, numvar, funit, nbyt
character (len=256) :: prefix, fname

fname = trim(prefix)//'.rec'

funit = 20

numvar = 2
nbyt   = numvar*4*nrec

open(funit, file=fname, form='unformatted', access='direct', recl=nbyt, status='old')

read(funit,rec=1) (xr(ir), yr(ir), ir = 1, nrec)

close(funit)

end subroutine read_receivers
