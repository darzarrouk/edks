subroutine write_disp(prefix, np, nrc, ux, uy, uz)

! prefix   prefix of data files
! nrc      number of recievers 
! np       number of patches 
! ux       m, displacements
! uy       m
! uz       m

implicit none

integer*4           :: id, irec, ip, np, nrc, nbyt
real*4              :: ux(np,nrc), uy(np,nrc), uz(np,nrc)
character (len=256) :: prefix, fname1, fname2, fname3

fname1  = trim(prefix)//'_ux.dis'
fname2  = trim(prefix)//'_uy.dis'
fname3  = trim(prefix)//'_uz.dis'

nbyt = nrc*4*np

open(20, file=fname1, form='unformatted', status='REPLACE', access='direct', recl=nbyt)
open(21, file=fname2, form='unformatted', status='REPLACE', access='direct', recl=nbyt)
open(22, file=fname3, form='unformatted', status='REPLACE', access='direct', recl=nbyt)

write(20,rec=1) ((ux(ip, irec), irec = 1, nrc), ip = 1, np)
write(21,rec=1) ((uy(ip, irec), irec = 1, nrc), ip = 1, np)
write(22,rec=1) ((uz(ip, irec), irec = 1, nrc), ip = 1, np)

close(20)
close(21)
close(22)

end subroutine write_disp
