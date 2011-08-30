	program f2
        integer i,j


        open(unit=7, file='sequen', form='unformatted')
        do i=1,20
                read(7, end=10) j
                write(*,*) j
        end do
        close(7)

 10     continue

        open(unit=8, file='direct', form='unformatted', recl=4, access='direct')
        do i=1,20
                read(8,rec=i,end=20) j
                write(*,*) j
        end do
        close(8)

 20     continue
        end

