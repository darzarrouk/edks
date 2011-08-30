	program bid
	
	use nrtype
	use nr

	implicit none
	integer			:: nk
	real(sp)		:: kval
	real(sp), allocatable	:: kv(:), j2(:)

	kval = 9.

	do nk=89835,90000
		allocate(j2(nk), kv(nk))
		kv    = kval
		print *, 'before j2, with: ', nk
		j2 = bessj(2,kv)
		print *, 'after j2'
		deallocate(j2, kv)
	end do
end program bid
