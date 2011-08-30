	program bid
	
	use nrtype
	use nr

	implicit none
	integer			:: nk
	real(sp)		:: kval
	real(sp), allocatable	:: kv(:), j1(:)

	kval = 8.

	do nk=89834,90000
		allocate(j1(nk), kv(nk))
		kv    = kval
		print *, 'before j1, with: ', nk
		j1 = bessj1(kv)
		print *, 'after j1'
		deallocate(j1, kv)
	end do
end program bid
