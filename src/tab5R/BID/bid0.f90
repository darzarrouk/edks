	program bid
	
	use nrtype
	use nr

	implicit none
	integer			:: nk
	real(sp)		:: kval
	real(sp), allocatable	:: kv(:), j0(:)

	kval = 8.

	do nk=89834,90000
		allocate(j0(nk), kv(nk))
		kv    = kval
		print *, 'before j0, with: ', nk
		j0 = bessj0(kv)
		print *, 'after j0'
		deallocate(j0, kv)
	end do
end program bid
