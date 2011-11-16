module defs_module

type edks_type
	character (len=128)	prefix
	character (len= 32)	date
	character (len= 32)	version
	character (len=128)	comment
	integer			nlayer
	real,pointer	::	rho(:)
	real,pointer	::	alpha(:)
	real,pointer	::	beta(:)
	real,pointer	::	thickness(:)
	integer			ndepth,   ndista
	real			depthmin, depthmax
	real			distamin, distamax
	real,pointer	::	depths(:)
	real,pointer	::	distas(:)
	real,pointer	::	zrtdsx(:,:,:)
end type edks_type		

end module defs_module
