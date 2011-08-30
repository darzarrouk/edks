MODULE model
	integer, parameter	:: MXMD=10
	integer			:: N, Np, M
	real, 	allocatable	:: vp0(:), vs0(:), rh0(:), th0(:)
	real, 	allocatable	::  vp(:),  vs(:),  rh(:),  th(:)
	real, 	allocatable	:: lambda(:), mu(:)
	real, 	allocatable	:: delta(:), gamma(:)
END MODULE model
	
MODULE parameters
	character (len=256)	:: model_filename
	real			:: htol
END MODULE parameters

MODULE  wavenumber
	real			:: facmin, facmax, HH
	real, 	allocatable	:: kv(:)
	real			:: kmin, dk, kmax
	integer			:: nk
END MODULE wavenumber

MODULE constants
	real, parameter		:: Id2(2,2) = reshape(source = (/ 1., 0., 0., 1. /), &
							shape = (/ 2, 2 /))
	real, parameter		:: Ze2(2,2) = reshape(source = (/ 0., 0., 0., 0. /), &
							shape = (/ 2, 2 /))
	real, parameter		:: DET_MIN  = 1.e-20
	real, parameter		:: M_PI = 3.14159265
END MODULE constants
