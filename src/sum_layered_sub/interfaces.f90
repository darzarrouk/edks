	include	'edks.inc'
	module  my_interfaces

	interface
        Subroutine src2mom                                      &
                 ( M, o_Mo, o_MomT, o_nu, o_mu,                 &
                 o_Area, o_strike, o_dip, o_rake, o_St, o_Sp, o_deltaV  )

        Use defs_Module
!       use f90_unix            !needed on the NAG compiler
        Implicit None
        Optional o_Mo, o_MomT, o_nu, o_mu,                              &
                 o_Area, o_strike, o_dip, o_rake, o_St, o_Sp, o_deltaV

!       Mandatory output arguments
        Real, Intent(out)               :: M(6)

!       Optional input arguments
        Real, Intent(in)                :: o_Mo, o_MomT(1:6), o_mu, o_nu
        Real, Intent(in)                :: o_strike, o_dip
        Real, Intent(inout)             :: o_rake
        Real, Intent(inout)             :: o_St, o_Sp
        Real, Intent(in)                :: o_Area, o_DeltaV
	end subroutine src2mom
	end interface

	interface
        Subroutine mom2disp_sub(Ux, Uy, Uz, Xr, Yr, Xs, Ys, Zs, edks, M, nsrc)
        !Subroutine mom2disp(Ux, Uy, Uz, Xr, Yr, Xs, Ys, Zs, edks, M)

        Use defs_Module
!       use f90_unix            !needed on the NAG compiler
        Implicit None

!       Mandatory output arguments
        Real, Intent(out)               :: Ux(:,:), Uy(:,:), Uz(:,:)    ! nsrc x nrec

!       Mandatory input arguments
        Real, Intent(in), target        :: Xr(:), Yr(:)                 ! nrec
        Real, Intent(in), target        :: Xs(:), Ys(:), Zs(:)          ! nsrc
        Real, Intent(in), target        :: M(6)                         ! 6
        Integer, Intent(in), target        :: nsrc
        Type (edks_Type), Intent(in)    :: edks

!        end subroutine mom2disp
	end subroutine mom2disp_sub
	end interface

end module my_interfaces
