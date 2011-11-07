subroutine patch2pts_sub(ntsp, nspp, ip, x, y, z, strike, dip, rake, area, slip, ident, xs, ys, zs, M, nsrc) 

use my_interfaces
use defs_module

implicit none

integer*4 :: ipt, ntsp, nspp, ip, isub, nsrc

real*4 :: strike(ntsp), dip(ntsp), rake(ntsp)
real*4 :: x(ntsp), y(ntsp), z(ntsp), slip(ntsp), M(6)
real*4 :: xs(nspp), ys(nspp), zs(nspp), area(ntsp)

real*4 :: strikes, dips, rakes, areas, slips

integer*4 :: ident(ntsp)


! needs to parse original monster file to extract the subtriangles for a given master triangle
xs = 0.0
ys = 0.0
zs = 0.0
ipt = 0
do isub = 1, ntsp
    if (ident(isub) == ip) then
        ipt = ipt+1
        xs(ipt) = x(isub)
        ys(ipt) = y(isub)
        zs(ipt) = z(isub)
        strikes = x(isub)
        dips    = dip(isub)
        rakes   = rake(isub)
        areas   = area(isub)
        slips   = slip(isub)
    endif
end do

nsrc = ipt

M = 0

! convert patch mechanism to moment for use in du_layer
call src2mom(M, o_Area=Areas, o_strike=strikes, o_dip=dips, o_rake=rakes, o_St=slips)

end subroutine patch2pts_sub
