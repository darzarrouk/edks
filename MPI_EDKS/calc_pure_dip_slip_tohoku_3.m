clear all;
close all;

km = 1.e3;
edks    = 'tohoku_3.edks'

%%%%% receivers info
nx    =    25;
ny    =    25;
xrmin =  -500*km;
yrmin = -1000*km;
xrmax =   500*km;
yrmax =     0*km;
xr = linspace(xrmin,xrmax,nx);
yr = linspace(yrmin,yrmax,ny);
[xrg, yrg] = meshgrid(xr,yr);

%%%%% source info

npw = 1;  % points per patch in the down dip direction
npy = 1;  % points per patch in the along strike direction

labelmec = 'Pure_strike_slip';

%%% Sources  %%%
%MM     = load('fault_geometry_my_format');
%th     = load('slip_model_my_format');
xs     = [0. 0.]*km;
ys     = [-500. -500.]*km;
zs     = [30. 30.]*km;
strike = [0. 0.];
dip    = [90. 90.];
area   = [1. 1.]*km*km;
np     = length(xs);
st_sl  = [2.e4 2.e4];
di_sl  = [0. 0.];
rake   = atan2(di_sl,st_sl)*180./pi;
slip   = sqrt(st_sl.^2 + di_sl.^2);
L      = sqrt(area);
W      = sqrt(area);
label = sprintf('TH3_%s', labelmec);

%slip = slip.*area;
%L    = ones(length(L),1); W=L;


[uxt, uyt, uzt] = layered_disloc(xs, ys, zs, strike, dip, rake, slip, L, W, npw, npy, xrg(:), yrg(:), edks);

print uxt
%if(np > 1)
%    EXT = reshape(sum(uxt'),ny,nx);
%    EYT = reshape(sum(uyt'),ny,nx);
%    EZT = reshape(sum(uzt'),ny,nx);
%else
%    EXT = reshape(uxt,ny,nx);
%    EYT = reshape(uyt,ny,nx);
%    EZT = reshape(uzt,ny,nx);  
%end

% Output
xfile = [label, '_Ux'];
yfile = [label, '_Uy'];
zfile = [label, '_Uz'];

fhx=fopen(xfile,'wt');
fhy=fopen(yfile,'wt');
fhz=fopen(zfile,'wt');
fprintf(fhx, '%6d %6d %12.6f %12.6f %12.6f %12.6f\n', ... 
        nx, ny, xrmin, xrmax,yrmin, yrmax);  
fprintf(fhy, '%6d %6d %12.6f %12.6f %12.6f %12.6f\n', ...
        nx, ny, xrmin, xrmax,yrmin, yrmax);
fprintf(fhz, '%6d %6d %12.6f %12.6f %12.6f %12.6f\n', ...
        nx, ny, xrmin, xrmax,yrmin, yrmax);

for jy=1:ny
    fprintf(fhx, '%15.6e', EXT(jy,:)); fprintf(fhx, '\n');
    fprintf(fhy, '%15.6e', EYT(jy,:)); fprintf(fhy, '\n');
    fprintf(fhz, '%15.6e', EZT(jy,:)); fprintf(fhz, '\n');
end
fclose(fhx);
fclose(fhy);
fclose(fhz);

fp = fopen([label, 'x_y_Ux_Uy_Uz.dat'],'wt');
for jx=1:nx
for jy=1:ny
    fprintf(fp, '%12.2f %12.2f %12.6f %12.6f %12.6f\n', ...
            xr(jx)/km, yr(jy)/km, EXT(jy,jx), EYT(jy,jx), EZT(jy,jx));
end
end
fclose(fp)

