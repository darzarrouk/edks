clear all;
close all;

km = 1.e3;
edks    = 'tohoku_3.edks'

%%%%% receivers info
nx    =    100;
ny    =    100;
xrmin =  -200*km;
yrmin = -200*km;
xrmax =   200*km;
yrmax =  200*km;
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
xs     = [0.]*km;
ys     = [0.]*km;
zs     = [30.]*km;
strike = [0.];
dip    = [90.];
area   = [1.]*km*km;
np     = length(xs);
st_sl  = [1];
di_sl  = [0.];
rake   = atan2(di_sl,st_sl)*180./pi;
slip   = sqrt(st_sl.^2 + di_sl.^2);
L      = sqrt(area);
W      = sqrt(area);
label = sprintf('TH3_%s', labelmec);

%slip = slip.*area;
%L    = ones(length(L),1); W=L;


[uxt, uyt, uzt] = layered_disloc(xs, ys, zs, strike, dip, rake, slip, L, W, npw, npy, xrg(:), yrg(:), edks);

if(np > 1)
    EXT = reshape(sum(uxt'),ny,nx);
    EYT = reshape(sum(uyt'),ny,nx);
    EZT = reshape(sum(uzt'),ny,nx);
else
    EXT = reshape(uxt,ny,nx);
    EYT = reshape(uyt,ny,nx);
    EZT = reshape(uzt,ny,nx);  
end

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures
label = 'test2';
figure(1, 'Position',[1,500,1,500]);
pcolor(xr/1e3,yr/1e3, EXT);
shading interp
colorbar
title([label ': East displacement [m]']);
xlabel('East distance [km]')
ylabel('North distance [km]')
print([label '_out_dE.png'])

figure(2, 'Position',[1,500,1,500]);
pcolor(xr/1e3,yr/1e3, EYT);
shading interp
colorbar
title([label ': North displacement [m]']);
xlabel('East distance [km]')
ylabel('North distance [km]')
print([label '_out_dN.png'])

figure(3, 'Position',[1,500,1,500]);
pcolor(xr/1e3,yr/1e3, EZT);
shading interp
colorbar
title([label ': Up displacement [m]']);
xlabel('East distance [km]')
ylabel('North distance [km]')
print([label '_out_dU.png'])

%%%
R = [xs/km, ys/km, zs/km, strike, dip, rake, slip, L/km, W/km];
save -ascii R R
% awk '{printf "%7.2f %7.2f %7.2f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n", $1, $2, $3, $4, $5, $6, $7, $8, $9}' R | head
