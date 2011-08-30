clear all;

%%%%% source info

npw = 10; % points per patch in the down dip direction
npy = 10; % points per patch in the along strike direction
np  = 1; % number of patches

strike = 0.0;   %deg
dip    = 45.0;  %deg
rake   = 0.0;   %deg
slip   = 2.0e2; %meters
L      = 0.5e3; %meters
W      = 0.5e3; %meters

xs     = 1.0e3; %middle of top edge of fault
ys     = 1.0e3; %middle of top edge of fault
zs_top = 1.0e3; %middle of top edge of fault

dipr    = dip * pi / 180;
striker = strike * pi / 180;

zs_mid = zs_top + 0.5*W*sin(dipr); %fault depth, meters (middle of fault patch)
zs_bot = zs_top + W*sin(dipr);     %fault depth, meters (bottom of fault patch)

%%%%% receiver info

nx    = 31;
ny    = 31;
xrmin = 0;
yrmin = 0;
xrmax = 3e3;
yrmax = 3e3;
dxr   = (xrmax - xrmin)/(nx-1);
dyr   = (yrmax - yrmin)/(ny-1);

xr = [xrmin:dxr:xrmax];
yr = [yrmin:dyr:yrmax];
[xrg, yrg]=meshgrid(xr,yr);
nrec = length(xrg(:));

%%%%% Okada test

ftype  = 1;    % 1=strike,2=dip3=tensile,4=mogi
nu     = 0.25; % Poisson's ratio 
xshift = -xs-W*cos(dipr)*cos(striker);
yshift = -ys+W*cos(dipr)*sin(striker);
[uxo, uyo, uzo] = calc_okada(slip, xrg(:)+xshift, yrg(:)+yshift, nu, dip, ...
			     zs_bot, L, W, ftype, strike);
%%%%% filenames

edks   = 'halfspace.edks'; % edks file, must have hdr.*.edks in directory
prefix = 'inverse'; %prefix to be used for IO matlab/fortran
			   
%%%%% write receiver location file (observation points)

temp = [xrg(:), yrg(:)]';
fid  = fopen([prefix '.rec'],'w');
fwrite(fid, temp, 'real*4');
fclose(fid);

%%%%% write fault patch information

temp = [xs(:), ys(:), zs_top(:), strike(:), dip(:), rake(:), W(:), L(:), slip(:)]';
fid  = fopen([prefix '.pat'],'w');
fwrite(fid, temp, 'real*4');
fclose(fid);

%%%%% call sum_layered to calculate Greens functions

time1  = cputime;
comstr = [' sum_layered.sgi ' edks ' ' prefix ' 'int2str(nrec) ' ' ...
	  int2str(np) ' ' int2str(npw) ' ' int2str(npy)]
unix(comstr);
dtime  = cputime-time1;
disp(['elapsed time in layered code: ' int2str(dtime) ' sec']);

%%%%% read sum_layered output Greens function

fid  = fopen([prefix '_ux.dis'],'r');
[ux, count] = fread(fid, [nrec, np], 'real*4');
fclose(fid);
fid  = fopen([prefix '_uy.dis'],'r');
[uy, count] = fread(fid, [nrec, np], 'real*4');
fclose(fid);
fid  = fopen([prefix '_uz.dis'],'r');
[uz, count] = fread(fid, [nrec, np], 'real*4');
fclose(fid);

%%%%% reformat, assumes one patch

ux = reshape(ux,nx,ny);
uy = reshape(uy,nx,ny);
uz = reshape(uz,nx,ny);
uxo = reshape(uxo,nx,ny);
uyo = reshape(uyo,nx,ny);
uzo = reshape(uzo,nx,ny);
armax  = num2str(max(abs( ux(:) +  i*uy(:))));
armaxo = num2str(max(abs(uxo(:) + i*uyo(:))));

%%%%% plotting

subplot(3,2,1)
pcolor(xr/1e3, yr/1e3, uz), hold on
shading interp
quiver(xr/1e3,yr/1e3,ux,uy), hold off
colormap(jet), 
xlabel('x, km'), ylabel('y, km'), title(['EDKS: max Uhoriz ' armax])
colorbar('v')
axis tight

subplot(3,2,2)
pcolor(xr/1e3, yr/1e3, uzo), hold on
shading interp
quiver(xr/1e3,yr/1e3,uxo,uyo), hold off
colormap(jet), 
xlabel('x, km'), ylabel('y, km'), title(['Okada: max Uhoriz ' armaxo]), colorbar('v')
axis tight

dux = ux-uxo;
duy = uy-uyo;
duz = uz-uzo;

subplot(3,2,3)
pcolor(xr/1e3, yr/1e3, duz), hold on
shading interp
quiver(xr/1e3,yr/1e3,dux,duy), hold off
colormap(jet), 
xlabel('x, km'), ylabel('y, km'), title('Difference, Uz'), colorbar('v')
axis tight

subplot(3,2,4)
pcolor(xr/1e3, yr/1e3, dux)
shading interp
colormap(jet), 
xlabel('x, km'), ylabel('y, km'), title('Difference, Ux'), colorbar('v')
axis tight

subplot(3,2,5)
pcolor(xr/1e3, yr/1e3, duy)
shading interp
colormap(jet), 
xlabel('x, km'), ylabel('y, km'), title('Difference, Uy'), colorbar('v')
axis tight

for ir = 1:nx
  uzp(ir) = uz(ir,ir);
  uzop(ir) = uzo(ir,ir);
  rd(ir) = abs(xr(ir)+i*yr(ir));
end 

rd = rd/1e3;

subplot(3,2,6)
ax1=plot(rd,uzp,'o',rd,uzop,'+'); hold on
ax2=plot(rd,uzp,rd,uzop);
ax3=plot(rd,(uzp-uzop)/max(abs(uzop(:)))*100); hold off
