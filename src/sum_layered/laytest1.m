clear all;

%%%%% source info

npw = 1; % points per patch in the down dip direction
npy = 1; % points per patch in the along strike direction

strike = [0.0];     %deg
dip    = [60.0];   %deg
rake   = [90.0];   %deg
slip   = [5e0]; %meters
L      = [1e1]; %meters, along strike
W      = [1e1]; %meters, down dip

xs     = [1.0e3]; %middle of top edge of fault
ys     = [1.0e3]; %middle of top edge of fault
zs     = [1.0e3]; %middle of top edge of fault

np  = length(zs);  % number of patches

%%%%% receiver info

nx    = 31;
ny    = 31;
xrmin = 0;
yrmin = 0;
xrmax = 5e3;
yrmax = 5e3;
dxr   = (xrmax - xrmin)/(nx-1);
dyr   = (yrmax - yrmin)/(ny-1);
xr    = [xrmin:dxr:xrmax];
yr    = [yrmin:dyr:yrmax];
[xrg, yrg] = meshgrid(xr,yr);

%%%%% filenames

edks   = 'halfspace.edks'; % edks file, must have hdr.*.edks in directory
			   
%%%%% call gateway routine to fortran 90 layered code, all input
%%%%% arrays must be 1D, output arrays are 2D (nreceivers, npatches)

[uxt, uyt, uzt] = layered_disloc(xs, ys, zs, strike, dip, rake, slip, ...
				 L, W, npw, npy, xrg(:), yrg(:), edks);

%%%%% postprocessing and benchmark against okada

for ip = 1:np
  
  %%%%% Okada test

  ftype  = 2;    % 1=strike,2=dip3=tensile,4=mogi
  nu     = 0.25; % Poisson's ratio 

  dipr    = dip * pi / 180;
  striker = strike * pi / 180;
  zs_bot  = zs + W.*sin(dipr); 

  xshift = -xs(ip)-W(ip)*cos(dipr(ip))*cos(striker(ip));
  yshift = -ys(ip)+W(ip)*cos(dipr(ip))*sin(striker(ip));
  
  [uxo, uyo, uzo] = calc_okada(slip(ip), xrg(:)+xshift, yrg(:)+yshift, ...
			       nu, dip(ip), zs_bot(ip), L(ip), W(ip), ...
			       ftype, strike(ip));

  %%%%% reformat, assumes one patch

  ux = reshape(uxt(:,ip),nx,ny);
  uy = reshape(uyt(:,ip),nx,ny);
  uz = reshape(uzt(:,ip),nx,ny);
  uxo = reshape(uxo,nx,ny);
  uyo = reshape(uyo,nx,ny);
  uzo = reshape(uzo,nx,ny);
  armax  = num2str(max(abs( ux(:) +  i*uy(:))));
  armaxo = num2str(max(abs(uxo(:) + i*uyo(:))));

  %%%%% plotting
  
  figure(ip)
  
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
  xlabel('x, km'), ylabel('y, km')
  title(['Okada: max Uhoriz ' armaxo]), colorbar('v')
  axis tight
  
  dux = ux-uxo;
  duy = uy-uyo;
  duz = uz-uzo;
  
  subplot(3,2,3)
  pcolor(xr/1e3, yr/1e3, duz), hold on
  shading interp
  quiver(xr/1e3,yr/1e3,dux,duy), hold off
  colormap(jet), 
  xlabel('x, km'), ylabel('y, km')
  title('Difference, Uz'), colorbar('v')
  axis tight
  
  subplot(3,2,4)
  pcolor(xr/1e3, yr/1e3, dux)
  shading interp
  colormap(jet), 
  xlabel('x, km'), ylabel('y, km')
  title('Difference, Ux'), colorbar('v')
  axis tight
  
  subplot(3,2,5)
  pcolor(xr/1e3, yr/1e3, duy)
  shading interp
  colormap(jet), 
  xlabel('x, km'), ylabel('y, km')
  title('Difference, Uy'), colorbar('v')
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
  
  orient tall
  wysiwyg
  
end


