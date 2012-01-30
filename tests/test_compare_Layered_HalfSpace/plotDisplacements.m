close all
clear all
colorLim = [-10, -1];
He = load('He.dat');
Hn = load('Hn.dat');
Hu = load('Hu.dat');

Te = load('Te.dat');
Tn = load('Tn.dat');
Tu = load('Tu.dat');

r = load('Dist.dat');
d = load('Dep.dat');


% Plot without normalization

%cmap = jet(256);


load cmap.mat
figure(1)
subplot(3,2,1);
imagesc(r,d,log10(abs((Te)))); colorbar();
set(gca, 'Xscale','log')
colormap(cmap)
caxis(colorLim)
ylabel('Depth')
title('log10(LayeredSpace pred) - Horizontal Disp')

subplot(3,2,3);
imagesc(r,d,log10(abs((He)))); colorbar();
set(gca, 'Xscale','log')
colormap(cmap)
caxis(colorLim)
ylabel('Depth')
title('log10(HomogeneousSpace pred) - Horizontal Disp')

subplot(3,2,5);
imagesc(r,d,log10(abs((Te - He)))); colorbar();
set(gca, 'Xscale','log')
colormap(cmap)
caxis(colorLim)
xlabel('Source - Receiver Horizontal Distance')
ylabel('Depth')
title('log10(LayeredSpace pred - HomogeneousSpace pred) - Horizontal Disp')


subplot(3,2,2);
imagesc(r,d,log10(abs((Tu)))); colorbar();
set(gca, 'Xscale','log')
colormap(cmap)
caxis(colorLim)
ylabel('Depth')
title('log10(LayeredSpace pred) - Vertical Disp')

subplot(3,2,4);
imagesc(r,d,log10(abs((Hu)))); colorbar();
set(gca, 'Xscale','log')
colormap(cmap)
caxis(colorLim)
ylabel('Depth')
title('log10(HomogeneousSpace pred) - Vertical Disp')



subplot(3,2,6);
imagesc(r,d,log10(abs((Tu - Hu)))); colorbar();
set(gca, 'Xscale','log')
colormap(cmap)
%caxis(colorLim)
xlabel('Source - Receiver Horizontal Distance')
ylabel('Depth')
title('log10(LayeredSpace pred - HomogeneousSpace pred) - Vertical Disp')


figure(2)
load 2signCmap.mat

subplot(1,2,1);
imagesc(r,d,Te - He); colorbar();
set(gca, 'Xscale','log')
colormap(cmap)
caxis([-1.0E-3, 1.0E-3])
xlabel('Source - Receiver Horizontal Distance')
ylabel('Depth')
title('Horizontal: LayeredSpace pred - HomogeneousSpace pred')

subplot(1,2,2);
imagesc(r,d,Tu - Hu); colorbar();
set(gca, 'Xscale','log')
colormap(cmap)
caxis([-1.0E-3, 1.0E-3])
xlabel('Source - Receiver Horizontal Distance')
ylabel('Depth')
title('Vertical: LayeredSpace pred - HomogeneousSpace pred')



% Plot normalizing the different depths by the closest receiver of the
% Homogeneous half space values.
%cmap = jet(256);

% get the normalization values
[Ndep, Nrec] = size(He);
for i = 1:Ndep
    Te(i,:) = Te(i,:)/max(abs(He(i,:)));
    Tn(i,:) = Tn(i,:)/max(abs(Hn(i,:)));
    Tu(i,:) = Tu(i,:)/max(abs(Hu(i,:)));  
    He(i,:) = He(i,:)/max(abs(He(i,:)));
    Hn(i,:) = Hn(i,:)/max(abs(Hn(i,:)));
    Hu(i,:) = Hu(i,:)/max(abs(Hu(i,:)));  
end


load 2signCmap.mat
colorLim = [-1,1];
figure(3)
subplot(3,2,1);
imagesc(r,d,(((Te)))); colorbar();
set(gca, 'Xscale','log')
colormap(cmap)
caxis(colorLim)
ylabel('Depth')
title('Horizontal: LayeredSpace pred (Normalized)')

subplot(3,2,3);
imagesc(r,d,(((He)))); colorbar();
set(gca, 'Xscale','log')
colormap(cmap)
caxis(colorLim)
ylabel('Depth')
title('Horizontal: HomogeneousSpace pred (Normalized)')

subplot(3,2,5);
imagesc(r,d,(((Te - He)))); colorbar();
set(gca, 'Xscale','log')
colormap(cmap)
caxis([-0.35, 0.35])
xlabel('Source - Receiver Horizontal Distance')
ylabel('Depth')
title('Horizontal: LayeredSpace pred - HomogeneousSpace pred (Normalized)')


subplot(3,2,2);
imagesc(r,d,(((Tu)))); colorbar();
set(gca, 'Xscale','log')
colormap(cmap)
caxis(colorLim)
ylabel('Depth')
title('Vertical: LayeredSpace pred (Normalized)')

subplot(3,2,4);
imagesc(r,d,(((Hu)))); colorbar();
set(gca, 'Xscale','log')
colormap(cmap)%
caxis(colorLim)
ylabel('Depth')
title('Vertical: HomogeneousSpace pred (Normalized)')



subplot(3,2,6);
imagesc(r,d,(((Tu - Hu)))); colorbar();
set(gca, 'Xscale','log')
colormap(cmap)
caxis([-0.35, 0.35])
xlabel('Source - Receiver Horizontal Distance')
ylabel('Depth')
title('Vertical: LayeredSpace pred - HomogeneousSpace pred (Normalized)')


