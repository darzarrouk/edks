close all
clear all

if 1 == 0
% load the Okada GFs
OeSS = load('./OKADA_GFs/O_GeSS.txt');
OnSS = load('./OKADA_GFs/O_GnSS.txt');
OuSS = load('./OKADA_GFs/O_GuSS.txt');
OeDS = load('./OKADA_GFs/O_GeDS.txt');
OnDS = load('./OKADA_GFs/O_GnDS.txt');
OuDS = load('./OKADA_GFs/O_GuDS.txt');

% load the EDKS GFs
EeSS = load('./EDKS_GF/GeSS.txt');
EnSS = load('./EDKS_GF/GnSS.txt');
EuSS = load('./EDKS_GF/GuSS.txt');
EeDS = load('./EDKS_GF/GeDS.txt');
EnDS = load('./EDKS_GF/GnDS.txt');
EuDS = load('./EDKS_GF/GuDS.txt');

save GFmatrices.mat
else 
    load GFmatrices.mat 
end

colorLim = [-8, -1];
load colormap.mat


% plot the strike slip Okada and EDKS GF
figure(1)
subplot(3,3,1); pcolor(log10(abs(OeSS))); shading flat; colorbar(); caxis(colorLim); title('Okada E StrikeSlip')
subplot(3,3,2); pcolor(log10(abs(OnSS))); shading flat; colorbar(); caxis(colorLim); title('Okada N StrikeSlip')
subplot(3,3,3); pcolor(log10(abs(OuSS))); shading flat; colorbar(); caxis(colorLim); title('Okada U StrikeSlip')

subplot(3,3,4); pcolor(log10(abs(EeSS))); shading flat; colorbar(); caxis(colorLim); title('EDKS E StrikeSlip')
subplot(3,3,5); pcolor(log10(abs(EnSS))); shading flat; colorbar(); caxis(colorLim); title('EDKS N StrikeSlip')
subplot(3,3,6); pcolor(log10(abs(EuSS))); shading flat; colorbar(); caxis(colorLim); title('EDKS U StrikeSlip')

subplot(3,3,7); pcolor(log10(abs((OeSS - EeSS)))); shading flat; colorbar(); caxis(colorLim); title('Diff E StrikeSlip')
subplot(3,3,8); pcolor(log10(abs((OnSS - EnSS)))); shading flat; colorbar(); caxis(colorLim); title('Diff N StrikeSlip')
subplot(3,3,9); pcolor(log10(abs((OuSS - EuSS)))); shading flat; colorbar(); caxis(colorLim); title('Diff U StrikeSlip')
colormap(cmap)

figure(2)
subplot(3,3,1); pcolor(log10(abs(OeDS))); shading flat; colorbar(); caxis(colorLim); title('Okada E DipSlip')
subplot(3,3,2); pcolor(log10(abs(OnDS))); shading flat; colorbar(); caxis(colorLim); title('Okada N DipSlip')
subplot(3,3,3); pcolor(log10(abs(OuDS))); shading flat; colorbar(); caxis(colorLim); title('Okada U DipSlip')

subplot(3,3,4); pcolor(log10(abs(EeDS))); shading flat; colorbar(); caxis(colorLim); title('EDKS E DipSlip')
subplot(3,3,5); pcolor(log10(abs(EnDS))); shading flat; colorbar(); caxis(colorLim); title('EDKS N DipSlip')
subplot(3,3,6); pcolor(log10(abs(EuDS))); shading flat; colorbar(); caxis(colorLim); title('EDKS U DipSlip')

subplot(3,3,7); pcolor(log10(abs((OeDS - EeDS)))); shading flat; colorbar(); caxis(colorLim); title('Diff E DipSlip')
subplot(3,3,8); pcolor(log10(abs((OnDS - EnDS)))); shading flat; colorbar(); caxis(colorLim); title('Diff N DipSlip')
subplot(3,3,9); pcolor(log10(abs((OuDS - EuDS)))); shading flat; colorbar(); caxis(colorLim); title('Diff U DipSlip')
colormap(cmap)