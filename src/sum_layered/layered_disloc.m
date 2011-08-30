function [ux, uy, uz] = layered_disloc(xs, ys, zs, strike, dip, ...
				       rake, slip, L, W, npw, npy, ...
				       xr, yr, edks)
%[ux, uy, uz] = layered_disloc(xs, ys, zs, strike, dip, rake,
%                              slip, L, W, npw, npy, xr, yr, edks)
%
% --- INPUT ---
% --- SOURCE INFO
% --- 1D arrays, length = number of fault patches
% xs     m, east coord to middle of top edge of fault 
% ys     m, north coord to middle of top edge of fault 
% zs m, depth coord to middle of top edge of fault (+ down) 
% strike deg, clockwise from north 
% dip    deg, 90 is vertical 
% rake   deg, 0 strike slip, 90 dip slip 
% slip   m, slip in the rake direction
% L      m, along strike fault length 
% W      m, down dip fault length 
% npw    scalar, # point sources per patch in the down dip direction
% npy    scalar, # point sources per patch in the along strike direction
% --- RECEIVER INFO
% 1D arrays, length = number of receivers
% xr     m, east coordinate of receivers 
% yr     m, north coordinate of receivers 
% --- ELASTIC STRUCTURE INFO
% edks   string, full name of edks file, e.g., halfspace.edks
% --- OUTPUT ---
% --- 2D arrays, (receivers, fault patches)
% ux     m, east displacement
% uy     m, west displacement
% uz     m, up displacement (+ up)

np   = length(xs); % number of patches
nrec = length(xr); % number of receivers

%%%%% filenames

prefix = 'snoopy'; %prefix to be used for IO matlab/fortran
			   
file_rec = [prefix '.rec'];
file_pat = [prefix '.pat'];
file_dux = [prefix '_ux.dis'];
file_duy = [prefix '_uy.dis'];
file_duz = [prefix '_uz.dis'];
unix(['rm -f ' file_rec ' ' file_pat ' ' file_dux ' ' file_duy ' ' file_duz]);

%%%%% write receiver location file (observation points)

temp = [xr(:), yr(:)]';
fid  = fopen(file_rec,'w');
fwrite(fid, temp, 'real*4');
fclose(fid);

%%%%% write fault patch information

temp = [xs(:), ys(:), zs(:), strike(:), dip(:), rake(:), ...
	W(:), L(:), slip(:)]';
fid  = fopen(file_pat,'w');
fwrite(fid, temp, 'real*4');
fclose(fid);

%%%%% call sum_layered to calculate Greens functions

comstr = [' sum_layered ' edks ' ' prefix ' 'int2str(nrec) ' ' ...
	  int2str(np) ' ' int2str(npw) ' ' int2str(npy)];
unix(comstr);

%%%%% read sum_layered output Greens function

fid  = fopen(file_dux,'r');
[ux, count] = fread(fid, [nrec, np], 'real*4');
fclose(fid);

fid  = fopen(file_duy,'r');
[uy, count] = fread(fid, [nrec, np], 'real*4');
fclose(fid);

fid  = fopen(file_duz,'r');
[uz, count] = fread(fid, [nrec, np], 'real*4');
fclose(fid);
