xr = [5.  10.  12. 40.];
yr = [40. -45. 45. 23.];


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

