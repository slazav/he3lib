% output of data

OUT=[t; chi; invte; invtn; q;];
%OUT=[E; DOS3]
fid=fopen('c:\matlab\tony\calcs\cusimdm3.dat','w') %open file for write

%fprintf(fid,'concentration   %6.4f,  channel %4.1f  nm\n',conc,10^9*a)
%fprintf(fid,'energy        density of states')
fprintf(fid,'alpha  %6.5f \n\n',alpha)
fprintf(fid,'time   chi    electron    nuclei  heatleak\n\n')

fprintf(fid,'%6.2f  %9.1f %9.1f %9.1f %11.3e \n',OUT)
%fprintf(fid,'%6.2f %9.5f \n',OUT)

fclose(fid)