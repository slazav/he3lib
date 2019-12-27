% output of data
%E=[1:1:length(DDD)];
% OUT=[T; MU; frac(1,:); N(1,:); N(2,:); CHI; BU; RAT;]
%OUT=[E; DDD];
fid=fopen('c:\matlab\tony\calcs\top.f','w') %open file for write

%fprintf(fid,'concentration   %6.4f,  channel %4.1f  nm\n',conc,10^9*a)
%fprintf(fid,'energy        density of states')
%fprintf(fid,'scale energy  %5.2f mK\n\n',A)
%fprintf(fid,'temp    chempot    fraction  mumber   shiftno    chiconf   chibulk   ratio\n\n')

%fprintf(fid,'%6.2f %9.5f %9.5f %9.5f %9.5f %11.3e %11.3e %9.5f\n',OUT)
fprintf(fid,'%9.5f \n',top)

fclose(fid)