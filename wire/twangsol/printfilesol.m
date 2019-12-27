% output of data

OUT=[T; F1; F2;]%; frac(1,:); N(1,:); N(2,:); CHI; BU; RAT;]

fid=fopen('c:\matlab\tony\calcs\bananah4a.dat','w') %open file for write

%fprintf(fid,'concentration   %6.4f,  channel %4.1f  nm\n',conc,10^9*a)
%fprintf(fid,'scale energy  %5.2f mK\n\n',A)
%fprintf(fid,'temp    chempot    fraction  number   shiftno    chiconf   chibulk   ratio\n\n')

fprintf(fid,'%9.4e %10.4e %10.4e\n',OUT)

fclose(fid)
