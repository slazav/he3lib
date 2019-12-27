% output of data

OUT=[T; MU; frac(1,:); N(1,:); N(2,:); CHI; BU; RAT;]

fid=fopen('c:\matlab\tony\calcs\qwd001aug10mK.f','w') %open file for write

fprintf(fid,'concentration   %6.4f,  channel %4.1f  nm\n',conc,10^9*a)
fprintf(fid,'scale energy  %5.2f mK\n\n',A)
fprintf(fid,'temp    chempot    fraction  number   shiftno    chiconf   chibulk   ratio\n\n')

fprintf(fid,'%6.2f %9.5f %9.5f %9.5f %9.5f %11.3e %11.3e %9.5f\n',OUT)

fclose(fid)