% output of data

fid=fopen('c:\matlab\tony\calcs\dosnewer2.f','r') %open file for read
%fprintf(fid,'concentration   %6.4f,  channel %4.1f  nm\n',conc,10^9*a)
%fprintf(fid,'scale energy  %5.2f mK\n\n',A)
%fprintf(fid,'temp    chempot    fraction  number   shiftno    chiconf   chibulk   ratio\n\n')

%fprintf(fid,'%6.2f %9.5f %9.5f %9.5f %9.5f %11.3e %11.3e %9.5f\n',OUT)
%D=fscanf(fid,'%6.2f %9.5f',inf);
[D,COUNT]=fscanf(fid,'%f',[2,inf]);
COUNT
fclose(fid)
DDD=D(2,:);
plot (DDD)
save
