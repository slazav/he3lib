% output of data

OUT=[t/3600; tn; te; th; the; b; koroverkap; boverth;  ss;];
%OUT=[E; DOS3]
fid=fopen('D:\AMG\matlab_files\TONY\demag\run15eddy02.dat','w'); %open file for write

%fprintf(fid,'concentration   %6.4f,  channel %4.1f  nm\n',conc,10^9*a)
%fprintf(fid,'energy        density of states')
%fprintf(fid,'heat leaks  %6.1f, %6.1f \n\n',qh,qhd)
fprintf(fid,' hours    Tn     Te     Th1    Th2    B    ratio  BoverT  entropy\n\n');

fprintf(fid,' %9.4f %9.5e %9.5e %9.5e %9.5e %9.5e %9.5e  %9.5e   %9.5e\n',OUT);
%fprintf(fid,'%6.2f %9.5f \n',OUT)

fclose(fid)

