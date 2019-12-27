function etar=redvis(tr)

%  a plot of CHH  reduced viscosity gives the following reasonable
%  form.   Basically, LOG( redvis) vz sqrt(1-T)^.5 reasonably slow
  
  
     if tr<0.6
       etar=0.11;
     elseif tr<0.7  
       SLOPE=-0.8562;
       D0=0.4;
       N0=-0.9586;
       etar=10^(N0+SLOPE*(sqrt(1-tr)-sqrt(D0)));
     elseif tr<0.8 
       SLOPE=-0.6183;
       D0=0.3;
       N0=-0.8861;
       etar=10^(N0+SLOPE*(sqrt(1-tr)-sqrt(D0)));
     elseif tr<0.9 
       SLOPE=-1.4172;
       D0=0.2;
       N0=-0.8239;
       etar=10^(N0+SLOPE*(sqrt(1-tr)-sqrt(D0)));
     elseif tr<0.95
       SLOPE=-1.7352;
       D0=0.1;
       N0=-0.6383;
       etar=10^(N0+SLOPE*(sqrt(1-tr)-sqrt(D0)));
     elseif tr<0.975
       SLOPE=-1.6177;
       D0=0.05;
       N0=-0.4776;
       etar=10^(N0+SLOPE*(sqrt(1-tr)-sqrt(D0)));
     elseif tr<1
       SLOPE=-2.3503;
       D0=0.025;
       N0=-0.3716;
       etar=10^(N0+SLOPE*(sqrt(1-tr)-sqrt(D0)));
     else 
         etar=1;
              
     end
     


     

