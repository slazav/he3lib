%function for viscosity - raw de Waele
%t in K

function eta = visc(t)
eta = .277e-7/t^2 +3.4e-7/t;  %original Zeegers et al

if t< .0165,
   eta=.305e-7/t^2+1.35e-7/t+2.2e-6; 
   % improved version 
			%allowing for mfp effects
         % reasonable fit up to 100mK
end
      
      
         
if t> .068,
  eta=.29e-7/t^2+1.65e-7/t+2.3e-6; % improved version 
			%allowing for mfp effects
         % reasonable fit up to 100mK
end
      

