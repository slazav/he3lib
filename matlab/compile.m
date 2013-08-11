function compile()

  func_i0 = {...
    'he3_gyro' 'he3_amass' 'he3_mmass'...
    'const_na' 'const_kb' 'const_r' 'const_h' 'const_hbar'...
    'const_pi'...
    'he3_pcr' 'he3_tcr' 'he3_pabn' 'he3_tabn'...
    'he3_psmin' 'he3_tsmin' 'he3_pa' 'he3_ta'...
    'he3_pb' 'he3_tb' 'he3_pneel' 'he3_tneel'...
  };

  func_i1 = {...
    'he3_pvap' 'he3_pmelt' 'he3_pmelt_gr' 'he3_tc' 'he3_tab' ...
    'he3_vm' 'he3_gammaf' 'he3_c1' 'he3_tmag' ...
    'he3_rho' 'he3_2n0' 'he3_mm' 'he3_meff' ...
    'he3_pf' 'he3_vf' 'he3_chi0'...
    'he3_f0a' 'he3_f0s' 'he3_f1a' 'he3_f1s' 'he3_z0' ...
    'he3_a' 'he3_gdk' 'he3_tfeff' 'he3_lscatt' ...
    'he3_bcsgap' ...
    'he3_tau_r' 'he3_tau_f'...
  };

  func_i2 = {...
    'he3_trivgap' 'he3_z3' 'he3_z5' 'he3_z7' 'he3_yosida0' ...
    'he3_flegg' 'he3_swvel' 'he3_swvel_par' 'he3_swvel_per' ...
    'he3_ds_exp' 'he3_dn_exp' 'he3_d_exp' 'he3_susept'...
  };

  function comp_octave(name, narg)
    mex(['-DFUNC=' name '_ -DNARGIN=' num2str(narg)],...
         '--output', name, '-lhe3',...
         ['-L' pwd ], ['-Wl,--rpath=' pwd ],...
         'mexfunc.c');
  end
  function comp_matlab(name, narg)
    mex(['-DFUNC=' name '_ -DNARGIN=' num2str(narg)],...
         '-output', name, 'mexfunc.c', './libhe3.so');
  end

  if length(ver('Octave')) comp=@comp_octave;
  else comp=@comp_matlab; end

  for i = 1:length(func_i0) comp(func_i0{i}, 0); end
  for i = 1:length(func_i1) comp(func_i1{i}, 1); end
  for i = 1:length(func_i2) comp(func_i2{i}, 2); end
end
%ln -s $< $@


