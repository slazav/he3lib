function compile_libhe3()

  func_i0 = {'he3_pa' 'he3_pabn' 'he3_tabn' 'he3_amass'};

  func_i1 = {'he3_pmelt' 'he3_tc' 'he3_tab' 'he3_vm' 'he3_meff' ...
             'he3_pf' 'he3_vf' 'he3_gammaf' 'he3_dnde' ...
             'he3_tau_r' 'he3_tau_f'...
             'he3_yosida' 'he3_z0' 'he3_f0a'};

  func_i2 = {'he3_flegg' 'he3_swvel' 'he3_swvel_par' 'he3_swvel_per' ...
             'he3_ds_exp' 'he3_dn_exp' 'he3_d_exp' 'he3_susept'};

  function comp(name, narg)
    mex(['-DFUNC=' name '_ -DNARGIN=' num2str(narg)],...
         '-output', name, 'mexfunc.c', './libhe3.so');
  end

  for i = 1:length(func_i0) comp(func_i0{i}, 0); end
  for i = 1:length(func_i1) comp(func_i1{i}, 1); end
  for i = 1:length(func_i2) comp(func_i2{i}, 2); end
end
%ln -s $< $@


