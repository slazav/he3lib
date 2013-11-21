function compile(v)

  func_i0 = {...
    'he3_gyro' 'he3_amass' 'he3_mmass' ...
    'const_na' 'const_kb' 'const_r' 'const_h' 'const_hbar' ...
    'const_pi' ...
    'he3_pcr' 'he3_tcr' 'he3_pabn' 'he3_tabn' ...
    'he3_pm' 'he3_tm' 'he3_pa' 'he3_ta' ...
    'he3_pb' 'he3_tb' 'he3_ps' 'he3_ts' ...
    'he3_pa_plts' 'he3_ta_plts' ...
    'he3_pb_plts' 'he3_tb_plts' 'he3_ps_plts' 'he3_ts_plts' ...
  };

  func_i1 = {...
    'he3_pvap' 'he3_pmelt' 'he3_tc' 'he3_tab' ...
    'he3_pmelt_plts', 'he3_plts2gr', 'he3_gr2plts'...
    'he3_vm' 'he3_gammaf' 'he3_c1' 'he3_tmag' ...
    'he3_rho' 'he3_2n0' 'he3_mm' 'he3_meff' ...
    'he3_pf' 'he3_vf' 'he3_chi_n' ...
    'he3_f0a' 'he3_f0s' 'he3_f1a' 'he3_f1s' ...
    'he3_a' 'he3_gdk' 'he3_tfeff' ...
    'he3_bcsgap' 'he3_bcsgap_fast' ...
    'he3_crsect_w' 'he3_crsect_wi' 'he3_crsect_wd' 'he3_crsect_wl'...
    'he3_scatt_l1a' 'he3_scatt_g0' 'he3_scatt_d0' 'he3_scatt_w0'...
    'he3_tau_r' 'he3_tau_f' ...
  };

  func_i2 = {...
    'he3_cv_n' ...
    'he3_trivgap' 'he3_yosida_par' 'he3_yosida_perp' ...
    'he3_z3' 'he3_z5' 'he3_z7' 'he3_lambda' ...
    'he3_rho_nb' 'he3_chi_b'...
    'he3_tau_n0' 'he3_tau_n_av' 'he3_tau_nd' ...
    'he3_tau0' 'he3_tau0lt' 'he3_tau_av' 'he3_fpath' ...
    'he3_tau_dperp' 'he3_tau_dpar'...
    'he3_diffn_hydr' 'he3_diff_hperp_zz' 'he3_diff_hpar_zz'...
    'he3_swvel' 'he3_swvel_par' 'he3_swvel_per' ...
    'he3_ds_exp' 'he3_dn_exp' 'he3_d_exp' ...
    'rota_c_ns' 'he3_nu_b' ...
  };

  func_i3 = { ...
    'he3_yosida' ...
    'he3_diffn_perp' ...
    'he3_diff_perp_xx' 'he3_diff_perp_zz' ...
    'he3_diff_perp_xx_im' 'he3_diff_perp_zz_im' ...
    'he3_diff_par_xx' 'he3_diff_par_zz' ...
    'rota_fork_cal' ...
  };

  func_i4 = { ...
  };

  func_i5 = { ...
    'he3_coll_int' 'he3_coll_int_lt' 'he3_coll_int_ht' ...
  };

  % note: old octave does not support ver('Octave') call
  if (nargin < 1) v=length(ver('Octave')); end
  if v comp=@comp_octave;
  else comp=@comp_matlab; end

  for i = 1:length(func_i0) comp(func_i0{i}, 0); end
  for i = 1:length(func_i1) comp(func_i1{i}, 1); end
  for i = 1:length(func_i2) comp(func_i2{i}, 2); end
  for i = 1:length(func_i3) comp(func_i3{i}, 3); end
  for i = 1:length(func_i4) comp(func_i4{i}, 4); end
  for i = 1:length(func_i5) comp(func_i5{i}, 5); end
end
%ln -s $< $@

function comp_octave(name, narg)
  mex(['-DFUNC=' name '_ -DNARGIN=' num2str(narg)],...
       '--output', name, '-lhe3',...
       ['-L' pwd ], ['-Wl,--rpath=' pwd ],...
       'mexfunc.c');
end
function comp_matlab(name, narg)
  mex(['-DFUNC=' name '_ -DNARGIN=' num2str(narg)],...
       '-output', name, '-lhe3',...
       ['-L' pwd ], ['-Wl,--rpath=' pwd ],...
       'mexfunc.c');
end

