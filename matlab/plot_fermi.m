function plot_pars()
  find_figure('He3 fermi-liquid pars'); clf;
  press=0:0.1:35;

  % fit Graywall 1983 gata
  pdat = [0   3   6   9   12  15  18  21  24  27  30  33  34.36];
  f0s = [915 1583 2222 2861 3497 4133 4803 5437 6102 6822 7560 8344 8709] / 100;
  f1s = [527 640 732 815 895 971 1047 1114 1180 1250 1320 1396 1428] / 100;
  f0a = [700 725 736 745 750 755 759 759 760 759 758 759 757] * -1e-3;
  f1a = [55 73 79 86 90 95 99 99 100 99 98 101 99] * -1e-2;
  vm  = [3684 3387 3207 3076 2971 2886 2813 2756 2706 2658 2614 2571 2554] /100;
  mm  = [276 313 344 372 398 424 449 471 493 517 540 565 576]/100;
  gg  = [274 295 312 328 343 358 373 387 400 413 427 442 449]/100;
  c1  = [182.9 227.5 259.7 285.9 308.0 327.1 345.0 360.5 375.1 389.3 403.0 415.9 421.7];
  tm  = [359 305 277 256 238 224 212 205 198 191 185 179 177]/1e3;


  hold on; title('Graywall-83 values');
  plot(pdat, f1s, 'r*');
  plot(pdat, f0s, 'g*');
  plot(pdat, f0a, 'b*');
  plot(pdat, f1a, 'm*');
  plot(pdat, vm, 'ro');
  plot(pdat, mm, 'bo');
  plot(pdat, gg, 'go');
  plot(pdat, c1, 'mo');

  legend('f1s', 'f0s', 'f0a', 'f1a', 'Vm', 'm_eff/m', 'gamma');

  p1s=polyfit(pdat, f1s, 6);
  fprintf('f1s: %e\n', p1s);
  plot(press, polyval(p1s, press), 'r-');

  p0s=polyfit(pdat, f0s, 6);
  fprintf('f0s: %e\n', p0s);
  plot(press, polyval(p0s, press), 'g-');

  p0a=polyfit(pdat, f0a, 6);
  fprintf('f0a: %e\n', p0a);
  plot(press, polyval(p0a, press), 'b-');

  p1a=polyfit(pdat, f1a, 6);
  fprintf('f1a: %e\n', p1a);
  plot(press, polyval(p1a, press), 'm-');

  pvm=polyfit(pdat, vm, 6);
  fprintf('vm: %e\n', pvm);
  plot(press, polyval(pvm, press), 'r-');

  pmm=polyfit(pdat, mm, 6);
  fprintf('mm: %e\n', pmm);
  plot(press, polyval(pmm, press), 'b-');

  pgg=polyfit(pdat, gg, 6);
  fprintf('gg: %e\n', pgg);
  plot(press, polyval(pgg, press), 'g-');

  pc1=polyfit(pdat, c1, 6);
  fprintf('c1: %e\n', pc1);
  plot(press, polyval(pc1, press), 'm-');

  ptm=polyfit(pdat, tm, 6);
  fprintf('tm: %e\n', ptm);
  plot(press, polyval(ptm, press), 'm-');

end


