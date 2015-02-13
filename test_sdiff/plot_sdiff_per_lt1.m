function plot_sdiff_per_lt1()

  addpath ../matlab

  figure; clf; hold on;

  f=833000;

  ttc=[0.14 0.15 0.16];
  col='rgbcmk';

  for i=1:length(ttc)
    p=0:0.05:30;
    Dxx=he3_diff_perp_xx(ttc(i), p, f);
    Dzz=he3_diff_perp_zz(ttc(i), p, f);
    plot(p, Dxx, [col(i) '-']);

    Dxx=he3_diff_perp_xx_im(ttc(i), p, f);
    Dzz=he3_diff_perp_zz_im(ttc(i), p, f);
    plot(p, Dxx, [col(i) '--']);
  end

end