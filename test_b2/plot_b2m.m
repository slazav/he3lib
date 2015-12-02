function plot_b2gap()
% b2 phase magnetization and susceptibility

  figure(1); clf; hold on;

  p=0;

  for ttc=0.1:0.1:0.9;
    gap=he3_gap(ttc,p);

    H  =0:100:4000;
    Hc = he3_b2hcr(ttc,p);

    He = he3_b2heff(ttc,p,H);
    M = (H-He)*he3_2n0(p)*(he3_gyro*const_hbar/2)^2/he3_f0a(p);
    chi  = (M(2:end)-M(1:end-1))./(H(2:end)-H(1:end-1));
    chi0 = he3_chi_b(ttc,p)*he3_chi_n(p);

    plot(H(1:end-1), chi/chi0, 'b-')
    if ttc==0.9;     t=sprintf('T/Tc = %.1f',ttc);
    else             t=sprintf('%.1f',ttc);
    end
    if ttc==0.2; t=[' ' t]; end
    if ttc==0.1; t=['   ' t]; end
    text(H(find(chi/chi0>3,1))-150, 3.02, t);
  end
  plot([H(1) H(end)], [1 1], 'k-')
  ylim([0 3])
  xlabel('H, G')
  ylabel('\chi(H)/\chi(H=0)')
end
