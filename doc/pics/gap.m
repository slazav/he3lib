#!/usr/bin/octave -qf

figure;
addpath ../../matlab
subplot(3,2,1); hold on;

  ttc=[0.001:0.001:1.1];


  gg00=he3_trivgap(ttc, 0);
  gg30=he3_trivgap(ttc, 30);

  plot(ttc, gg00, 'k-', 'linewidth', 2)
  plot(ttc, gg30, 'k-')
  plot(ttc, he3_bcsgap(ttc), 'r-')

  plot(ttc, he3_yosida(ttc,gg00,0), 'r', 'linewidth', 2)
  plot(ttc, he3_yosida(ttc,gg00,2), 'g', 'linewidth', 2)
  plot(ttc, he3_yosida(ttc,gg00,4), 'b', 'linewidth', 2)

  plot(ttc, he3_yosida(ttc,gg30,0), 'r')
  plot(ttc, he3_yosida(ttc,gg30,2), 'g')
  plot(ttc, he3_yosida(ttc,gg30,4), 'b')

  legend('P=0', 'P=30', 'location', 'southwest')

%  plot(ttc, he3_z3(ttc,gg), 'r')
%  plot(ttc, he3_z5(ttc,gg), 'g')
%  plot(ttc, he3_z7(ttc,gg), 'b')
  text(0.5, 1.5, '\Delta_{\rm BCS}', 'fontsize', 8, 'color', 'r');
  text(0.8, 1.5, '\Delta', 'fontsize', 8, 'color', 'k');
  text(0.73, 0.7,  'Y_0', 'fontsize', 8, 'color', 'r');
  text(0.83, 0.5, 'Y_2', 'fontsize', 8, 'color', 'g');
  text(0.93, 0.2, 'Y_4', 'fontsize', 8, 'color', 'b');
  xlim([0 1.02]);

subplot(3,2,2); hold on;
  p=0;

  plot(ttc, he3_chi_b(ttc,00), 'm', 'linewidth', 2)
  plot(ttc, he3_chi_b(ttc,30), 'm')
  xlim([0 1.02]);
  ylim([0 3e-8]);
  legend('P=0', 'P=30')
  text(0.2, 2e-8, '\chi_B', 'fontsize', 10);

subplot(3,2,4:6); hold on;

  % data from Hakonen, et al., LJTP 76, p.225-283 (1989)
  ttc005=[762 769 791 837 847 902 919 941]/1000;
  nub005=[614 588 509 386 365 214 169 118] * 1e7;

  ttc050=[627 637 656 679 697 715 734 743 765 ...
          810 835 848 864 881 893 905 916 927 940 950]/1000;
  nub050=[1959 1870 1766 1636 1544 1438 1351 1299 1176 ...
          896 741 661 579 494 420 350 300 249 200 151] * 1e7;

  ttc102=[598 613 621 625 635 645 655 663 676 ...
          687 697 703 712 720 728 734 759 768 ...
          780 792  ]/1000;
  nub102=[2982 2895 2837 2760 2698 2613 2528 2445 2358 ...
          2263 2210 2127 2083 2003 1947 1877 1708 1618 ...
          1515 1420  ] * 1e7;

  ttc155=[475 480 483 487 495 506 517 527 540 ...
          553 568 584 602 617 634 654 675 704 ...
          725  ]/1000; %add more
  nub155=[5167 5144 5075 5104 5043 4967 4890 4797 4691 ...
          4541 4407 4227 3994 3872 3659 3426 3143 2835 ...
          2572  ] * 1e7; %add more

  ttc250=[502 514 521 528 540 555 578 600 618 ...
          647 670 687 718 752 782 826 846 877 ...
          896 907 917 930]/1000;
  nub250=[6695 6635 6522 6398 6300 6070 5820 5615 5275 ...
          4895 4604 4326 3791 3223 2833 2305 1925 1513 ...
          1220 1095 985 790] * 1e7;

  % our data - HPD oscillations
  % 19.5 bar data?

  plot(ttc, he3_nu_b(ttc,0.5).^2, 'r')
  plot(ttc, he3_nu_b(ttc,5).^2, 'g')
  plot(ttc, he3_nu_b(ttc,10.2).^2, 'b')
  plot(ttc, he3_nu_b(ttc,15.5).^2, 'm')
  plot(ttc, he3_nu_b(ttc,25.0).^2, 'c')

  xlim([0 1]);
  ylim([0 1e+11]);
  legend('0.5 bar', '5 bar', '10.2 bar', '15.5 bar', '25 bar')

  plot(ttc005, nub005, 'ro')
  plot(ttc050, nub050, 'go')
  plot(ttc102, nub102, 'bo')
  plot(ttc155, nub155, 'mo')
  plot(ttc250, nub250, 'co')

print gap.eps -deps -color

