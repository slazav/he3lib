#!/usr/bin/octave -qf

addpath ../matlab

  figure; clf; hold on;
  subplot(2,1,1); hold on; title('gap vs ttc')

    ttc=[0.001:0.001:1];
    leg={};

    plot(ttc, he3_bcsgap(ttc), 'k-')
    leg{end+1}=sprintf('BCS');

    c='rgbm';
    p=[0 10 20 30];
    for i=1:length(p);
      plot(ttc, he3_trivgap(ttc,p(i)), [ c(i) '-']);
      leg{end+1}=sprintf('WC+ %d bar', p(i));
    end

    legend(leg, 'location', 'southwest')
    xlim([0 1]);
    ylim([0 2]);
    set(gca,'yTick', 0:0.2:2);
    grid on;

  subplot(2,1,2); hold on;

    p=0:1:34;
    leg={};

    c='rgbm';
    ttc=[0 0.5];
    for i=1:length(ttc);
      plot(p, he3_bcsgap(ttc(i))*ones(size(p)), 'k-')
      plot(p, he3_trivgap(ttc(i),p), 'r-');
      plot(p, he3_todogap(ttc(i),p), 'b-');
      leg{end+1}=sprintf('BCS');
      leg{end+1}=sprintf(' WC+ %d bar', p(i));
      leg{end+1}=sprintf('TODO %d bar', p(i));
    end

    legend(leg, 'location', 'northwest')
%    xlim([0 1]);
%    ylim([0 2]);
%    set(gca,'yTick', 0:0.2:2);
    grid on;


  print gap1.eps -deps "-S740,1160" -color

