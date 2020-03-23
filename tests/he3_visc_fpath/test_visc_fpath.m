function test_gap()

  find_figure('test_visc_fpath')
  clf; hold on;

  pp = [00 10 20 30];
  cols = 'rgbm';
  ttc = 0.3:0.003:1;
  for i = 1:length(pp)
    ff = fopen(sprintf('1997_einzel_parp_%02d.txt', pp(i)));
    data = textscan(ff, '%f %f');
    fclose(ff);
    plot(data{1}, data{2}, [cols(i) '*'])

    mfp = he3_visc_fpath(ttc,pp(i))*10000 % cm -> um
    plot(ttc, mfp, [cols(i) '-'])
  end
  xlim([0.2 1.0])
  ylim([0 100])
end

