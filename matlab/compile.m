function compile()
% Parse he3.h file and compile matlab interfaces
% for all constants and functions.

  f=fopen('../he3.h');
  while ~feof(f)
    str=fgets(f);

    % constants: extern double name_;
    s=regexp(str,...
      '^\s*extern\s+double\s+(\w+)_\s*;',...
      'tokens','once');
    if length(s)>0; comp(s{1}, 0); end

    % functions: double name_(double *arg1, ...);
    % 1..8 double* aguments
    a1='(\s*double\s+\*\w+\s*)';
    a2='(,\s*double\s+\*\w+\s*)?';
    s=regexp(str, ...
      ['^\s*double\s+(\w+)_\(' a1 a2 a2 a2 a2 a2 a2 a2 '\)\s*;'],...
      'tokens', 'once');

    % for n<nmax arguments octave have less cells in the s
    % and matlab have empty cells in the end of s
    narg=0;
    for i=2:length(s); if length(s{i}); narg=i-1; end; end
    if narg>0; comp(s{1}, narg); end

  end
  fclose(f);
end

function comp(name, narg)
  fprintf('>> compiling %s (%d args)\n', name, narg);
  if iscell(name); name=name{1}; end

  if length(regexp(version, '^3\.\d+\.\d+'))
    % old matlab + octave
    mex(['-DFUNC=' name '_ -DNARGIN=' num2str(narg)],...
         '-o', name, '-lhe3',...
        ['-L' pwd ], ['-Wl,--rpath=' pwd ],...
        'mexfunc.c');
    return
  end
  if length(regexp(version, '^8\.')) || length(regexp(version, '^9\.'))
    % new matlab
    mex(['-DFUNC=' name '_'], ['-DNARGIN=' num2str(narg)],...
        '-output', name, '-lhe3',...
        ['-L' pwd ], ['LDFLAGS="$LDFLAGS -Wl,-rpath=' pwd '"'],...
        'mexfunc.c');
    return
  end
  fprintf('Unknown matlab version: %s\n', version);
  exit
end



