all:
clean:
	rm -f *.mexglx *.mex

octave:
	rm -f *.mex
	octave -q --eval compile

matlab:
	rm -f *.mexglx
	matlab -nojvm -nosplash -r 'compile;exit'

matlab64:
	matlab64 -nojvm -nosplash -r 'compile; exit'
