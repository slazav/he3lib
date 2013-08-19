all: he3lib.pdf spin_diff.pdf

he3lib.dvi: ptable.tex he3lib.tex

ptable.tex: ptable.m
	octave -q $< || rm -f ptable.tex

%.pdf: %.ps
	ps2pdf $<
%.ps: %.dvi
	dvips $<
%.dvi: %.tex
	make -C pics
	latex $<

clean:
	rm -f *.aux *.log *.nav *.out *.snm *.toc *.dvi *.eps
	rm -f ptable.tex
