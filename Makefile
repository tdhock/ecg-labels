HOCKING-ecg.pdf: HOCKING-ecg.tex figure-two-ecg-graphs.tex
	pdflatex HOCKING-ecg
figure-two-ecg-graphs.tex: figure-two-ecg-graphs.R
	R --vanilla < $<
