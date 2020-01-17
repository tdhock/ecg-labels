HOCKING-ecg.pdf: HOCKING-ecg.tex figure-two-ecg-graphs.tex figure-one-ecg-graph.tex
	pdflatex HOCKING-ecg
figure-two-ecg-graphs.tex: figure-two-ecg-graphs.R
	R --vanilla < $<
figure-one-ecg-graph.tex: figure-one-ecg-graph.R
	R --vanilla < $<
