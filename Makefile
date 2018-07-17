all:
	@pdflatex 2019_01-Zwanenburg_PhD_Thesis.tex
	@bibtex 2019_01-Zwanenburg_PhD_Thesis
	@pdflatex 2019_01-Zwanenburg_PhD_Thesis.tex
	@pdflatex 2019_01-Zwanenburg_PhD_Thesis.tex

pdf:
	@pdflatex 2019_01-Zwanenburg_PhD_Thesis.tex
