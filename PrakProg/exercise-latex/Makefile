CFLAGS = -Wall -std=gnu99 `gsl-config --cflags`
LDLIBS = `gsl-config --libs`
CC = gcc


report.pdf: report.tex figure.pdf
	pdflatex $<
	pdflatex $<

figure.pdf: plot.gpi data
	gnuplot plot.gpi

data: main
	./main -3 3 0.1 > data

.PHONEY:clean
clean:
	$(RM) main data *.pdf *.log *.aux *.dvi *.out
