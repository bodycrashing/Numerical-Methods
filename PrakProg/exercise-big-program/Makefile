CFLAGS = -Wall -std=gnu99 `gsl-config --cflags`
LDLIBS = `gsl-config --libs`
CC = gcc
OBJ = main.o log_int.o


.PHONEY:clean all

report.pdf: report.tex out.txt figure.pdf
	pdflatex $<
	pdflatex $<


out.txt: main
	./$< > $@


figure.pdf: plot.gpi out.txt
		gnuplot plot.gpi

plot.gpi:Makefile
	echo 'set term pdf' >$@
	echo 'set out "figure.pdf"' >> $@
	echo 'set title  "log(x) integral definition"' >> $@
	echo 'set xlabel "x"' >> $@
	echo 'set ylabel "data"' >> $@
	echo 'set grid' >> $@
	echo 'set key left' >> $@
	echo 'plot \' >>$@
	echo '"out.txt" index 0 using 1:2 with points pointtype 6 title "log(x) eq. 1.0.1"\
				,"out.txt" index 0 using 1:3 with lines title "math.h log(x)";\' >>$@

main: $(OBJ)

$(OBJ): log_int.h

clean:
	$(RM) main data *.pdf *.log *.aux *.dvi *.out *.svg
