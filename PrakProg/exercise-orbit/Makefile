CFLAGS = -Wall -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`
PLOTTER = gnuplot

.PHONEY:default

default:logistic.svg out.txt out1.txt out2.txt

logistic.svg:logistic.gpi main
	gnuplot $<

out.txt: main Makefile
	./$< 0 0 > $@

out1.txt: main Makefile
	./$< 0 -0.4 > $@

out2.txt: main Makefile
	./$< 0.01 -0.4 > $@


orbit1.svg: orbit1.gpi orbit1.txt
	gnuplot $<

orbit2.svg: orbit2.gpi orbit2.txt
	gnuplot $<

orbit3.svg: orbit3.gpi orbit3.txt
	gnuplot $<



.PHONEY:clean
clean:
	find . -type f -executable -delete
	rm -f *.o *.txt *.data *.svg
