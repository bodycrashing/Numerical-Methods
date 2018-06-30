CFLAGS = -Wall -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`
PLOTTER = gnuplot

.PHONEY:default

default:logistic.svg out1.txt out2.txt out3.txt plot1.svg plot2.svg plot3.svg

logistic.svg:logistic.gpi logistic.out.txt
	gnuplot $<

logistic.out.txt: main
	./$< > $@


out1.txt: main Makefile
	./$< 0 0 > $@

out2.txt: main Makefile
	./$< 0 -0.4 > $@

out3.txt: main Makefile
	./$< 0.01 -0.4 > $@

plot1.svg: plot1.gpi out1.txt
	gnuplot $<

plot2.svg: plot2.gpi out2.txt
	gnuplot $<

plot3.svg: plot3.gpi out3.txt
	gnuplot $<



.PHONEY:clean
clean:
	find . -type f -executable -delete
	rm -f *.o *.txt *.data *.svg
