CFLAGS = -Wall -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`
PLOTTER = gnuplot

.PHONEY:default

default:logistic.svg  plot1.svg plot2.svg plot3.svg

logistic.svg:logistic.gpi logistic.out.txt
	gnuplot $<

logistic.out.txt: main_logistic
	./$< > $@

main_logistic: main_logistic.o

plot1.svg: plot1.gpi out1.txt
	gnuplot $<

plot2.svg: plot2.gpi out2.txt
	gnuplot $<

plot3.svg: plot3.gpi out3.txt
	gnuplot $<

out1.txt: main Makefile
	./$< 0 0 > $@

out2.txt: main Makefile
	./$< 0 -0.4 > $@

out3.txt: main Makefile
	./$< 0.01 -0.4 > $@

main: main.o

.PHONEY:clean
clean:
	find . -type f -executable -delete
	rm -f *.o *.txt *.data *.svg
