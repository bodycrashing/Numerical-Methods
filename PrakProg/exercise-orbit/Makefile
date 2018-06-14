CFLAGS = -Wall -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`
PLOTTER = gnuplot

.PHONEY:default

default:logistic.svg

logistic.svg:logistic.gpi logistic.out.txt
	gnuplot $<

logistic.out.txt:orbit
	./$< > $@

.PHONEY:clean
clean:
	find . -type f -executable -delete
	rm -f *.o *.txt *.data *.svg
