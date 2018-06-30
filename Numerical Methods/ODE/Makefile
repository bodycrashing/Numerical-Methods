CFLAGS += $$(gsl-config --cflags) -std=gnu1x -Wall
LDLIBS += $$(gsl-config --libs)

OBJ = main.o runge_kutta.o

.PHONEY: all

all: out.txt harmonic.svg logistic.svg
	cat $<

harmonic.svg: plot.harmonic.gpi data.txt
	gnuplot $<

logistic.svg: plot.logistic.gpi data.txt
	gnuplot $<

 out.txt data.txt: main
	./$< > $@

main: $(OBJ)



$(OBJ): runge_kutta.h

clean:
	rm --force *.o *.svg *.txt main data
