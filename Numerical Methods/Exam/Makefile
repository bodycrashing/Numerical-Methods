CFLAGS += $$(gsl-config --cflags) -std=gnu1x -Wall
LDLIBS += $$(gsl-config --libs)

OBJ = main.o rk5.o

.PHONEY: all

all: out.txt harmonic.svg
	cat $<

harmonic.svg: plot.gpi data.txt
	gnuplot $<

 out.txt data.txt: main
	./$< > $@

main: $(OBJ)

$(OBJ): rk5.h

clean:
	rm --force *.o *.svg *.txt main data
