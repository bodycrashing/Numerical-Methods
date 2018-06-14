CFLAGS = -Wall -std=gnu1x -O2
CFLAGS += $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)
PLOTTER = gnuplot
OBJ = main.o linspace.o linterp.o qspline.o cspline.o

.PHONEY:all clean lindisp qdisp

all: linterp.svg qspline.svg cspline.svg data.out.txt

linterp.svg: linterp.gpi data.out.txt
		$(PLOTTER) $<

qspline.svg: qspline.gpi data.out.txt
		$(PLOTTER) $<

cspline.svg: cspline.gpi data.out.txt
		$(PLOTTER) $<

data.out.txt: main
	./$< > $@

main: $(OBJ) interp.h


lindisp:linterp.svg
	display $< &

qdisp:qspline.svg
	display $< &

clean:
	find . -type f -executable -delete
	rm -f *.o *.txt *.data *.svg
