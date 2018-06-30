CFLAGS = $$(gsl-config --cflags) -std=gnu1x -Wall
LDLIBS = $$(gsl-config --libs)

.PHONEY: clean

plot.svg: plot.gpi out.txt
	gnuplot $<

out.txt: main
	./$< > $@

main: F.o

clean:
	$(RM) *.o main *.svg *.txt
