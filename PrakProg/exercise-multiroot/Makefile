CFLAGS = $$(gsl-config --cflags) -std=gnu1x -Wall
LDLIBS = $$(gsl-config --libs)

plot.svg: plot.gpi out.txt
	gnuplot $<

out.txt: main
	./$< > $@

main: F.o

.PHONEY: clean test
clean:
	$(RM) *.o main


test:
	echo $(CFLAGS)
	echo $(LDLIBS)
