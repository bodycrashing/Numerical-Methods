CFLAGS = -Wall -std=gnu99 `gsl-config --cflags`
LDLIBS = `gsl-config --libs`

.PHONEY:all

all:integ.gpi integ.out.txt
	gnuplot $<

integ.out.txt:integ
	./$< > $@

.PHONEY: disp
disp:integ.svg
	display $< &

.PHONEY:clean
clean:
	find . -type f -executable -delete
	rm -f *.o *.txt *.data *.svg
