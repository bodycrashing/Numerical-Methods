CFLAGS += -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)

.PHONEY:all

all:airy.svg out.txt

out.txt: main
	./$< > $@
	cat $@

airy.svg:airy.gpi airy.out.txt
	gnuplot $<

airy.out.txt:airy
	./$< > $@

.PHONEY:clean
clean:
	find . -type f -executable -delete
	rm -f *.o *.txt *.data *.svg
