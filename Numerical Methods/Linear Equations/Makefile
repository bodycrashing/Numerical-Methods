CFLAGS += -Wall -O -std=gnu1x `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

all: out.txt
	cat $<

out.txt: main
	./$< > $@

main: main.o qr_gs_decomp.o
main.o qr_gs_decomp.o:qr.h

clean:
	find . -type f -executable -delete
	rm -f *.o out*
