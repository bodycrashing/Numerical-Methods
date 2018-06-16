CFLAGS += $$(gsl-config --cflags) -std=gnu1x -Wall
LDLIBS += $$(gsl-config --libs)
obj = newton.o qr_gs_decomp.o root.h

.PHONEY: test all clean

all: out.txt
	cat $<

out.txt: main
	./$< > $@

main: $(obj)

clean:
	$(RM) *.o main *.txt

test:
	@echo $(LDLIBS)
	@echo $(CFLAGS)
