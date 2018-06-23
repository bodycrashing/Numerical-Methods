CFLAGS +=  -Wall -std=gnu1x $$(gsl-config --cflags) -fopenmp
LDLIBS += -lm -lgomp $$(gsl-config --libs)
OBJ =  smp.o main.o


.PHONEY: all clean

all: out.txt
	cat $<

out.txt: main
	./$< > $@

main: $(OBJ)

$(OBJ): smp.h


clean:
	$(RM) *.o main *.txt *.svg *.log
