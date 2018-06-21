CFLAGS += $$(gsl-config --cflags) -std=gnu1x -Wall
LDLIBS += $$(gsl-config --libs)

OBJ = main.o runge_kutta.o

.PHONEY: all
all: out.txt

out.txt: main
	./$< > $@

main: $(OBJ)

$(OBJ): runge_kutta.h

clean:
	rm --force *.o *.svg *.txt main data
