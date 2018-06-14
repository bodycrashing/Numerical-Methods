CFLAGS += -Wall -std=c11 `gsl-config --cflags`
LDLIBS += `gsl-config --libs`
OBJ = main.o runge_kutta.o

.PHONEY: all
all: out.txt
	cat $<

out.txt: main
	./$< > $@

main: $(OBJ)

$(OBJ): runge_kutta.h

clean:
	rm --force *.o *.svg *.txt main data
