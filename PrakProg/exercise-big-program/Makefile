CFLAGS = -Wall -std=gnu99 `gsl-config --cflags`
LDLIBS = `gsl-config --libs`
CC = gcc
OBJ = main.o log_int.o


.PHONEY:clean all

all: out.txt
	cat $<




out.txt: main
	./$< > $@


main: $(OBJ)

$(OBJ): log_int.h

clean:
	$(RM) main data *.pdf *.log *.aux *.dvi *.out
