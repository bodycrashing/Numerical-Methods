CFLAGS = -Wall -std=c99
LDLIBS = -lm
OBJ = main.o nvector.o


.PHONEY: all clean

all: out.txt
	cat $<

out.txt: main
	./$< > $@

main: $(OBJ)

$(OBJ): nvector.h

clean:
	-rm -f out.txt main main.o nvector.o
