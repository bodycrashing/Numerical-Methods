CFLAGS = -Wall -std=gnu99 -O1
LSLIBS = -lm
CC = gcc
OBJ = main.o equal.o


.PHONEY: default clean

default: out.txt
	cat $<

out.txt: main
	./$< > $@

main: $(OBJ)

clean:
	$(RM) *.o main out.txt
