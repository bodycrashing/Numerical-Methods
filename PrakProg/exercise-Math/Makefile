CLFLAGS = -Wall -std=gnu99
CC = gcc
LDFLAGS = -lm


.PHONEY: all

all: out.txt
	cat $<

out.txt: Math
	./$< > $@


.PHONEY: clean
clean:
	$(RM) -f math mathout.o
