CFLAGS = -Wall -std=gnu99 -O1
LSLIBS = -lm
CC = gcc



default: out.txt
	cat $<

out.txt: exercise2
	./$< > $@

exercise2: exercise2.c

clean:
	rm -f out.txt exercise2 exercise2.o
