CFLAGS += $$(gsl-config --cflags) -std=c99 -Wall
LDLIBS += $$(gsl-config --libs)
OBJ = main.o integration.o
.PHONEY: all
all: out.txt
	cat $<

out.txt: main Makefile
	./$< 1>$@

main: $(OBJ)

$(OBJ): adaptive.h

clean: ; $(RM)\
		main *.o \
		Integration *.o\
