CFLAGS += -Wall -std=c11 `gsl-config --cflags`
LDLIBS += `gsl-config --libs`
OBJ = main.o newton.o qr_gs_decomp.o

all:out.txt
	cat $<

out.txt: main
	./$< > $@

main: $(OBJ)

$(OBJ): qr.h

clean:
	rm -f main *.o out*
q
