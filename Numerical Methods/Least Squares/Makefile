CFLAGS += -Wall -O -std=gnu1x `gsl-config --cflags`
LDLIBS += `gsl-config --libs`
OBJ = main.o qr_gs_decomp.o least_squares.o jacobi.o

all: out.txt plot.svg

plot.svg: out.txt Makefile
	echo '\
	set term svg background rgb "white";\
	set out "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set tics out;\
	set key out;\
	set title "least-squares fit of F_c(x)=c_0/x+c_1+c_2*x to some data" ;\
	plot "$<" index 0 with errorbars title "data"\
	,"$<" index 1 using 1:2 with lines title "F_c(x)"\
	,"$<" index 1 using 1:3 with lines title "F_C(x)_{c0+dc0}"\
	,"$<" index 1 using 1:4 with lines title "F_C(x)_{c0-dc0}"\
	,"$<" index 2 using 1:3 with lines title "F_C(x)_{c1+dc1}"\
	,"$<" index 2 using 1:4 with lines title "F_C(x)_{c1-dc1}"\
	,"$<" index 3 using 1:3 with lines title "F_C(x)_{c2+dc2}"\
	,"$<" index 3 using 1:4 with lines title "F_C(x)_{c2-dc2}"\
	'|gnuplot

out.txt: main
	./$< > $@

main: $(OBJ)
$(OBJ): qr.h least_squares.h

clean:
	find . -type f -executable -delete
	rm -f *.o out*
