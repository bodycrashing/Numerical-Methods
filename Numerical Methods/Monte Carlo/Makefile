CFLAGS += $$(gsl-config --cflags) -std=gnu1x -Wall
LDLIBS += $$(gsl-config --libs)
OBJ =  Monte_Carlo.o main.o integration.o


.PHONEY: test all clean

all: plot.gpi out.txt error.svg
	gnuplot $<

error.svg: main data.txt

data.txt: main

out.txt: main
	./$< > $@
	cat $@

main: $(OBJ)

$(OBJ): Monte_Carlo.h

plot.gpi: Makefile
		echo 'set term svg size 800,600 background rgb "white"' > $@
		echo 'set out  "error.svg"'>> $@
		echo 'set xlabel "Number of sampling points"' >> $@
		echo 'set ylabel "Error"' >> $@
		echo 'set key right top font ",14"' >> $@
		echo 'set title  "Error of Plain Monte Carlo Integration as Function of N"' >> $@
		echo 'f(x)=b+(a/sqrt(x))' >>$@
		echo 'b=0.1; a=10' >>$@
		echo 'fit f(x) "data.txt" via b,a' >>$@
		echo 'F = sprintf("fit: %4.1f+(%3.0f/sqrt(N))",b,a)' >>$@
		echo 'plot "data.txt" using 1:2 with points  title "Error vs. N",f(x) with lines title F \' >> $@ 


clean:
	$(RM) *.o main *.txt *.svg *.log

test:
	@echo $(LDLIBS)
	@echo $(CFLAGS)
