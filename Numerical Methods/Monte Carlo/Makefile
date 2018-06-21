CFLAGS += $$(gsl-config --cflags) -std=gnu1x -Wall
LDLIBS += $$(gsl-config --libs)
obj =  Monte_Carlo.o main.o integration.o Monte_Carlo.h


.PHONEY: test all clean

all: plot.gpi out.txt error.svg
	gnuplot $<

error.svg: main data.txt

data.txt: main

out.txt: main
	./$< > $@
	cat $@

main: $(obj)


plot.gpi: Makefile
		echo 'set term svg size 800,600 background rgb "white"' > $@
		echo 'set out  "error.svg"'>> $@
		echo 'set xlabel "N"' >> $@
		echo 'set ylabel "Error"' >> $@
		echo 'set key right top font ",12"' >> $@
		echo 'set title  "Error of Plain Monte Carlo Integration on Func. x0*x0+x1*x1+x2*x2"' >> $@
		echo 'f(x)=b+(a/sqrt(x))' >>$@
		echo 'b=0.1; a=10' >>$@
		echo 'fit f(x) "data.txt" via b,a' >>$@
		echo 'F = sprintf("fit: error=%4.1f+(%3.0f/sqrt(N))",b,a)' >>$@
		echo 'plot "data.txt" using 1:2 with points  title "Data of Error vs. N",f(x) with lines title F \' >> $@ \


clean:
	$(RM) *.o main *.txt *.svg *.log

test:
	@echo $(LDLIBS)
	@echo $(CFLAGS)
