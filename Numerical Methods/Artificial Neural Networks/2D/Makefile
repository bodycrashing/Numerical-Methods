CFLAGS += -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)
OBJ = main.o ann.o

Interpolation.svg: plot.gpi out.txt Makefile
	gnuplot $<

out.txt:main
	./$< 1> $@ 2>log

main: $(OBJ)

$(OBJ): ann.h

plot.gpi: Makefile
	echo 'set term svg size 800,600 background rgb "white"' > $@
	echo 'set out  "Neural_Network_Interpolation.svg"'>> $@
	echo 'set xlabel "x"' >> $@
	echo 'set ylabel "y"' >> $@
	echo 'set key right bottom font ",10"' >> $@
	echo 'set grid' >> $@
	echo 'set view 70,50,1,1' >> $@
	echo 'set title  "Neural Network Interpolation of cos(5*x-1)*exp(-y*y)*exp(-x*x) "' >> $@
	echo 'splot cos(5*x-1)*exp(-y*y)*exp(-x*x) title "Data","out.txt" index 1 with points pointtype 6 title "Neural Interpolation" ' >> $@

clean:
		rm -f main *.o
