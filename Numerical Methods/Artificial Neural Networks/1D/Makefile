CFLAGS += $$(gsl-config --cflags) -std=gnu1x -Wall
LDLIBS += $$(gsl-config --libs)
obj =  ann.o minimization.o qr.o minimization.h
datafile = plot_data.txt

.PHONEY: all clean

all: out.txt plot.svg
	cat $<
	touch plot.svg


plot.svg: $(datafile)
	@echo ' \
    	set term $(subst .,,$(suffix $@)) background rgb "white"; \
        set out "$@"; \
        set grid; \
        set xlabel "x"; \
        set ylabel "y"; \
        set title "ANN interpolation"; \
        plot "$<" index 0 with points title "cos(x)sin(2x) data points" \
        ,"$<" index 1 with lines title "ANN interpolation" ; \
        ' | gnuplot


out.txt: main
	./$< > $@

main: $(obj)

clean:
	$(RM) *.o main *.txt *.svg *.log
