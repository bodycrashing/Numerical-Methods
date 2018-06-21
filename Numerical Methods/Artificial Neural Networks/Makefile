CFLAGS += $$(gsl-config --cflags) -std=gnu1x -Wall
LDLIBS += $$(gsl-config --libs)
obj =  ann.c minimization.o qr.o minimization.h
datafile = plot_data.txt

.PHONEY: test all clean

all: out.txt plot.svg plot2.svg
	cat $<
	touch plot.svg plot2.svg

plot2.svg: $(datafile)
	@echo ' \
        set term $(subst .,,$(suffix $@)) background rgb "white"; \
        set out "$@"; \
        set grid; \
        set xlabel "x"; \
        set ylabel "y"; \
        set title "ANN interpolation"; \
        splot "$<" index 2 with points title " cos(x+y)*sin(2*(x+y)) data points" \
	,"$<" index 3  using 1:2:3 with point title "ANN interpolation" ; \
        ' | gnuplot


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

test:
	@echo $(LDLIBS)
	@echo $(CFLAGS)
