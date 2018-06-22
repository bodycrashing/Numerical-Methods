CFLAGS +=  -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS += -lm $$(gsl-config --libs)
OBJ = main.o minimization.o qr_gs_decomp.o
plot_obj = fit_data.txt
.PHONEY: all

all: out.txt plot.svg
	cat $<
	touch $^


plot.svg: $(plot_obj)
	echo ' \
	set term $(subst .,,$(suffix $@)) background rgb "white"; \
	set out "$@"; \
	set grid; \
	set xlabel "t"; \
	set ylabel "y"; \
	set title "Least Squares Fit Using Quasi-Newton Minimization"; \
	plot "$<" index 0 with errorbars title "data" \
	,"$<" index 1 using 1:2 with lines title "fit" ; \
	' | gnuplot

out.txt: main
	./$< >$@

main: $(OBJ)

$(OBJ): minimization.h

clean:
	$(RM) *.o main *.txt *.svg
