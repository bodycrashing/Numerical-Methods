CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

fit.svg: plot.gpi out.txt
	gnuplot $<

out.txt: main
	./$< > $@

main:

.PHONEY: clean test disp
disp: fit.svg
	display $^ &


test:
	echo $(LDLIBS)
	echo $(CFLAGS)

clean:
	$(RM) main *.o
