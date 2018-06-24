CFLAGS = -Wall -std=gnu99
LDLIBS = -lm
SCRIPT = log.gpi
FORMAT = svg
PLOTTER = gnuplot

.PHONEY: all
all: plot.$(FORMAT)

plot.$(FORMAT) : out-stdio.txt out-cmdline.txt Makefile
	echo 'set term $(FORMAT) background rgb "white"' > $(SCRIPT)
	echo 'set out "$@"' >> $(SCRIPT)
	echo 'set key out' >> $(SCRIPT)
	echo 'set xlabel "x"' >> $(SCRIPT)
	echo 'set ylabel "y"' >> $(SCRIPT)
	echo 'plot \' >> $(SCRIPT)
	echo '"$(word 1,$^)" with points title "stdio" \' >> $(SCRIPT)
	echo ',"$(word 2,$^)" with points title "cmdline" \' >> $(SCRIPT)
	echo ',sin(x) with lines \' >> $(SCRIPT)
	echo ',cos(x) with lines' >> $(SCRIPT)
	$(PLOTTER) $(SCRIPT)

out-stdio.txt: in.txt main-stdio
	./$(word 2,$^) < $< > $@

out-cmdline.txt: in.txt main-cmdline
	./$(word 2,$^) `cat $<` > $@

in.txt: Makefile
	seq 1 0.5 4 > in.txt


.PHONEY: clean
clean :
	find ./ -type f -executable -delete
	$(RM) *.txt $(SCRIPT) *.svg *.png *.pdf log*

.PHONEY: test
test:
	find ./ -type f -executable -exec echo {} \;
