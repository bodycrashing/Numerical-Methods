CFLAGS +=  $$(gsl-config --cflags) -Wall -std=gnu1x
LDLIBS += $$(gsl-config --libs)


all: out.txt times.svg
	cat $<

N=5
out.txt: main Makefile
	./main $(N) 1> $@

main: main.o jacobi.o


times.svg: out.times.txt
	echo '\
	set term svg background rgb "white"; \
	set out "$@"; \
	set xlabel "matrix size"; \
	set ylabel "diagonalization time, sec"; \
	set title "diagonalization time as function of matrix size (`uname -m`)"; \
	set key left; \
	f(x)=b+(x/a)**c; \
	b=0; a=100; c=3;\
	fit f(x) "$<" via a,c; \
	T = sprintf("fit: (n/%3.0f)^{%3.1f}",a,c); \
	plot "$<" title "measurement", f(x) title T; \
	'|tee log.gpi|gnuplot 2> /dev/null


out.times.txt: main Makefile
	cat /dev/null > $@
	for n in `seq 100 15 150`; do time --format "$$n %U" --append --output $@ ./main $$n >/dev/null; done

clean:
	rm -f main *.o out* *.png *.svg log* *log
