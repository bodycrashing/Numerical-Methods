CFLAGS +=  $$(gsl-config --cflags) -Wall -std=gnu1x
LDLIBS += $$(gsl-config --libs)
OBJ = jacobi.o main_cyclic.o main_eig_by_eig.o

all: out.txt times.cyclic.svg times.eigbyeig.svg
	cat $<

N=5
out.txt: main Makefile
	./main $(N) 1> $@

main: main.o jacobi.o


times.cyclic.svg: out.times.cyclic.txt
	echo '\
	set term svg background rgb "white"; \
	set out "$@"; \
	set xlabel "matrix size"; \
	set ylabel "diagonalization time, sec"; \
	set title "diagonalization time as function of matrix size (`uname -m`)"; \
	set key left; \
	f(x)=b+(x/a)**c; \
	b=0; a=100; c=3;\
	fit f(x) "out.times.cyclic.txt" via a,c; \
	T = sprintf("fit: (n/%3.0f)^{%3.1f}",a,c); \
	plot "out.times.cyclic.txt" title "measurement cyclic", f(x) title T \
	'|tee log.gpi|gnuplot 2> /dev/null


times.eigbyeig.svg: out.times.eigbyeig.txt
	echo '\
	set term svg background rgb "white"; \
	set out "$@"; \
	set xlabel "matrix size"; \
	set ylabel "diagonalization time, sec"; \
	set title "diagonalization time as function of matrix size (`uname -m`)"; \
	set key left; \
	f1(x)=b1+(x/a1)**c1; \
	b1=0; a1=100; c1=3;\
	fit f1(x) "out.times.eigbyeig.txt" via a1,c1; \
	T1 = sprintf("fit: (n/%3.0f)^{%3.1f}",a1,c1); \
	plot "out.times.eigbyeig.txt" title "measurement eig by eig", f1(x) title T1 \
	'|tee log.gpi|gnuplot 2> /dev/null



out.times.cyclic.txt: main_cyclic Makefile
	cat /dev/null > $@
	for n in `seq 50 10 150`; do time --format "$$n %U" --append --output $@ ./main_cyclic $$n >/dev/null; done

main_cyclic: jacobi.o main_cyclic.o



out.times.eigbyeig.txt: main_eig_by_eig Makefile
	cat /dev/null > $@
	for n in `seq 30 10 130`; do time --format "$$n %U" --append --output $@ ./main_eig_by_eig $$n >/dev/null; done

main_eig_by_eig: jacobi.o main_eig_by_eig.o


$(OBJ): jacobi.h


clean:
	rm -f main *.o out* *.png *.svg log* *log *.txt
