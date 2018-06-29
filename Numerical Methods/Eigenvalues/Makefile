CFLAGS +=  $$(gsl-config --cflags) -Wall -std=gnu1x
LDLIBS += $$(gsl-config --libs)
OBJ = jacobi.o main_cyclic.o main_eig_by_eig.o

all: out.txt times.svg
	cat $<

N=5
out.txt: main Makefile
	./main $(N) 1> $@

main: main.o jacobi.o


times.svg: out.times.cyclic.txt out.times.eigbyeig.txt
	echo '\
	set term svg background rgb "white"; \
	set out "$@"; \
	set xlabel "matrix size"; \
	set ylabel "diagonalization time, sec"; \
	set title "diagonalization time as function of matrix size (`uname -m`)"; \
	set key left; \
	f1(x)=b1+(x/a1)**c1; \
	b1=0; a1=100; c1=3;\
	fit f1(x) "out.times.cyclic.txt" via a1,c1; \
	T1 = sprintf("fit: (n/%3.0f)^{%3.1f}",a1,c1); \
	f2(x)=b2+(x/a2)**c2; \
	b2=0; a2=100; c2=3;\
	fit f2(x) "out.times.eigbyeig.txt" via a2,c2; \
	T2 = sprintf("fit: (n/%3.0f)^{%3.1f}",a2,c2); \
	plot "out.times.cyclic.txt" title "measurement cyclic", f1(x) title T1 \
		,"out.times.eigbyeig.txt" title "measurement cyclic", f2(x) title T2; \
	'|tee log.gpi|gnuplot 2> /dev/null



out.times.cyclic.txt: main_cyclic Makefile
	cat /dev/null > $@
	for n in `seq 50 10 60`; do time --format "$$n %U" --append --output $@ ./main_cyclic $$n >/dev/null; done

main_cyclic: jacobi.o main_cyclic.o



out.times.eigbyeig.txt: main_eig_by_eig Makefile
	cat /dev/null > $@
	for n in `seq 50 10 60`; do time --format "$$n %U" --append --output $@ ./main_eig_by_eig $$n >/dev/null; done

main_eig_by_eig: jacobi.o main_eig_by_eig.o


$(OBJ): jacobi.h


clean:
	rm -f main *.o out* *.png *.svg log* *log *.txt
