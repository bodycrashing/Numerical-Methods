.PHONEY: all clean

all:test_cmd.txt test_stdio.txt plot
	cat test.io.out.txt
	cat test.cmd.out.txt

plot:input.data.txt out.txt
	gnuplot plot.gpi

out.txt:main_stdio main_cmdline
	cat input.data.txt | ./main_stdio > out.io.txt
	./main_cmdline `cat input.data.txt` > out.cmd.txt

input.data.txt:
	echo `seq 0 0.2 5` > input.data.txt

test_cmd:test_cmd.txt
	cat test.cmd.out.txt

test_stdio:test_stdio.txt
	cat test.io.out.txt

test_stdio.txt:main_stdio
	echo 1 2 3 4 5 | ./main_stdio -> test.io.out.txt
	echo `seq 1 0.5 4` | ./main_stdio >> test.io.out.txt

test_cmd.txt:main_cmdline
	./main_cmdline 1 2 3 4 5 6 > test.cmd.out.txt
	./main_cmdline `seq 0 0.5 2` >> test.cmd.out.txt


main_stdio:
	gcc main_stdio.c -o main_stdio -lm


main_cmdline:
	gcc main_cmdline.c -o main_cmdline -lm

clean:
	$(RM) main_stdio main_stdio.o out.txt test.io.out.txt test.cmd.out.txt \
	 main_cmdline main_cmdline.o input.data.txt out.cmd.txt out.io.txt plot.svg
