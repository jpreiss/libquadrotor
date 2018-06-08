libquadtest: src/libquadrotorcontrol.c include/libquadrotorcontrol.h include/cmath3d/math3d.h
	mkdir -p bin
	gcc -std=c11 -Iinclude -Wdouble-promotion -DTEST -g -Wall -Wpedantic -o bin/libquadtest src/libquadrotorcontrol.c -lm

test: libquadtest
	./bin/libquadtest

debug:
	gdb -tui ./bin/libquadtest
