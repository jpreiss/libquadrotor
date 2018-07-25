libquadtest: src/*.c include/*.h
	mkdir -p bin
	gcc -std=c99 -Iinclude -Iinclude/cmath3d -Wdouble-promotion -g -Wall -Wpedantic -o bin/libquadtest src/libquadrotor_test.c -lm

test: libquadtest
	./bin/libquadtest

debug: libquadtest
	gdb -tui ./bin/libquadtest

clean:
	rm -rf bin/*
