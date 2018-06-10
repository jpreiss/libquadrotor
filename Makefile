libquadtest: src/*.c include/*.h
	mkdir -p bin
	gcc -std=c11 -Iinclude -Wdouble-promotion -g -Wall -Wpedantic -o bin/libquadtest src/test.c -lm

test: libquadtest
	./bin/libquadtest

debug: libquadtest
	gdb -tui ./bin/libquadtest
