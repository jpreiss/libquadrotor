libquadtest: src/*.c include/*.h
	mkdir -p bin
	gcc -std=c99 -Iinclude -Wdouble-promotion -g -Wall -Wpedantic -o bin/libquadtest src/test.c -lm

test: libquadtest
	./bin/libquadtest

debug: libquadtest
	gdb -tui ./bin/libquadtest

clean:
	rm -rf bin/*
