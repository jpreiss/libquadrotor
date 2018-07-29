libquadtest: src/*.c src/polysolve.c include/*.h
	mkdir -p bin
	gcc -std=c99 -Iinclude -Iinclude/cmath3d -Wdouble-promotion -g -Wall -Wpedantic -o bin/libquadtest src/libquadrotor_test.c -lm

src/polysolve.c: codegen/polynomial_boundary_solver.py
	python3 codegen/polynomial_boundary_solver.py 7 > src/polysolve.c

test: libquadtest
	./bin/libquadtest

debug: libquadtest
	gdb -tui ./bin/libquadtest

clean:
	rm -rf bin/*
