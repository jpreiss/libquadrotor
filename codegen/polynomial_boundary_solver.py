#!/usr/bin/env python
"""Polynomial boundary condition code generator.

Generates C code of closed form solution for polynomial coefficients
given the values and derivatives at domain boundaries.
The first boundary is fixed at 0, the second boundary is user input.
The generated code is printed to stdout.
The polynomial degree is read from the first command line argument.
The polynomial degree must be odd, so there is exactly one solution.
The code will expect values and (degree - 1) / 2 derivatives at boundaries.
"""

import sys

import numpy as np
import sympy as sym
from sympy.solvers.solveset import linsolve


def poly_deriv_matrix(degree):
	d = range(1, degree + 2)
	D = sym.diag(*d)
	D.col_del(-1)
	D = D.col_insert(0, sym.zeros(degree + 1, 1))
	return D

def time_vec(t, degree):
	return sym.Matrix([t ** i for i in range(degree + 1)])

def boundary_solve(degree):

	size = degree + 1
	assert(size % 2 == 0)

	# vectors of successive powers of time argument at boundaries.
	t0_vec = time_vec(0, degree)
	t1 = sym.Symbol("T1")
	t1_vec = time_vec(t1, degree)
	n_derivs = size // 2

	# matrix representing derivative operator for polynomials
	# represented as vectors of coefficients in increasing order.
	D = poly_deriv_matrix(degree)
	deriv_mats = [sym.eye(size)]
	for i in range(n_derivs - 1):
		deriv_mats.append(D * deriv_mats[-1])

	# left side of the linear system.
	lhs = sym.Matrix(
		[t0_vec.T * D for D in deriv_mats] +
		[t1_vec.T * D for D in deriv_mats])

	# right side of the linear system.
	rhs = sym.Matrix(
		[sym.Symbol("d0_" + str(i)) for i in range(n_derivs)] +
		[sym.Symbol("d1_" + str(i)) for i in range(n_derivs)])

	# solve the linear system symbolically using Gaussian elimination.
	system = (lhs, rhs)
	solution = list(linsolve(system))[0]
	solution = sym.simplify(solution)

	# in the generated code, we want to replace powers of T
	# t1^1, t1^2, ... with named variables T1, T2, ...
	# so we can compute them with successive multiplication in C.
	# the substitution must happen in reverse,
	# otherwise subs() will factor t1**k into t1*t1**(k-1), etc.
	for i in range(degree, 1, -1):
		tpow = t1 ** i
		symver = sym.Symbol("T" + str(i))
		solution = solution.subs(tpow, symver)

	# generate the C code.

	# declare the function.
	head = f"""void solve_poly{degree}_boundary(
	float const vals_init[{n_derivs}],
	float const vals_end[{n_derivs}],
	float const time_end,
	float coefs[{size}])
{{
	float const T1 = time_end;"""

	# compute the successive powers of end time value.
	time_powers = "\n".join([
		f"""	float const T{i} = T{i-1} * T1;"""
		for i in range(2, size)])

	# unpack the boundary values/derivs into scalar variables.
	derivs_begin = "\n".join([
		f"""	float const d0_{i} = vals_init[{i}];"""
		for i in range(n_derivs)])

	derivs_end = "\n".join([
		f"""	float const d1_{i} = vals_end[{i}];"""
		for i in range(n_derivs)])

	# compute the polynomial coefficients using our exact solution.
	coefs = "\n".join([
		f"""	coefs[{i}] = {soln};"""
		for i, soln in enumerate(solution)])

	tail = "}"

	# collect all the steps into one string.
	return "\n\n".join([
		head,
		time_powers,
		derivs_begin,
		derivs_end,
		coefs,
		tail,
	])


def main():
	degree = int(sys.argv[1])
	soln = boundary_solve(degree)
	sys.stdout.write(soln)


if __name__ == "__main__":
	main()

