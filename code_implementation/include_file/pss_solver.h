#ifndef PSS_SOLVER_H_
#define PSS_SOLVER_H_

// head file
#include "../include_file/linear_system.h"

// function prototypr
/*
 * PSS iteration method
 * intput: matrix, vector, solution, max_iter_pss, relative tolerance
 * */
int pss_solver_iteration( const CSR_MAT *, const RHS_VEC *, const double, 
	double *, const int, const double );

#endif
