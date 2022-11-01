#ifndef LINEAR_SYSTEM_SOLVER_H_
#define LINEAR_SYSTEM_SOLVER_H_

// head file
/*
 * linear_system.h
 * */
#include "../include_file/linear_system.h"

// function prototype
/*
 * dividing line
 * -------dividing line------
 * */
void dividing_line_1( const char * );

/*
 * dividing line
 * ======linear system solver======
 * */
void dividing_line_2( const char * );

/*
 * dividing line
 * ======PSS solver======
 * */
void dividing_line_3( const char * );

/*
 * PSS solver
 * input: matrix, rhs vector, max_iteration, relative tolerance
 * return 0 means success
 * computing solution
 * */
int pss_solver( const CSR_MAT *, const RHS_VEC *, RHS_VEC *, const int, const double );

#endif
