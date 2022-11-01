#ifndef PSS_SOLVER_ITERATION_H_
#define PSS_SOLVER_ITERATION_H_

// head file
#include "../include_file/linear_system.h"

// function prototype
/*
 * computing L2 norm of vector
 * */
double computing_vector_norm_2( const double *, int );

/*
 * computing residual vector
 * input: matrix, rhs vector, temporary vector, residual vector
 * return: residual vector
 * */
int computing_residual_vector( const CSR_MAT *, const RHS_VEC *, const double *, double * );

/*
 * computing matrix-by-vector product
 * input: matrix A, vector b
 * output: computing A x b
 * */
int matrix_by_vector_product( const CSR_MAT *, const double *, double * );

/*
 * block matrix-by-vector product
 * input: block matrix T, vector b
 * output: T x b
 * */
int block_matrix_by_vector_product_complex( const CSR_MAT *, const double *, double * );

/*
 * block matrix-by-vector product
 * input: block matrix W, vector b
 * output: W x b
 * */
int block_matrix_by_vector_product_real( const CSR_MAT *, const double *, double * );

/*
 * PSS iterative scheme
 * input: mat, rhs_vector, alpha_pss
 * ouput: iteration solution
 * return: EXIT_SUCCESS
 * */
int pss_solver_iterative_scheme( const CSR_MAT *, const RHS_VEC *, const double, double * );

#endif
