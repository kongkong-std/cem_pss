#ifndef PSS_CG_SOLVER_H_
#define PSS_CS_SOLVER_H_

// head file
//#include "../include_file/linear_system.h"
#include "../include_file/pss_solver_iteration.h"

// function prototype
/*
 * ( alpha I + T ) * u( 1 : n ) = cg_rhs( 1 : n )
 * residual_vector = cg_rhs( 1 : n ) - ( alpha I + T ) * u( 1 : n )
 * */
int cg_solver_residual_vector( const CSR_MAT *, const double,
	const double *, const double *, double * );

/*
 * vector inner product
 * input: vec_1, vec_2, size
 * return: value = ( vec_1, vec_2 )
 * */
double cg_vector_inner_product( const double *, const double *, const int );

/*
 * computing matrix-by-vector product
 * A p_j = ( alpha I + T ) p_j
 * input: mat, alpha, p_j
 * output: temp_p_j
 * return EXIT_SUCCESS
 * */
int cg_matrix_by_vector_product( const CSR_MAT *, const double,
	const double *, double * );

/*
 * updating solution
 * x_{ j + 1 } = x_j + alpha_j p_j
 * */
int cg_updating_solution( double *, const double *, const double, const int );

/*
 * updating gradient vector
 * p_{ j + 1 } = r_{ j + 1 } + beta_j p_j
 * */
int cg_updating_gradient_vector( double *, const double *, const double, const int );

/*
 * updating residual vector
 * r_{ j + 1 } = r_j - alpha_j A p_j
 * */
int cg_updating_residual_vector( double *, const double, const double *, const int );

#endif
