#ifndef PSS_BICGSTAB_SOLVER_H_
#define PSS_BICGSTAB_SOLVER_H_

// head file
#include "../include_file/pss_solver_iteration.h"

// function prototype
/*
 * computing residual vector
 * res = b - A * u
 * input: matrix, solution, rhs_vector
 * output: residual_vector
 * return: EXIT_SUCCESS
 * */
int bicgstab_computing_residual_vector( const CSR_MAT * , const double,
	const double *, const double *, double * );

/*
 * computing residual vector part 1
 * temp_vec_1 = alpha I * u( 1 : n ) - W * u( n + 1 : 2n )
 * input: matrix, solution
 * output: residual vector
 * return: EXIT_SUCCESS
 * */
int bicgstab_computing_residual_vector_part_1( const CSR_MAT *, const double,
	const double *, double * );

/*
 * computing residual vector part 2
 * temp_vec_2 = W * u( 1 : n ) + ( alpha I + T ) * u( n + 1 : 2n )
 * input: matrix, solution
 * output: residual vector
 * return: EXIT_SUCCESS
 * */
int bicgstab_computing_residual_vector_part_2( const CSR_MAT *, const double,
	const double *, double * );

/*
 * computing vector inner product
 * input: vec_1, vec_2
 * output: value = ( vec_1, vec_2 )
 * */
double bicgstab_vector_inner_product( const double *, const double *, const int );

/*
 * computing A p_j
 * input: matrix, gradient vector
 * output: temp = A p_j
 * return: EXIT_SUCCESS
 * */
int bicgstab_computing_matrix_gradient_vector_product( const CSR_MAT *, const double,
	const double *, double * );

/*
 * computing s_j = r_j - alpha_j A p_j
 * input: r_j, alpha_j, A p_j
 * output: s_j
 * return: EXIT_SUCCESS
 * */
int bicgstab_computing_vector_s( double *, const double *, const double *,
	const double, const int );

/*
 * updating solution x_{ j + 1 } = x_j + alpha_j p_j + omega_j s_j
 * */
int bicgstab_updating_solution( double *, const double *, const double *,
	const double, const double, const int );

/*
 * updating residual vector
 * r_{ j + 1 } = s_j - omega_j A s_j
 * */
int bicgstab_updating_residual_vector( double *,
	const double *,
	const double *,
	const double,
	const int );

/*
 * updating gradient vector
 * p_{ j + 1 } = r_{ j + 1 } + beta_j ( p_j - omega_j A p_j )
 * */
int bicgstab_updating_gradient_vector( double *,
	const double *,
	const double *,
	const double,
	const double,
	const int );

#endif
