#ifndef PSS_GMRES_SOLVER_H_
#define PSS_GMRES_SOLVER_H_

// include file
#include "../include_file/pss_solver_iteration.h"

// global variants
#define GMRES_RESTART 30

// function prototype
/*
 * computing residual vector
 * input: gmres_mat, gmres_rhs, temp_solution
 * ouput: residual_vector
 * return: EXIT_SUCCESS
 * */
int gmres_computing_residual_vector( double *,
	const CSR_MAT *, const double, const double *, const double *);

/*
 * gmres iteration with restart
 * input: gmres_mat, gmres_rhs, restart_num
 * output: temp_solution
 * return: EXIT_SUCCESS
 * */
int gmres_restart_iteration( double *, 
	const CSR_MAT *, const double, const double *, const int );

/*
 * computing residual vector part 1
 * input: matrix, solution
 * output: temp_vec_1
 * */
int gmres_computing_residual_vector_part_1( double *,
	const CSR_MAT *, const double, const double * );

/*
 * computing residual vector part 2
 * input: matrix, solution
 * output: temp_vec_2
 * */
int gmres_computing_residual_vector_part_2( double *,
	const CSR_MAT *, const double, const double * );

#endif
