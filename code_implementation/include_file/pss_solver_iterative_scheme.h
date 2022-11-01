#ifndef PSS_SOLVER_ITERATIVE_SCHEME_H_
#define PSS_SOLVER_ITERATIVE_SCHEME_H_

// head file
#include "../include_file/linear_system.h"

// global variants
#define CG_MAX_ITER 2000
#define CG_RELA_TOL 1e-8
#define CG_RES_TOL 1e-8
#define GMRES_MAX_ITER 2000
#define GMRES_RELA_TOL 1e-8
#define GMRES_RES_TOL 1e-8
#define BiCGSTAB_MAX_ITER 2000
#define BiCGSTAB_RELA_TOL 1e-8

// function prototype
/*
 * computing rhs vector part 1
 * input: mat, alpha, temp_solution
 * output: matrix-by-vector product
 * return: EXIT_SUCCESS
 * */
int cg_rhs_mvp_1( const CSR_MAT *, const double, const double *, double * );

/*
 * computing rhs vector part 2
 * input: mat, alpha, temp_solution
 * output: matrix-by-vector product
 * return: EXIT_SUCCESS
 * */
int cg_rhs_mvp_2( const CSR_MAT *, const double, const double *, double * );

/*
 * CG solver
 * input: matrix, rhs_vector, max_iter, rela_tol, res_tol
 * output: solution
 * return: EXIT_SUCCESS
 * */
int pss_cg_solver( const CSR_MAT *, const double,
	const double *, double *,
	const int, const double, const double );

/*
 * computing rhs vector part 1
 * input: matrix, alpha, temp_solution
 * output: matrix-by-vector product
 * return EXIT_SUCCESS
 * */
int gmres_rhs_mvp_1( const CSR_MAT *, const double, const double *, double * );

/*
 * computing rhs vector part 2
 * intput: matrix, alpha, temp_solution
 * output: matrix-by-vector product
 * return: EXIT_SUCCESS
 * */
int gmres_rhs_mvp_2( const CSR_MAT *, const double, const double *, double * );

/*
 * GMRES solver
 * input: matrix, rhs_vector, max_iter, rela_tol, res_tol
 * output: solution
 * return: EXIT_SUCCESS
 * */
int pss_gmres_solver( const CSR_MAT *, const double, 
	const double *, double *,
	const int, const double, const double );

/*
 * computing rhs vector part 1
 * input: matrix, alpha, temp_solution
 * output: matrix-by-vector product
 * return: EXIT_SUCCESS
 * */
int bicgstab_rhs_mvp_1( const CSR_MAT *, const double, const double *, double * );

/*
 * computing rhs vector part 2
 * input: matrix, alpha, temp_solution
 * output: matrix-by-vector product
 * return: EXIT_SUCCESS
 * */
int bicgstab_rhs_mvp_2( const CSR_MAT *, const double, const double *, double * );

/*
 * BiCGSTAB solver
 * input: matrix, rhs vector, max_iter, rela_tol
 * output: solution
 * return: EXIT_SUCCESS
 * */
int pss_bicgstab_solver( const CSR_MAT *, const double,
	const double *, double *,
	const int, const double );

#endif
