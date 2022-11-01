#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../include_file/pss_cg_solver.h"

/*
 * ( alpha I + T ) * u( 1 : n ) = cg_rhs( 1 : n )
 * alpha I * u( n + 1 : 2n )    = cg_rhs( n + 1 : 2n )
 * */
int pss_cg_solver( const CSR_MAT * cg_mat, const double alpha_pss,
	const double * cg_rhs, double * solution,
	const int cg_max_iter, const double cg_rela_tol, const double cg_res_tol )
{
    int index;
    int iter;

    // alpha I * u( n + 1 : 2n ) = cg_rhs( n + 1 : 2n )
    for( index = 0; index < cg_mat->n_size; index++ )
    {
	*( solution + cg_mat->n_size + index ) = 
	    ( *( cg_rhs + cg_mat->n_size + index ) ) / alpha_pss;
    }

    double * temp_solution;    // iterative solution u_k
    double * cg_residual_vector;
    double * cg_gradient_vector, * temp_cg_gradient_vector;
    double cg_residual_vector_norm_2;
    double cg_rhs_vector_norm_2;
    double cg_rela_residual_norm_2;

    // memory allocation for temp_solution
    if( ( temp_solution = ( double * ) malloc( ( cg_mat->n_size ) * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // initialize for temp_solution
    for( index = 0; index < cg_mat->n_size; index++ )
    {
	*( temp_solution + index ) = 0;
    }

    // memory allocation for residual vector
    if( ( cg_residual_vector = ( double * ) malloc( ( cg_mat->n_size ) * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // computing residual vector
    cg_solver_residual_vector( cg_mat, alpha_pss, 
	    cg_rhs, temp_solution, cg_residual_vector );

    // memory allocation for gradient vector
    if( ( cg_gradient_vector = ( double * ) malloc( ( cg_mat->n_size ) * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // initialize for gradient vector
    for( index = 0; index < cg_mat->n_size; index++ )
    {
	*( cg_gradient_vector + index ) = *( cg_residual_vector + index );
    }

    // initialize for temporary gradient vector
    if( ( temp_cg_gradient_vector = ( double * ) malloc( ( cg_mat->n_size ) * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    cg_residual_vector_norm_2 = computing_vector_norm_2( cg_residual_vector, cg_mat->n_size );
    cg_rhs_vector_norm_2 = computing_vector_norm_2( cg_rhs, cg_mat->n_size );
    cg_rela_residual_norm_2 = cg_residual_vector_norm_2 / cg_rhs_vector_norm_2;

    iter = 0;
    printf( "iter\t||b||\t||r||\t||r||/||b||\n" );
    printf( "%d\t%.12le\t%.12le\t%.12le\n", iter,
	    cg_rhs_vector_norm_2,
	    cg_residual_vector_norm_2,
	    cg_rela_residual_norm_2 );

    double cg_alpha, cg_beta;    // alpha = ( r_j, r_j ) / ( A p_j, p_j )
    double cg_alpha_1, cg_alpha_2, cg_alpha_3;    // alpha = alpha_1 / alpha_2
    while( iter < cg_max_iter
	    && cg_residual_vector_norm_2 > cg_res_tol )
    {
	iter++;

	// cg_alpha_1 = ( r_j, r_j )
	cg_alpha_1 = cg_vector_inner_product( cg_residual_vector, cg_residual_vector, cg_mat->n_size );
#if 0
	printf( "\ncg_alpha_1 : %.6f\n", cg_alpha_1 );
#endif

	// cg_alpha_2 = ( A p_j, p_j )
	// temp = A p_j = ( alpha I + T ) p_j
	cg_matrix_by_vector_product( cg_mat, alpha_pss, cg_gradient_vector, temp_cg_gradient_vector );
	cg_alpha_2 = cg_vector_inner_product( temp_cg_gradient_vector, cg_gradient_vector, cg_mat->n_size );
#if 0
	printf( "\ncg_alpha_2 : %.6f\n", cg_alpha_2 );
#endif

	cg_alpha = cg_alpha_1 / cg_alpha_2;

	// updating solution
	/*
	 * x_{j + 1} = x_j + cg_alpha * gradient_vector
	 * */
	cg_updating_solution( temp_solution, cg_gradient_vector, cg_alpha, cg_mat->n_size );

	// updating residual vector
	/*
	 * residual_vector = rhs_vec( 1 : n ) - ( alpha I + T ) u( 1 : n )
	 *                 = r - alpha_j ( alpha I + T ) A p_j
	 * */
#if 0
	cg_solver_residual_vector( cg_mat, alpha_pss,
		cg_rhs, solution, cg_residual_vector );
#endif

	cg_updating_residual_vector( cg_residual_vector, cg_alpha, 
		temp_cg_gradient_vector, cg_mat->n_size );

	// conputing ( r_{ j + 1 }, r_{ j + 1 } )
	cg_alpha_3 = cg_vector_inner_product( cg_residual_vector, cg_residual_vector, cg_mat->n_size );

	// computing beta = ( r_{ j + 1 }, r_{ j + 1 } ) / ( r_j, r_j )
	cg_beta = cg_alpha_3 / cg_alpha_1;

	// updating gradient vector
	cg_updating_gradient_vector( cg_gradient_vector, cg_residual_vector,
		cg_beta, cg_mat->n_size);

	// information of iteration
	cg_residual_vector_norm_2 = computing_vector_norm_2( cg_residual_vector, cg_mat->n_size );
	cg_rela_residual_norm_2 = cg_residual_vector_norm_2 / cg_rhs_vector_norm_2;
	printf( "%d\t%.12le\t%.12le\t%.12le\n", iter,
		cg_rhs_vector_norm_2,
		cg_residual_vector_norm_2,
		cg_rela_residual_norm_2 );
    }

    // updating solution
    for( index = 0; index < cg_mat->n_size; index++ )
    {
	*( solution + index ) = *( temp_solution + index );
    }

    // free memory
    free( temp_solution );
    free( cg_residual_vector );
    free( cg_gradient_vector );
    free( temp_cg_gradient_vector );

    return 0;
}

int cg_updating_residual_vector( double * residual, const double alpha,
	const double * temp_gradient, const int size )
{
    int index;

    for( index = 0; index < size; index++ )
    {
	*( residual + index ) = *( residual + index ) 
	    - alpha * ( *( temp_gradient + index ) );
    }

    return 0;
}

int cg_updating_gradient_vector( double * gradient_vector, const double * residual_vector,
	const double beta, int size )
{
    int index;

    for( index = 0; index < size; index++ )
    {
	*( gradient_vector + index ) = ( *( residual_vector + index ) ) 
	    + ( beta * ( *( gradient_vector + index ) ) );
    }

    return 0;
}

int cg_updating_solution( double * solution, const double * gradient_solution,
	const double alpha, const int size )
{
    int index;

    for( index = 0; index < size; index++ )
    {
	*( solution + index ) += alpha * ( *( gradient_solution + index ) );
    }

    return 0;
}

int cg_matrix_by_vector_product( const CSR_MAT * cg_mat, const double alpha,
	const double * gradient_vector, double * temp_gradient_vector )
{
    int index;

    // computing alpha I p_j
    for( index = 0; index < cg_mat->n_size; index++ )
    {
	*( temp_gradient_vector + index ) = alpha * 
	    ( *( gradient_vector + index ) );
    }

    int index_start, index_end, index_j;
    double temp_value;

    // computing T p_j
    for( index = 0; index < cg_mat->n_size; index++ )
    {
	temp_value = 0;
	index_start = *( cg_mat->row_pointers + index );
	index_end = *( cg_mat->row_pointers + index + 1 );
	for( index_j = index_start; index_j < index_end; index_j++ )
	{
	    temp_value += ( ( cg_mat->nz_values + index_j )->im_part ) 
		* ( *( gradient_vector + ( *( cg_mat->column_pointers + index_j ) ) ) );
	}
	*( temp_gradient_vector + index ) += temp_value;
    }

    return 0;
}

double cg_vector_inner_product( const double * vec_1, const double * vec_2, const int size )
{
    int index;
    double value;

    // initialize
    value = 0;

    for( index = 0; index < size; index++ )
    {
	value += ( *( vec_1 + index ) ) * ( *( vec_2 + index ) );
    }

    return value;
}

int cg_solver_residual_vector( const CSR_MAT * cg_mat, const double alpha_pss,
	const double * cg_rhs, const double * solution, double * cg_residual_vector )
{
    int index, index_j;
    int index_start, index_end;
    double temp_value;

    // cg_residual_vector = cg_rhs( 1 : n ) - alpha I * u( 1 : n )
    for( index = 0; index < cg_mat->n_size; index++ )
    {
	*( cg_residual_vector + index ) = *( cg_rhs + index )
	    - alpha_pss * ( *( solution + index ) ) ;
    }

    // cg_residual_vector = cg_residual_vector - T * u( 1 : n )
    for( index = 0; index < cg_mat->n_size; index++ )
    {
	temp_value = 0;
	index_start = *( cg_mat->row_pointers + index );
	index_end = *( cg_mat->row_pointers + index + 1 );
	for( index_j = index_start; index_j < index_end; index_j++ )
	{
	    temp_value += ( ( cg_mat->nz_values + index )->im_part ) 
		* ( *( solution + ( *( cg_mat->column_pointers + index_j ) ) ) );
	}
	*( cg_residual_vector + index ) -= temp_value;
    }

    return 0;
}
