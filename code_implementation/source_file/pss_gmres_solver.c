#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../include_file/pss_gmres_solver.h"

int pss_gmres_solver( const CSR_MAT * gmres_mat, const double alpha_pss,
	const double * gmres_rhs, double * solution,
	const int gmres_max_iter, const double gmres_rela_tol, const double gmres_res_tol )
{
    int index;
    int iter_gmres;

    double * temp_solution;    // iterative solution u_k
    double * temp_residual_vector;    // b - A u_k
    double gmres_rhs_norm, temp_residual_vector_norm;

    gmres_rhs_norm = computing_vector_norm_2( gmres_rhs, 2 * gmres_mat->n_size );

    // memory allocation for temp_solution
    if( ( temp_solution = ( double * ) malloc( 2 * ( gmres_mat->n_size ) * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // initialize temp_solution
    for( index = 0; index < gmres_mat->n_size; index++ )
    {
	*( temp_solution + index ) = 0;
	*( temp_solution + gmres_mat->n_size + index ) = 0;
    }

    // memory allocation for residual vector
    if( ( temp_residual_vector = ( double * ) malloc( 2 * ( gmres_mat->n_size ) * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // computing residual vector
    /*
     * gmres_mat * temp_solution = gmres_rhs
     * temp_residual_vector = gmres_rhs - gmres_mat * temp_solution
     * */
    gmres_computing_residual_vector( temp_residual_vector,
	    gmres_mat, alpha_pss, gmres_rhs, temp_solution );

    temp_residual_vector_norm = computing_vector_norm_2( temp_residual_vector, 
	    2 * ( gmres_mat->n_size ) );

    iter_gmres = 0;
    printf( "\niter_gmres\t||b||\t||r||\t||r||/||b||\n" );
    printf( "%d\t%.12le\t%.12le\t%.12le\n", iter_gmres,
	    gmres_rhs_norm,
	    temp_residual_vector_norm,
	    temp_residual_vector_norm / gmres_rhs_norm );

    // GMRES iteration
    while( iter_gmres < ( gmres_max_iter / GMRES_RESTART )
	    && temp_residual_vector_norm > gmres_res_tol )
    {
	iter_gmres++;

	gmres_restart_iteration( temp_solution,
		gmres_mat, alpha_pss, gmres_rhs, GMRES_RESTART );

	gmres_computing_residual_vector( temp_residual_vector,
		gmres_mat, alpha_pss, gmres_rhs, temp_solution );
	temp_residual_vector_norm = computing_vector_norm_2( temp_residual_vector,
		2 * ( gmres_mat->n_size ) );

	printf( "%d\t%.12le\t%.12le\t%.12le\n", iter_gmres,
		gmres_rhs_norm,
		temp_residual_vector_norm,
		temp_residual_vector_norm / gmres_rhs_norm );
    }

    // updating solution
    for( index = 0; index < gmres_mat->n_size; index++ )
    {
	*( solution + index ) = *( temp_solution + index );
	*( solution + gmres_mat->n_size + index ) = *( temp_solution + gmres_mat->n_size + index );
    }

    // free memory
    free( temp_solution );
    free( temp_residual_vector );

    return 0;
}

int gmres_computing_residual_vector( double * residual_vector,
	const CSR_MAT * mat, const double alpha,
	const double * rhs, const double * solution )
{
    int index;

    // residual_vector = rhs - gmres_mat * solution
    /*
     * part 1: residual_vector = rhs
     * */
    for( index = 0; index < mat->n_size; index++ )
    {
	*( residual_vector + index ) = *( rhs + index );
	*( residual_vector + mat->n_size + index ) = *( rhs + mat->n_size + index );
    }

    /*
     * part 2: residual_vector -= gmres_mat * solution
     * computing gmres_mat * solution
     *     [ alpha I         -W ] u( 1 : n )
     *     [ W      alpha I + T ] u( n + 1 : 2n )
     * */
    double * temp_vec_1, * temp_vec_2;

    // memory allocation
    if( ( temp_vec_1 = ( double * ) malloc( mat->n_size * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    if( ( temp_vec_2 = ( double * ) malloc( mat->n_size * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // temp_vec_1 = alpha I * u( 1 : n ) - W * u( n + 1 : 2n )
    gmres_computing_residual_vector_part_1( temp_vec_1,
	    mat, alpha, solution );
    
    // temp_vec_2 = W * u( 1 : n ) + ( alpha I + T ) * u( n + 1 : 2n )
    gmres_computing_residual_vector_part_2( temp_vec_2,
	    mat, alpha, solution );

    /*
     * part 2: residual_vector( 1 : n )      -= temp_vec_1
     *         residual_vector( n + 1 : 2n ) -= temp_vec_2
     * */
    for( index = 0; index < mat->n_size; index++ )
    {
	*( residual_vector + index ) -= *( temp_vec_1 + index );
	*( residual_vector + mat->n_size + index ) -= *( temp_vec_2 + index );
    }

    // free memory
    free( temp_vec_1 );
    free( temp_vec_2 );

    return 0;
}

int gmres_computing_residual_vector_part_1( double * temp_vector,
	const CSR_MAT * mat, const double alpha, const double * solution )
{
    // temp_vec_1 = alpha I * u( 1 : n ) - W * u( n + 1 : 2n )
    int index, index_j;
    int index_start, index_end;
    double temp_value;

    for( index = 0; index < mat->n_size; index++ )
    {
	*( temp_vector + index ) = alpha * ( *( solution + index ) );
    }

    for( index = 0; index < mat->n_size; index++ )
    {
	temp_value = 0;
	index_start = *( mat->row_pointers + index );
	index_end = *( mat->row_pointers + index + 1 );
	for( index_j = index_start; index_j < index_end; index_j++ )
	{
	    temp_value += ( ( mat->nz_values + index_j )->re_part )
		* ( *( solution + mat->n_size + ( *( mat->column_pointers + index_j ) ) ) );
	}
	*( temp_vector + index ) -= temp_value;
    }

    return 0;
}

int gmres_computing_residual_vector_part_2( double * temp_vector,
	const CSR_MAT * mat, const double alpha, const double * solution )
{
    // temp_vec_2 = W * u( 1 : n ) + ( alpha I + T ) * u( n + 1 : 2n )
    int index, index_j;
    int index_start, index_end;
    double temp_value;

    // alpha I * u( n + 1 : 2n )
    for( index = 0; index < mat->n_size; index++ )
    {
	*( temp_vector + index ) = alpha * ( *( solution + mat->n_size + index ) );
    }

    // T * u( n + 1 : 2n )
    for( index = 0; index < mat->n_size; index++ )
    {
	temp_value = 0;
	index_start = *( mat->row_pointers + index );
	index_end = *( mat->row_pointers + index + 1 );
	for( index_j = index_start; index_j < index_end; index_j++ )
	{
	    temp_value += ( ( mat->nz_values + index_j )->im_part ) 
		* ( *( solution + mat->n_size + ( *( mat->column_pointers + index_j ) ) ) );
	}
	*( temp_vector + index ) += temp_value;
    }

    // W * u( 1 : n )
    for( index = 0; index < mat->n_size; index++ )
    {
	temp_value = 0;
	index_start = *( mat->row_pointers + index );
	index_end = *( mat->row_pointers + index + 1 );
	for( index_j = index_start; index_j < index_end; index_j++ )
	{
	    temp_value += ( ( mat->nz_values + index_j )->re_part ) 
		* ( *( solution + ( *( mat->column_pointers + index_j ) ) ) );
	}
	*( temp_vector + index ) += temp_value;
    }

    return 0;
}
