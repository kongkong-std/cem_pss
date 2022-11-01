#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "../include_file/pss_solver_iteration.h"

int pss_solver_iteration( const CSR_MAT * mat, const RHS_VEC * b_vec, const double alpha_pss,
	double * solution, const int max_iter, const double rela_tol )
{
    double res_norm_2, res_norm_2_0, rela_res_norm_2;
    double * temp_solution;
    double * residual_vector, * rhs_vector;
    int index, iter;

    // memory allocation for temporary solution
    if( ( temp_solution = ( double * ) malloc( 2 * ( b_vec->n_size ) * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // initialize
    for( index = 0; index < 2 * ( b_vec->n_size ); index++ )
    {
	*( temp_solution + index ) = 0;
    }

    // memory allocation for residual vector
    if( ( residual_vector = ( double * ) malloc( 2 * ( b_vec->n_size ) * sizeof( double ) ) )
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    if( ( rhs_vector = ( double * ) malloc( 2 * ( b_vec->n_size ) * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    for( index = 0; index < b_vec->n_size; index++ )
    {
	*( rhs_vector + index ) = ( b_vec->vec_values + index )->im_part;
	*( rhs_vector + b_vec->n_size + index ) = ( b_vec->vec_values + index )->re_part;
    }

    //=====computing residual vector
    computing_residual_vector( mat, b_vec, temp_solution, residual_vector );
    //=====end computing residual vector
    
    iter = 0;

    //=====computing L2 norm
    res_norm_2_0 = computing_vector_norm_2( rhs_vector, 2 * ( b_vec->n_size ) );
    res_norm_2 = computing_vector_norm_2( residual_vector, 2 * ( b_vec->n_size ) );
    rela_res_norm_2 = res_norm_2 / res_norm_2_0;
    //=====end computing L2 norm
    
    // norm of residual vector
    FILE * fp;
    if( ( fp = fopen( "../result_file/residual_norm.txt", "w" ) ) 
	    == NULL )
    {
	fprintf( stdout, "Cannot write to file!\n" );
	exit( EXIT_FAILURE );
    }

    fprintf( fp, "iter\t||b||\t||r||\t||r||/||b||\n" );
    fprintf( fp, "%d\t%.12le\t%.12le\t%.12le\n", iter, res_norm_2_0, res_norm_2,
	    rela_res_norm_2 );

    // iteration scheme
    while( rela_res_norm_2 > rela_tol &&
	    iter < max_iter )
    {
	iter++;
	printf( "\ninformation iter %d:\n", iter );
	pss_solver_iterative_scheme( mat, b_vec, alpha_pss, temp_solution );

	// computing norm of residual vector
        computing_residual_vector( mat, b_vec, temp_solution, residual_vector );
	res_norm_2 = computing_vector_norm_2( residual_vector, 2 * ( b_vec->n_size ) );
	rela_res_norm_2 = res_norm_2 / res_norm_2_0;

	fprintf( fp, "%d\t%.12le\t%.12le\t%.12le\n", iter, res_norm_2_0,
		res_norm_2, rela_res_norm_2 );
    }
    

    fclose( fp );

#if 0
    // test
    for( index = 0; index < 2 * b_vec->n_size; index++ )
    {
	*( temp_solution + index ) = index / 10.;
    }
#endif

    // updating solution
    for( index = 0; index < 2 * b_vec->n_size; index++ )
    {
	*( solution + index ) = *( temp_solution + index );
    }

    // free memory
    free( temp_solution );
    free( residual_vector );
    free( rhs_vector );

    return 0;
}

int computing_residual_vector( const CSR_MAT * mat, const RHS_VEC * b_vec,
	const double * temp, double * residual )
{
    matrix_by_vector_product( mat, temp, residual );

    int index;
    for( index = 0; index < mat->n_size; index++ )
    {
	*( residual + index ) = ( b_vec->vec_values + index )->im_part
	    - ( *( residual + index ) );
	*( residual + mat->n_size + index ) = ( b_vec->vec_values + index )->re_part 
	    - ( *( residual + mat->n_size + index ) );
    }

    return 0;
}

int matrix_by_vector_product( const CSR_MAT * mat, const double * temp, double * residual )
{
    int index;
    double * temp_1, * temp_2, * temp_3, * temp_4;

    // memory allocation
    if( ( temp_1 = ( double * ) malloc( mat->n_size * sizeof( double ) ) ) == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    if( ( temp_2 = ( double * ) malloc( mat->n_size * sizeof( double ) ) ) == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    if( ( temp_3 = ( double * ) malloc( mat->n_size * sizeof( double ) ) ) == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    if( ( temp_4 = ( double * ) malloc( mat->n_size * sizeof( double ) ) ) == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // temp_1 = T u_re
    block_matrix_by_vector_product_complex( mat, temp, temp_1 );

    // temp_2 = W u_im
    block_matrix_by_vector_product_real( mat, temp + mat->n_size, temp_2 );
    
    // temp_3 = W u_re
    block_matrix_by_vector_product_real( mat, temp, temp_3 );

    // temp_4 = T u_im
    block_matrix_by_vector_product_complex( mat, temp + mat->n_size, temp_4 );

    for( index = 0; index < mat->n_size; index++ )
    {
	// T u_re - W u_im = temp_1 - temp_2
	*( residual + index ) = ( *( temp_1 + index ) ) - ( *( temp_2 + index ) );

	// W u_re + T u_im = temp_3 + temp_4
	*( residual + mat->n_size + index ) = ( *( temp_3 + index ) ) + ( *( temp_4 + index ) );
    }

    // free memory
    free( temp_1 );
    free( temp_2 );
    free( temp_3 );
    free( temp_4 );

    return 0;
}

int block_matrix_by_vector_product_real( const CSR_MAT * mat, const double * temp, double * result )
{
    int index, index_j, index_start, index_end;
    double sum;

    // computing spMVP
    for( index = 0; index < mat->n_size; index++ )
    {
	sum = 0;
	index_start = *( mat->row_pointers + index );
	index_end = *( mat->row_pointers + index + 1 );
	for( index_j = index_start; index_j < index_end; index_j++ )
	{
	    sum += ( ( mat->nz_values + index_j )->re_part ) 
		* ( *( temp + ( *( mat->column_pointers + index_j ) ) ) );
	}
	*( result + index ) = sum;
    }

    return 0;
}

int block_matrix_by_vector_product_complex( const CSR_MAT * mat, const double * temp, double * result )
{
    int index, index_j, index_start, index_end;
    double sum;

    for( index = 0; index < mat->n_size; index++ )
    {
	sum = 0;
	index_start = *( mat->row_pointers + index );
	index_end = *( mat->row_pointers + index + 1 );
	for( index_j = index_start; index_j < index_end; index_j++ )
	{
	    sum += ( (mat->nz_values + index_j )->im_part )
		*( *( temp + ( *( mat->column_pointers + index_j ) ) ) );
	}
	*( result + index ) = sum;
    }

    return 0;
}

double computing_vector_norm_2( const double * vec_pointer, int size )
{
    int index;
    double sum;

    // initialize
    sum = 0;

    for( index = 0; index < size; index++ )
    {
	sum += ( *( vec_pointer + index ) ) * ( *( vec_pointer + index ) );
    }

    return sqrt( sum );
}
