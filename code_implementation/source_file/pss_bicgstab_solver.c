#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../include_file/pss_bicgstab_solver.h"

int pss_bicgstab_solver( const CSR_MAT * bicgstab_mat, const double alpha_pss,
	const double * bicgstab_rhs, double * solution,
	const int bicgstab_max_iter, const double bicgstab_rela_tol )
{
    int index;
    int iter;

    double * temp_solution;    // iterative solution u_k

    if( ( temp_solution = ( double * ) malloc( 2 * ( bicgstab_mat->n_size ) 
		    * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // initialize temp_solution
    for( index = 0; index < 2 * ( bicgstab_mat->n_size ); index++ )
    {
	*( temp_solution + index ) = 0;
    }

    // initialize residual vector
    /*
     * residual_vector = rhs - mat * solution
     * */
    double * bicgstab_residual_vector;
    double bicgstab_residual_vector_norm_2_init;
    double bicgstab_residual_vector_norm_2;
    double bicgstab_relative_residual_norm_2;

    if( ( bicgstab_residual_vector = ( double * ) malloc( 2 * ( bicgstab_mat->n_size ) 
		    * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // computing residual vector BiCGSTAB
    bicgstab_computing_residual_vector( bicgstab_mat, alpha_pss,
	    bicgstab_rhs, temp_solution, bicgstab_residual_vector );

    // computing norm_2 of residual vector
    bicgstab_residual_vector_norm_2_init = computing_vector_norm_2( bicgstab_residual_vector, 
	    2 * ( bicgstab_mat->n_size ) );
    bicgstab_residual_vector_norm_2 = computing_vector_norm_2( bicgstab_residual_vector,
	    2 * ( bicgstab_mat->n_size ) );
    bicgstab_relative_residual_norm_2 = bicgstab_residual_vector_norm_2 /
	bicgstab_residual_vector_norm_2_init;

    // initialize
    iter = 0;
    printf( "\niter\t||r_0||\t||r_k||\t||r_k||/||r_0||\n" );
    printf( "%d\t%.12le\t%.12le\t%.12le\n", iter,
	    bicgstab_residual_vector_norm_2_init,
	    bicgstab_residual_vector_norm_2,
	    bicgstab_relative_residual_norm_2 );

#if 1
    double * bicgstab_gradient_vector;    // p_j
    double * bicgstab_gradient_vector_temp;    // A p_j
    double * bicgstab_residual_vector_arbitrary;    // arbitrary r0*
    
    if( ( bicgstab_gradient_vector = ( double * ) malloc( 2 * ( bicgstab_mat->n_size )
		    * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // initialize gradient vector p_j
    for( index = 0; index < 2 * ( bicgstab_mat->n_size ); index++ )
    {
	*( bicgstab_gradient_vector + index ) = 
	    *( bicgstab_residual_vector + index );
    }

#if 1
    if( ( bicgstab_residual_vector_arbitrary = ( double * ) malloc( 2 * ( bicgstab_mat->n_size )
		    * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

#if 1
    // initialize arbitrary ro*
    for( index = 0; index < 2 * ( bicgstab_mat->n_size ); index++ )
    {
	*( bicgstab_residual_vector_arbitrary + index ) = 
	    *( bicgstab_residual_vector + index );
    }
#endif
#endif

    // memory allocation for A p_j
    if( ( bicgstab_gradient_vector_temp = ( double * ) malloc( 2 * ( bicgstab_mat->n_size )
		    * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }
#endif

    double bicgstab_alpha;
    double bicgstab_alpha_1, bicgstab_alpha_2;

    double * bicgstab_residual_vector_temp;    // s_j
    
    double bicgstab_omega;
    double bicgstab_omega_1, bicgstab_omega_2;

    double * bicgstab_residual_vector_temp_product;    // A s_j

    double bicgstab_beta;
    double bicgstab_beta_1;

    if( ( bicgstab_residual_vector_temp = ( double * ) malloc( 2 * ( bicgstab_mat->n_size )
		    * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    if( ( bicgstab_residual_vector_temp_product = ( double * ) malloc( 2 * ( bicgstab_mat->n_size )
		    * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }
    
    while( iter < bicgstab_max_iter &&
	    bicgstab_relative_residual_norm_2 > bicgstab_rela_tol )
    {
	iter++;
#if 1
	// alpha_1 = ( r_j, r_0* )
	bicgstab_alpha_1 = bicgstab_vector_inner_product( bicgstab_residual_vector,
		bicgstab_residual_vector_arbitrary, 2 * ( bicgstab_mat->n_size ) );

	printf( "\nbicgstab_alpha_1: %.12le\n", bicgstab_alpha_1 );

	// computing A p_j
	bicgstab_computing_matrix_gradient_vector_product( bicgstab_mat, alpha_pss,
		bicgstab_gradient_vector, bicgstab_gradient_vector_temp );

	// alpha_2 = ( A p_j, r_0* )
	bicgstab_alpha_2 = bicgstab_vector_inner_product( bicgstab_gradient_vector_temp,
		bicgstab_residual_vector_arbitrary, 2 * ( bicgstab_mat->n_size ) );

	printf( "\nbicgstab_alpha_2: %.12le\n", bicgstab_alpha_2 );

	// alpha = alpha_1 / alpha_2
	bicgstab_alpha = bicgstab_alpha_1 / bicgstab_alpha_2;

	printf( "\nbicgstab_alpha: %.12le\n", bicgstab_alpha );

	// s_j = r_j - alpha A p_j
	bicgstab_computing_vector_s( bicgstab_residual_vector_temp, bicgstab_residual_vector, 
		bicgstab_gradient_vector_temp, bicgstab_alpha, 2 * ( bicgstab_mat->n_size ) );

	// computing A s_j
	bicgstab_computing_matrix_gradient_vector_product( bicgstab_mat, alpha_pss,
		bicgstab_residual_vector_temp, bicgstab_residual_vector_temp_product );

	// omega_1 = ( A s_j, s_j )
	bicgstab_omega_1 = bicgstab_vector_inner_product( bicgstab_residual_vector_temp_product,
		bicgstab_residual_vector_temp, 2 * ( bicgstab_mat->n_size ) );

	printf( "\nbicgstab_omega_1: %.12le\n", bicgstab_omega_1 );

	// omega_2 = ( A s_j, A s_j )
	bicgstab_omega_2 = bicgstab_vector_inner_product( bicgstab_residual_vector_temp_product,
		bicgstab_residual_vector_temp_product, 2 * ( bicgstab_mat->n_size ) );

	// omega = omega_1 / omega_2
	bicgstab_omega = bicgstab_omega_1 / bicgstab_omega_2;

	printf( "\nbicgstab_omega_2 : %.12le\n", bicgstab_omega_2 );
	printf( "\nbicgstab_omega: %.12le\n", bicgstab_omega );

	// updating solution
	/*
	 * x_{ j + 1 } = x_j + alpha_j p_j + omega_j s_j
	 * */
	bicgstab_updating_solution( temp_solution,
		bicgstab_gradient_vector,
		bicgstab_residual_vector_temp,
		bicgstab_alpha, bicgstab_omega,
		2 * ( bicgstab_mat->n_size ) );

	// updating residual
	/*
	 * r_{ j + 1 } =  s_j - omega_j A s_j
	 * */
	bicgstab_updating_residual_vector( bicgstab_residual_vector,
		bicgstab_residual_vector_temp,
		bicgstab_residual_vector_temp_product,
		bicgstab_omega,
		2 * ( bicgstab_mat->n_size ) );

	// bicgstab_beta_1 = ( r_{ j + 1 }, r0* )
	bicgstab_beta_1 = bicgstab_vector_inner_product( bicgstab_residual_vector,
		bicgstab_residual_vector_arbitrary,
		2 * ( bicgstab_mat->n_size ) );

	// bicgstab_beta = bicgstab_beta_1 / bicgstab_alpha_1 * ( alpha_j / omega_j )
	bicgstab_beta = bicgstab_beta_1 * bicgstab_alpha / bicgstab_alpha_1 / bicgstab_omega;

	printf( "\nbicgstab_beta_1 : %.12le\n", bicgstab_beta_1 );
	printf( "\nbicgstab_beta : %.12le\n", bicgstab_beta );
#if 0
	bicgstab_beta = ( bicgstab_beta_1 / bicgstab_alpha_1 ) 
	    * ( bicgstab_alpha / bicgstab_omega );
#endif
	
	// updating gradient vector
	bicgstab_updating_gradient_vector( bicgstab_gradient_vector,
		bicgstab_residual_vector,
		bicgstab_gradient_vector_temp,
		bicgstab_beta,
		bicgstab_omega,
		2 * ( bicgstab_mat->n_size ) );

	// computing ||r||
	bicgstab_residual_vector_norm_2 = computing_vector_norm_2( bicgstab_residual_vector,
		2 * ( bicgstab_mat->n_size ) );

	printf( "\nbicgstab_residual_vector_norm_2 : %.6f\n", bicgstab_residual_vector_norm_2 );

	bicgstab_relative_residual_norm_2 = bicgstab_residual_vector_norm_2 /
	    bicgstab_residual_vector_norm_2_init;

	printf( "%d\t%.12le\t%.12le\t%.12le\n", iter,
		bicgstab_residual_vector_norm_2_init,
		bicgstab_residual_vector_norm_2,
		bicgstab_relative_residual_norm_2 );
#endif
    }

    // updating solution
    for( index = 0; index < 2 * ( bicgstab_mat->n_size ); index++ )
    {
	*( solution + index ) = *( temp_solution + index );
    }

    // free memory
    free( temp_solution );
    free( bicgstab_residual_vector );
    free( bicgstab_residual_vector_arbitrary );
    free( bicgstab_residual_vector_temp );
    free( bicgstab_residual_vector_temp_product );
    free( bicgstab_gradient_vector );
    free( bicgstab_gradient_vector_temp );

    return 0;
}

int bicgstab_updating_gradient_vector( double * gradient_vector,
	const double * residual_vector,
	const double * gradient_vector_temp,
	const double beta,
	const double omega,
	const int size )
{
    int index;

    for( index = 0; index < size; index++ )
    {
	*( gradient_vector + index ) = ( *( residual_vector + index ) ) 
	    + beta * ( ( *( gradient_vector + index ) ) 
		    - omega * ( *( gradient_vector_temp + index ) ) );
    }

    return 0;
}

int bicgstab_updating_residual_vector( double * residual_vector,
	const double * s_vector,
	const double * s_vector_temp,
	const double omega,
	const int size )
{
    int index;

    for( index = 0; index < size; index++ )
    {
	*( residual_vector + index ) = ( *( s_vector + index ) )
	    - omega * ( *( s_vector_temp + index ) );
    }

    return 0;
}

int bicgstab_updating_solution( double * solution,
	const double * gradient_vector,
	const double * s_vector,
	const double alpha,
	const double omega,
	const int size )
{
    int index;

    for( index = 0; index < size; index++ )
    {
	*( solution + index ) += alpha * ( *( gradient_vector + index ) )
	    + omega * ( *( s_vector + index ) );
    }

    return 0;
}

int bicgstab_computing_vector_s( double * temp_vec,
	const double * residual_vector, 
	const double * gradient_vector,
	const double alpha,
	const int size )
{
    int index;

    for( index = 0; index < size; index++ )
    {
	*( temp_vec + index ) = ( *( residual_vector + index ) ) 
	    - alpha * ( *( gradient_vector + index ) );
    }

    return 0;
}

int bicgstab_computing_matrix_gradient_vector_product( const CSR_MAT * mat, const double alpha,
	const double * gradient_vector, double * gradient_vector_temp )
{
    int index;

    double * temp_vec_1, * temp_vec_2;

    if( ( temp_vec_1 = ( double * ) malloc( ( mat->n_size ) * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    if( ( temp_vec_2 = ( double * ) malloc( ( mat->n_size ) * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // computing temp_vec_1
    /*
     * temp_vec_1 = alpha I * p( 1 : n ) - W * p( n + 1 : 2n )
     * */
    bicgstab_computing_residual_vector_part_1( mat, alpha, gradient_vector, temp_vec_1 );

    // computing temp_vec_2
    /*
     * temp_vec_2 = W * p( 1 : n ) + ( alpha I + T ) * p( n + 1 : 2n )
     * */
    bicgstab_computing_residual_vector_part_2( mat, alpha, gradient_vector, temp_vec_2 );

    for( index = 0; index < mat->n_size; index++ )
    {
	*( gradient_vector_temp + index ) = *( temp_vec_1 + index );
	*( gradient_vector_temp + mat->n_size + index ) = *( temp_vec_2 + index );
    }

    // free memory
    free( temp_vec_1 );
    free( temp_vec_2 );

    return 0;
}

double bicgstab_vector_inner_product( const double * vec_1, const double * vec_2, const int size )
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

int bicgstab_computing_residual_vector( const CSR_MAT * mat, const double alpha,
	const double * rhs, const double * solution, double * residual_vector )
{
    int index;

    // block linear system
    /*
     * [ alpha I          -W ] u( 1 : n )      = rhs( 1 : n )
     * [ W       alpha I + T ] u( n + 1 : 2n ) = rhs( n + 1 : 2n )
     * */

    double * temp_vec_1, * temp_vec_2;

    if( ( temp_vec_1 = ( double * ) malloc( ( mat->n_size ) * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }
    
    if( ( temp_vec_2 = ( double * ) malloc( ( mat->n_size ) * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // temp_vec_1 = alpha I * u( 1 : n ) - W * u( n + 1 : 2n )
    bicgstab_computing_residual_vector_part_1( mat, alpha, solution, temp_vec_1 );

    // temp_vec_2 = W * u( 1 : n ) + ( alpha I + T ) * u( n + 1 : 2n )
    bicgstab_computing_residual_vector_part_2( mat, alpha, solution, temp_vec_2 );

    // residual_vec = [ temp_vec_1, temp_vec_2 ] ^ T
    for( index = 0; index < mat->n_size; index++ )
    {
	*( residual_vector + index ) = ( *( rhs + index ) ) - ( *( temp_vec_1 + index ) );
	*( residual_vector + mat->n_size + index ) = ( *( rhs + mat->n_size + index ) ) 
	    - ( *( temp_vec_2 + index ) );
    }

    // free memory
    free( temp_vec_1 );
    free( temp_vec_2 );

    return 0;
}

int bicgstab_computing_residual_vector_part_1( const CSR_MAT * mat, const double alpha,
	const double * solution, double * temp_vector )
{
    int index, index_j;
    int index_start, index_end;
    double temp_value;

    // alpha I * u( 1 : n )
    for( index = 0; index < mat->n_size; index++ )
    {
	*( temp_vector + index ) = alpha * ( *( solution + index ) );
    }

    // W * u( n + 1 : 2n )
    for( index = 0; index < mat->n_size; index++ )
    {
	temp_value = 0;
	index_start = *( mat->row_pointers + index );
	index_end = *( mat->row_pointers + index + 1 );
	for( index_j = index_start; index_j < index_end; index_j++ )
	{
	    temp_value += ( ( mat->nz_values + index_j )->re_part ) 
		* ( *( solution + mat->n_size + ( *( mat->column_pointers + index_j ) ) ) ) ;
	}
	*( temp_vector + index ) -= temp_value;
    }

    return 0;
}

int bicgstab_computing_residual_vector_part_2( const CSR_MAT * mat, const double alpha,
	const double * solution, double * temp_vector )
{
    int index, index_j;
    int index_start, index_end;
    double temp_value;

    // alpha I * u( n + 1 : 2n )
    for( index = 0; index < mat->n_size; index++ )
    {
	*( temp_vector + index ) = alpha * 
	    ( *( solution + mat->n_size + index ) );
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

    return 0;
}
