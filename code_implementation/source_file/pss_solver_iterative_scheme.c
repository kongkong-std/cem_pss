#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../include_file/pss_solver_iterative_scheme.h"

int pss_solver_iterative_scheme( const CSR_MAT * mat, const RHS_VEC * b_vec,
	const double alpha_pss, double * temp_solution )
{
    // block linear system
    /*
     * [ T  -W ] u( 1 : n )      = im( b )
     * [ W   T ] u( n + 1 : 2n ) = re( b )
     * T = im( mat ), W = re( mat )
     * */

    // matrix splitting iteration
    /*
     * part 1:
     * [ alpha I + T  O ] u( 1 : n ) 	  = [ alpha I        W ] u( 1 : n ) + im( b )
     * [ O      alpha I ] u( n + 1 : 2n ) = [ -W   alpha I - T ] u( n + 1 : 2n ) + re( b )
     * */

    int index;
    double * rhs_b_cg;
    double * temp_vec_1, * temp_vec_2;

    // memory allocation for CG solver rhs vector
    if( ( rhs_b_cg = ( double * ) malloc( 2 * ( b_vec->n_size ) * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // initialize for rhs vector
    for( index = 0; index < b_vec->n_size; index++ )
    {
	*( rhs_b_cg + index ) = ( b_vec->vec_values + index )->im_part;
	*( rhs_b_cg + b_vec->n_size + index ) = ( b_vec->vec_values + index )->re_part;
    }

    if( ( temp_vec_1 = ( double * ) malloc( ( b_vec->n_size ) * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    if( ( temp_vec_2 = ( double * ) malloc( ( b_vec->n_size ) * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // computing block rhs vector
    // temp_vec_1 = alpha_pss I u( 1 : n ) + W u( n + 1 : 2n )
    cg_rhs_mvp_1( mat, alpha_pss, temp_solution, temp_vec_1 );

    // temp_vec_2 = -W u( 1 : n ) + ( alpha_pss I - T ) u( n + 1 : 2n )
    cg_rhs_mvp_2( mat, alpha_pss, temp_solution, temp_vec_2 );

    // computing rhs vector 
    for( index = 0; index < b_vec->n_size; index++ )
    {
	*( rhs_b_cg + index ) += *( temp_vec_1 + index );
	*( rhs_b_cg + b_vec->n_size + index ) += *( temp_vec_2 + index );
    }

#if 0
    // print CG rhs_vector
    FILE * fp;

    if( ( fp = fopen( "../result_file/cg_rhs_vector.txt", "w" ) ) 
	    == NULL )
    {
	fprintf( stdout, "Cannot open file!\n" );
	exit( EXIT_FAILURE );
    }

    fprintf( fp, "%d\n", 2 * b_vec->n_size );

    for( index = 0; index < 2 * ( b_vec->n_size ); index++ )
    {
	fprintf( fp, "%.12le\n", *( rhs_b_cg + index ) );
    }

    fclose( fp );
#endif
    
    // CG solver parameter
    puts( "\n==========CG solver==========" );
    
    puts( "parameters for CG:" );
    printf( "CG_MAX_ITER\tCG_RELA_TOL\tCG_RES_TOL\n" );
    printf( "%d\t%.12le\t%.12le\n", CG_MAX_ITER, CG_RELA_TOL, CG_RES_TOL );
    
    pss_cg_solver( mat, alpha_pss, rhs_b_cg, temp_solution, CG_MAX_ITER, CG_RELA_TOL, CG_RES_TOL );
    
    puts( "\n==========end CG solver==========" );

    // matrix splitting iteration
    /*
     * part 2:
     * [ alpha I        -W ] u( 1 : n )      = [ alpha I - T    O ] u( 1 : n ) + im( b )
     * [ W     alpha I + T ] u( n + 1 : 2n ) = [ O        alpha I ] u( n + 1 : 2n ) + re( b )
     * */
    double * rhs_b_gmres;

    // memory allocation for GMRES solver rhs vector
    if( ( rhs_b_gmres = ( double * ) malloc( 2 * ( b_vec->n_size ) * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // initialize rhs_b_gmres
    for( index = 0; index < b_vec->n_size; index++ )
    {
	*( rhs_b_gmres + index ) = ( b_vec->vec_values + index )->im_part;
	*( rhs_b_gmres + b_vec->n_size + index ) = ( b_vec->vec_values + index )->re_part;
    }

    // computing block linear system rhs vector
    /*
     * temp_vec_1 = ( alpha I - T ) u( 1 : n )
     * */
    gmres_rhs_mvp_1( mat, alpha_pss, temp_solution, temp_vec_1 );

    /*
     * temp_vec_2 = alpha I * u( n + 1 : 2n )
     * */
    gmres_rhs_mvp_2( mat, alpha_pss, temp_solution, temp_vec_2 );

    // computing gmres rhs vector
    for( index = 0; index < b_vec->n_size; index++ )
    {
	*( rhs_b_gmres + index ) += *( temp_vec_1 + index );
	*( rhs_b_gmres + b_vec->n_size + index ) += *( temp_vec_2 + index );
    }

    // BiCGSTAB solver parameter
    puts( "\n==========GMRES solver==========" );
    
    puts( "parameters for GMRES:" );
    printf( "GMRES_MAX_ITER\tGMRES_RELA_TOL\tGMRES_RES_TOL\n" );
    printf( "%d\t%.12le\t%.12le\n", GMRES_MAX_ITER, GMRES_RELA_TOL, GMRES_RES_TOL );
    
    pss_gmres_solver( mat, alpha_pss, rhs_b_gmres, temp_solution, 
	    GMRES_MAX_ITER, GMRES_RELA_TOL, GMRES_RES_TOL );
    
    puts( "\n==========end GMRES solver==========" );

#if 0
    double * rhs_b_bicgstab;

    // memory allocation for BiCGSTAB solver rhs vector
    if( ( rhs_b_bicgstab = ( double * ) malloc( 2 * ( b_vec->n_size ) * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // initialize rhs_b_bicgstab
    for( index = 0; index < b_vec->n_size; index++ )
    {
	*( rhs_b_bicgstab + index ) = ( b_vec->vec_values + index )->im_part;
	*( rhs_b_bicgstab + b_vec->n_size + index ) = ( b_vec->vec_values + index )->re_part;
    }

    // computing block linear system rhs vector
    /*
     * temp_vec_1 = ( alpha I - T ) u( 1 : n )
     * */
    bicgstab_rhs_mvp_1( mat, alpha_pss, temp_solution, temp_vec_1 );

    /*
     * temp_vec_2 = alpha I u( n + 1 : 2n )
     * */
    bicgstab_rhs_mvp_2( mat, alpha_pss, temp_solution, temp_vec_2 );

    // computint bicgstab rhs vector
    for( index = 0; index < b_vec->n_size; index++ )
    {
	*( rhs_b_bicgstab + index ) += *( temp_vec_1 + index );
	*( rhs_b_bicgstab + b_vec->n_size + index ) += *( temp_vec_2 + index );
    }

#if 0
    // print BiCGSTAB rhs vector
    FILE * fp;

    if( ( fp = fopen( "../result_file/bicgstab_rhs_vector.txt", "w" ) ) 
	    == NULL )
    {
	fprintf( stdout, "Cannot write file!\n" );
	exit( EXIT_FAILURE );
    }

    for( index = 0; index < b_vec->n_size; index++ )
    {
	fprintf( fp, "%.12le\t%.12le\n", *( rhs_b_bicgstab + index ),
		*( rhs_b_bicgstab + b_vec->n_size + index ) );
    }

    fclose( fp );
#endif
    
    // BiCGSTAB solver parameter
    puts( "\n==========BiCGSTAB solver==========" );
    
    puts( "parameters for BiCGSTAB:" );
    printf( "BiCGSTAB_MAX_ITER\tBiCGSTAB_RELA_TOL\n" );
    printf( "%d\t%.12le\n", BiCGSTAB_MAX_ITER, BiCGSTAB_RELA_TOL );
    
    pss_bicgstab_solver( mat, alpha_pss, rhs_b_cg, temp_solution, 
	    BiCGSTAB_MAX_ITER, BiCGSTAB_RELA_TOL );
    
    puts( "\n==========end BiCGSTAB solver==========" );
#endif

    // free memory
    free( rhs_b_cg );
    free( rhs_b_gmres );
    //free( rhs_b_bicgstab );
    free( temp_vec_1 );
    free( temp_vec_2 );

    return 0;
}

int gmres_rhs_mvp_1( const CSR_MAT * mat, const double alpha, const double * temp_solution,
	double * temp_vector )
{
    int index, index_j;
    int index_start, index_end;
    double temp_value;

    // temp_vector = [ alpha I - T ] u( 1 : n ) = alpha I * u( 1 : n ) - T * u( 1 : n )
    for( index = 0; index < mat->n_size; index++ )
    {
	*( temp_vector + index ) = alpha * ( *( temp_solution + index ) );
    }
    for( index = 0; index < mat->n_size; index++ )
    {
	temp_value = 0;
	index_start = *( mat->row_pointers + index );
	index_end = *( mat->row_pointers + index + 1 );
	for( index_j = index_start; index_j < index_end; index_j++ )
	{
	    temp_value += ( ( mat->nz_values + index_j )->im_part ) 
		* ( *( temp_solution + ( *( mat->column_pointers + index_j ) ) ) );
	}
	*( temp_vector + index ) -= temp_value;
    }

    return 0;
}

int gmres_rhs_mvp_2( const CSR_MAT * mat, const double alpha, const double * temp_solution,
	double * temp_vector )
{
    int index;

    // temp_vector = alpha I * u( n + 1 : 2n )
    for( index = 0; index < mat->n_size; index++ )
    {
	*( temp_vector + index ) = alpha * ( *( temp_solution + mat->n_size + index ) );
    }

    return 0;
}

int bicgstab_rhs_mvp_1( const CSR_MAT * mat, const double alpha, const double * temp_solution,
	double * temp )
{
    int index, index_j;
    int index_start, index_end;
    double temp_value;

    // temp = ( alpha I - T ) u( 1 : n ) = alpha I u( 1 : n ) - T u( 1 : n )
    for( index = 0; index < mat->n_size; index++ )
    {
	*( temp + index ) = alpha * ( *( temp_solution + index ) );
    }
    for( index = 0; index < mat->n_size; index++ )
    {
	temp_value = 0;
	index_start = *( mat->row_pointers + index );
	index_end = *( mat->row_pointers + index + 1 );
	for( index_j = index_start; index_j < index_end; index_j++ )
	{
	    temp_value += ( ( mat->nz_values + index_j )->im_part ) 
		* ( *( temp_solution + ( *( mat->column_pointers + index_j ) ) ) );
	}
	*( temp + index ) -= temp_value;
    }

    return 0;
}

int bicgstab_rhs_mvp_2( const CSR_MAT * mat, const double alpha, const double * temp_solution,
	double * temp )
{
    int index;

    // temp = alpha I u( n + 1 : 2n )
    for( index = 0; index < mat->n_size; index++ )
    {
	*( temp + index ) = alpha * ( *( temp_solution + mat->n_size + index ) );
    }

    return 0;
}

int cg_rhs_mvp_1( const CSR_MAT * mat, const double alpha, const double * temp_solution, 
	double * temp )
{
    int index, index_j;
    int index_start, index_end;
    double temp_value;

    // temp = alpha I * u( 1 : n ) + W * u( n + 1 : 2n )
    for( index = 0; index < mat->n_size; index++ )
    {
	*( temp + index ) = alpha * ( *( temp_solution + index ) );
    }

    for( index = 0; index < mat->n_size; index++ )
    {
	temp_value = 0;
	index_start = *( mat->row_pointers + index );
	index_end = *( mat->row_pointers + index + 1 );
	for( index_j = index_start; index_j < index_end; index_j++ )
	{
	    temp_value += ( ( mat->nz_values + index_j )->re_part ) 
		* ( *( temp_solution + mat->n_size + ( *( mat->column_pointers + index_j ) ) ) );
	}
	*( temp + index ) += temp_value;
    }

    return 0;
}

int cg_rhs_mvp_2( const CSR_MAT * mat, const double alpha, const double * temp_solution, 
	double * temp )
{
    int index, index_j;
    int index_start, index_end;
    double temp_value;

    // -W * u( 1 : n ) + ( alpha I - T ) * u( n + 1 : 2n )
    for( index = 0; index < mat->n_size; index++ )
    {
	*( temp + index ) = alpha * ( *( temp_solution + mat->n_size + index ) );
    }

    for( index = 0; index < mat->n_size; index++ )
    {
	temp_value = 0;
	index_start = *( mat->row_pointers + index );
	index_end = *( mat->row_pointers + index + 1 );
	for( index_j = index_start; index_j < index_end; index_j++ )
	{
	    temp_value += ( ( mat->nz_values + index_j )->re_part ) 
		* ( *( temp_solution + ( *( mat->column_pointers + index_j ) ) ) );
	}
	*( temp + index ) -= temp_value;
    }

    for( index = 0; index < mat->n_size; index++ )
    {
	temp_value = 0;
	index_start = *( mat->row_pointers + index );
	index_end = *( mat->row_pointers + index + 1 );
	for( index_j = index_start; index_j < index_end; index_j++ )
	{
	    temp_value += ( ( mat->nz_values + index_j )->im_part ) 
		* ( *( temp_solution + mat->n_size + ( *( mat->column_pointers + index_j ) ) ) );
	}
	*( temp + index ) -= temp_value;
    }

    return 0;
}
