#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../include_file/linear_system.h"

void linear_system( void )
{
    puts( "loading matrix..." );

    FILE * fp;
    CSR_MAT A_sparse;
    RHS_VEC RHS_b;

    int index;

    //==========matrix file==========//
    if( ( fp = fopen( "../linear_system/dof_48192/A_maxwell_lumpedport.txt", "r" ) ) == NULL )
    {
	fprintf( stdout, "Cannot open file!\n" );
	exit( EXIT_FAILURE );
    }

    fscanf( fp, "%d", &A_sparse.n_size );
    fscanf( fp, "%d", &A_sparse.n_nz );

    // row of matrix
    if( ( A_sparse.n_row = ( int * ) malloc( A_sparse.n_size * sizeof( int ) ) ) == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }
    for( index = 0; index < A_sparse.n_size; index++ )
    {
	fscanf( fp, "%d", A_sparse.n_row + index );
    }

    // count of non-zero elements in every row
    if( ( A_sparse.row_pointers = ( int * ) malloc( ( A_sparse.n_size + 1 ) * sizeof( int ) ) ) == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }
    for( index = 0; index < A_sparse.n_size + 1; index++ )
    {
	fscanf( fp, "%d", A_sparse.row_pointers + index );
    }

    // index of non-zero elements column
    if( ( A_sparse.column_pointers = ( int * ) malloc( A_sparse.n_nz * sizeof( int ) ) ) == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }
    for( index = 0; index < A_sparse.n_nz; index++ )
    {
	fscanf( fp, "%d", A_sparse.column_pointers + index );
    }
#if 1
    // value of non-zero elements
    if( ( A_sparse.nz_values = ( COMPLEX_VALUE * ) malloc( A_sparse.n_nz * sizeof( COMPLEX_VALUE ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }
    for( index = 0; index < A_sparse.n_nz; index++ )
    {
	fscanf( fp, "%lg %lg", &( ( A_sparse.nz_values + index)->re_part ),
		&( ( A_sparse.nz_values + index )->im_part ) );
    }
#endif

    fclose( fp );
    puts( "\nmatrix loaded !!!" );
    //==========end matrix file==========//
#ifdef CHECK_MAT_DEBUG
    print_check_matrix( &A_sparse );    // input check for matrix
#endif

    putchar( '\n' );
    for( index = 0; index < 40; index++ )
    {
	putchar( '-' );
    }
    printf( "%s", "dividing line" );
    for( index = 0; index < 40; index++ )
    {
	putchar( '-' );
    }
    putchar( '\n' );

    //==========rhs vector file==========//
    puts( "\nloading vector ..." );
    if( ( fp = fopen( "../linear_system/dof_48192/RHS_maxwell_lumpdeport.txt", "r" ) ) == NULL )
    {
	fprintf( stdout, "Cannot open file!\n" );
	exit( EXIT_FAILURE );
    }

    fscanf( fp, "%d", &RHS_b.n_size );

    // memory allocation for right-hand side vector
    if( ( RHS_b.vec_values = ( COMPLEX_VALUE * ) malloc( RHS_b.n_size * sizeof( COMPLEX_VALUE ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    for( index = 0; index < RHS_b.n_size; index++ )
    {
	fscanf( fp, "%le %le", &( ( RHS_b.vec_values + index )->re_part ),
		&( ( RHS_b.vec_values + index )->im_part ) );
    }

    fclose( fp );
    puts( "\nvector loaded !!!" );
    //==========end rhs vector file==========//
#if 0
    print_check_vector( &RHS_b );    // input check for vector
#endif

    putchar( '\n' );
    for( index = 0; index < 10; index++ )
    {
	putchar( '=' );
    }
    printf( "%s", "linear system loaded" );
    for( index = 0; index < 10; index++ )
    {
	putchar( '=' );
    }
    putchar( '\n' );

    // solving linear system
    linear_system_solver( &A_sparse, &RHS_b );

    //==========free memory==========//
    free_memory_matrix( &A_sparse );
    free_memory_vector( &RHS_b );
}

int print_check_vector( const RHS_VEC * vec )
{
    int index;

    printf( "\nn_size\n" );
    printf( "%d\n", vec->n_size );

    printf( "\nre_part\tim_part\n" );
    for( index = 34632; index < 34732; index++ )
    {
	printf( "%.12e\t%.12e\n", ( vec->vec_values + index )->re_part,
		( vec->vec_values + index )->im_part );
    }

    return EXIT_SUCCESS;
}

int print_check_matrix( const CSR_MAT * mat )
{
    printf( "n_size\tn_nz\n" );
    printf( "%d\t%d\n", mat->n_size, mat->n_nz );

    int index;
    printf( "re_part\tim_part\n" );
    for( index = 0; index < 100; index++ )
    {
	printf( "%.12e\t%.12e\n", ( mat->nz_values + index)->re_part,
		( mat->nz_values + index )->im_part );
    }

    return EXIT_SUCCESS;
}

int free_memory_matrix( const CSR_MAT * mat )
{
    free( mat->n_row );
    free( mat->row_pointers );
    free( mat->column_pointers );
    free( mat->nz_values );

    puts( "!!!! Free Memory Successfully !!!!" );

    return EXIT_SUCCESS;
}

int free_memory_vector( const RHS_VEC * vec )
{
    free( vec->vec_values );

    puts( "!!!! Free Memory Successfully !!!!" );
    
    return EXIT_SUCCESS;
}
