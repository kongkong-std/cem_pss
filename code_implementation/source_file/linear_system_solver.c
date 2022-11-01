#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../include_file/linear_system_solver.h"

void linear_system_solver( const CSR_MAT * A_mat, const RHS_VEC * b_vec )
{
    dividing_line_1( "dividing line" );
    dividing_line_2( "linear system solver" );

    RHS_VEC * solution;    // solution for complex system
    if( ( solution = ( RHS_VEC * ) malloc( sizeof( RHS_VEC ) ) ) == NULL )
    {
	fprintf( stdout, "Memory failed!\n" );
	exit( EXIT_FAILURE );
    }
#if 1
    // memory allocation for solution
    solution->n_size = b_vec->n_size;
    printf( "==%d==\n", solution->n_size );
    if( ( solution->vec_values = ( COMPLEX_VALUE * ) malloc( solution->n_size * sizeof( COMPLEX_VALUE ) ) )
	    == NULL )
    {
	fprintf( stdout, "Memory failed!\n" );
	exit( EXIT_FAILURE );
    }
#endif

    //==========pss interation==========//
    int max_iter_pss;
    double rela_tol_pss;

    puts( "\nParameters for PSS interation:" );
    printf( "max iterations: " );
    scanf( "%d", &max_iter_pss );
    printf( "relative tolerance: " );
    scanf( "%lg", &rela_tol_pss );

    dividing_line_3( "PSS solver" );
    
    printf( "max_iter\tcriteria\n" );
    printf( "%d\t%.2le\n", max_iter_pss, rela_tol_pss );

    pss_solver( A_mat, b_vec, solution, max_iter_pss, rela_tol_pss );

    //==========solution to file==========//
    FILE * fp;
    int index;
    if( ( fp = fopen( "../result_file/solution_vector.txt", "w" ) ) == NULL )
    {
	fprintf( stdout, "Cannot write file!\n" );
	exit( EXIT_FAILURE );
    }

    fprintf( fp, "%d\n", solution->n_size );
    for( index = 0; index < solution->n_size; index++ )
    {
	fprintf( fp, "%.12le\t%.12le\n", (solution->vec_values + index )->re_part,
		( solution->vec_values + index )->im_part );
    }

    fclose( fp );
    //==========end solution to file==========//

    free_memory_vector( solution );
    free( solution );

    dividing_line_3( "end PSS solver" );
    //==========end pss iteration==========//

    dividing_line_2( "end linear system solver" );
}

void dividing_line_1( const char * str_pointer )
{
    int index;

    putchar( '\n' );
    for( index = 0; index < 40; index++ )
    {
	putchar( '-' );
    }
    printf( "%s", str_pointer );
    for( index = 0; index < 40; index++ )
    {
	putchar( '-' );
    }
    putchar( '\n' );
}

void dividing_line_2( const char * str_pointer )
{
    int index;

    putchar( '\n' );
    for( index = 0; index < 10; index++ )
    {
	putchar( '=' );
    }
    printf( "%s", str_pointer );
    for( index = 0; index < 10; index++ )
    {
	putchar( '=' );
    }
    putchar( '\n' );
}

void dividing_line_3( const char * str_pointer )
{
    int index;

    putchar( '\n' );
    for( index = 0; index < 5; index++ )
    {
	putchar( '=' );
    }
    printf( "%s", str_pointer );
    for( index = 0; index < 5; index++ )
    {
	putchar( '=' );
    }
    putchar( '\n' );
}
