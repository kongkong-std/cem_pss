#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../include_file/pss_solver.h"

int pss_solver( const CSR_MAT * mat, const RHS_VEC * b_rhs,
	RHS_VEC * solution, const int max_iter, const double rela_tol )
{
    int index;
    double alpha_pss;

    printf( "\ninput alpha: " );
    scanf( "%lg", &alpha_pss );
    printf( "\nalpha_pss = %lg\n", alpha_pss );

    double * solution_pss;    // solution for PSS solver
    
    // memory allocation for PSS solution
    if( ( solution_pss = ( double * ) malloc( 2 * ( solution->n_size ) * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // initialize solution_pss
    for( index = 0; index < 2 * ( solution->n_size ); index++ )
    {
	*( solution_pss + index ) = 0;
    }

    pss_solver_iteration( mat, b_rhs, alpha_pss, solution_pss, max_iter, rela_tol );
    
    // solution
    for( index = 0; index < solution->n_size; index++ )
    {
	( solution->vec_values + index )->re_part = *( solution_pss + index );
	( solution->vec_values + index )->im_part = -( *( solution_pss + index + solution->n_size ) );
    }

    // free memory
    free( solution_pss );
    
    return 0;
}
