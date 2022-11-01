#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../include_file/gmres_restart_iteration.h"

// static function prototype
static int CopyItemToNode( ARNOLDI_VEC, ARNOLDI_VEC_NODE * );

int gmres_restart_iteration( double * solution,
	const CSR_MAT * mat, const double alpha_pss,
	const double * rhs_vector, const int restart_number )
{
    int index, index_i;
    double * residual_solution;
    double norm_residual_solution;

    if( ( residual_solution = ( double * ) malloc( 2 * ( mat->n_size ) * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // computing residual_solution = gmres_rhs - gmres_mat * solution
    gmres_computing_residual_vector( residual_solution,
	    mat, alpha_pss,
	    rhs_vector, solution );
    norm_residual_solution = computing_vector_norm_2( residual_solution, 2 * mat->n_size );

    // linked list
    /*
     * Arnoldi decomposition
     * A Q_k = Q_{k+1} \hat H_{ k + 1 }
     * A [ q1, q2, ..., q_k ] = [ q_1, ..., q_k, q_{k + 1} ] \hat H_{k + 1}
     * */
    ARNOLDI_VEC_LIST arnoldi_matrix;
    ARNOLDI_VEC vector_temp;

    // initialize linked list
    InitializeLinkedList( &arnoldi_matrix );

    // determine whether linked list is full or not
    if( IsLinkedListFull( &arnoldi_matrix ) )
    {
	fprintf( stdout, "Linked list is full!\n" );
	exit( EXIT_FAILURE );
    }

    // create linked list element
    vector_temp.n_size = 2 * mat->n_size;
    if( ( vector_temp.vec_values = ( double * ) malloc( vector_temp.n_size * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stdout, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    for( index = 0; index < vector_temp.n_size; index++ )
    {
	*( vector_temp.vec_values + index ) = ( *( residual_solution + index ) ) 
	    / norm_residual_solution ;
    }

    // add element to linked list
    AddItemToLinkedList( vector_temp, &arnoldi_matrix );

    // GMRES restart iteration
    ARNOLDI_VEC_NODE * v_vector = arnoldi_matrix.head_node;
    ARNOLDI_VEC_NODE * v_vector_temp = arnoldi_matrix.head_node;
    double * * gmres_up_hessenberg;    // up hessenberg matrix, ( m + 1 ) * m

    // memory allocation for Hessenberg matrix
    if( ( gmres_up_hessenberg = ( double * * ) malloc( ( restart_number + 1 ) * sizeof( double * ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }
    for( index = 0; index < restart_number + 1; index++ )
    {
	if( ( *( gmres_up_hessenberg + index ) = ( double * )
		    malloc( restart_number * sizeof( double ) ) ) 
		== NULL )
	{
	    fprintf( stderr, "Memory allocation failed!\n" );
	    exit( EXIT_FAILURE );
	}
    }

    // initialize for Hessenberg matrix
    for( index_i = 0; index_i < restart_number + 1; index_i++ )
    {
	for( index = 0; index < restart_number; index++ )
	{
	    *( *( gmres_up_hessenberg + index_i ) + index ) = 0;
	}
    }

    for( index = 0; index < restart_number; index++ )
    {
	// omega_j = gmres_mat * v_j
	/*
	 * re-use residual_solution
	 * */
	gmres_restart_computing_mvp( residual_solution,
		mat, alpha_pss, v_vector );

	/*
	 * h_ij = ( v_i, w_j )
	 * w_j = w_j - h_ij v_i
	 * */
	v_vector_temp = arnoldi_matrix.head_node;
	for( index_i = 0; index_i < index; index_i++ )
	{
	    *( *( gmres_up_hessenberg + index_i ) + index ) =
		gmres_restart_inner_product( residual_solution,
			v_vector_temp );
	    gmres_restart_updating_w_vector( residual_solution,
		    v_vector_temp,
		    *( *( gmres_up_hessenberg + index_i ) + index ) );
	    v_vector_temp = v_vector_temp->next;
	}

	// h_{ j + 1, j } = || w_vector ||
	*( *( gmres_up_hessenberg + index + 1 ) + index ) = 
	    computing_vector_norm_2( residual_solution, 2 * mat->n_size );

	// vector_temp = w_vector / h _ { j + 1, j }
	for( int index_temp = 0; index_temp < vector_temp.n_size; index_temp++ )
	{
	    *( vector_temp.vec_values + index_temp ) = ( *( residual_solution + index_temp ) ) 
		/ ( *( *( gmres_up_hessenberg + index + 1 ) + index ) );
	}

	// add vector_temp to linked list
        AddItemToLinkedList( vector_temp, &arnoldi_matrix );
	v_vector = v_vector->next;
    }

    // computing least-square equation
    /*
     * y = argmin || \beta e_1 - \hat H_m y_m || _ 2
     * */
    double * lse_para;

    // memory allocation
    if( ( lse_para = ( double * ) malloc( restart_number * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    gmres_restart_lse( lse_para, gmres_up_hessenberg,
	    restart_number, norm_residual_solution );

    // check up Hessenberg matrix
#if 0
    FILE * fp;
    if( ( fp = fopen( "../result_file/up_hessenberg.txt", "w" ) ) 
	    == NULL )
    {
	fprintf( stderr, "Cannot open file!\n" );
	exit( EXIT_FAILURE );
    }

    for( index_i = 0; index_i < restart_number + 1; index_i++ )
    {
	for( index = 0; index < restart_number; index++ )
	{
	    fprintf( fp, "%.2f\t", *( *( gmres_up_hessenberg + index_i ) + index ) );
	}
	fprintf( fp, "\n" );
    }

    fclose( fp );
#endif

    // updating solution
    /*
     * solution = solution + [ v_1, ..., v_m ] y
     *          = solution + ( y_1 v_1 + y_2 v_2 + ... + y_m v_m )
     * */
    gmres_restart_updating_solution_linear_system( solution, lse_para, &arnoldi_matrix,
	    2 * ( mat->n_size ), restart_number );

    // empty linked lise
    EmptyLinkedList( &arnoldi_matrix );

    // free memory
    free( residual_solution );
    free( vector_temp.vec_values );
    for( index = 0; index < restart_number + 1; index++ )
    {
	free( *( gmres_up_hessenberg + index ) );
    }
    free( gmres_up_hessenberg );
    free( lse_para );

    return 0;
}

int EmptyLinkedList( ARNOLDI_VEC_LIST * list_pointer )
{
    ARNOLDI_VEC_NODE * node_current;

    while( list_pointer->head_node != NULL )
    {
	node_current = ( list_pointer->head_node )->next;
	free( list_pointer->head_node );
	list_pointer->head_node = node_current;
    }

    list_pointer->list_size = 0;

    return 0;
}

int gmres_restart_updating_solution_linear_system( double * solution, const double * combination,
	const ARNOLDI_VEC_LIST * list_pointer, int solution_size, int combination_size )
{
    ARNOLDI_VEC_NODE * vector_temp = list_pointer->head_node;
    int iter, index_i;

    for( iter = 0; iter < combination_size; iter++ )
    {
	for( index_i = 0; index_i < solution_size; index_i++ )
	{
	    *( solution + index_i ) += ( *( combination + iter ) ) 
		* ( *( ( vector_temp->vec_item ).vec_values + index_i ) );
	}
	vector_temp = vector_temp->next;
    }

    return 0;
}

int gmres_restart_lse( double * lse_solution, const double * * h_mat,
	const int size, const double beta )
{
    int index_i, index_j;
    double * vec_lse;

    if( ( vec_lse = ( double * ) malloc( ( size + 1 ) * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // set value for vec_lse
    /*
     * vec_lse( 1 ) = beta
     * vec_lse( i ) = 0, i \ne 1
     * */
    for( index_i = 0; index_i < size + 1; index_i++ )
    {
	if( index_i == 0 )
	{
	    *( vec_lse + index_i ) = beta;
	}
	else
	{
	    *( vec_lse + index_i ) = 0;
	}
    }

    // householder transformation
    /*
     * h_mat = q_mat * r_mat
     * size h_mat = size + 1 * size
     * size q_mat = size + 1 * size + 1, q_mat ' x q_mat = identity matrix
     * size r_mat = size + 1 x size, r_mat is upper triangular matrix
     * || vec_lse - h_mat y || = || q_mat' vec_lse - r_mat y ||
     * */
    gmres_restart_householder_lse( lse_solution, vec_lse, h_mat, size );

    // free memory
    free( vec_lse );
}

int gmres_restart_householder_lse( double * vec_y,
	const double * vec_lse, const double * * mat, const int size )
{
    int m_size, n_size;
    int iter;
    int index_i, index_j, index_k;
    double temp_value;

    // initialize m_size, n_size
    /*
     * m_size = size + 1
     * n_size = size
     * */
    m_size = size + 1;
    n_size = size;

    // double * * q_mat;
    double * * r_mat;

#if 0
    // size q_mat = m_size x m_size
    if( ( q_mat = ( double * * ) malloc( m_size * sizeof( double * ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }
    for( index_i = 0; index_i < m_size; index_i++ )
    {
	if( ( *( q_mat + index_i ) = ( double * ) malloc( m_size * sizeof( double ) ) ) 
		== NULL )
	{
	    fprintf( stderr, "Memory allocation failed!\n" );
	    exit( EXIT_FAILURE );
	}
    }

    // initialize q_mat
    /*
     * q_mat = identity matrix
     * q_mat( i, i ) = 1
     * q_mat( i, j ) = 0, i \ne j
     * */
    for( index_i = 0; index_i < m_size; index_i++ )
    {
	for( index_j = 0; index_j < m_size; index_j++ )
	{
	    if( index_i == index_j )
	    {
		*( *( q_mat + index_i ) + index_j ) = 1;
	    }
	    else
	    {
		*( *( q_mat + index_i ) + index_j ) = 0;
	    }
	}
    }
#endif

    // size r_mat = m_size x n_size
    if( ( r_mat = ( double * * ) malloc( m_size * sizeof( double * ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }
    for( index_i = 0; index_i < m_size; index_i++ )
    {
	if( ( *( r_mat + index_i ) = ( double * ) malloc( n_size * sizeof( double ) ) ) 
		== NULL )
	{
	    fprintf( stderr, "Memory allocation failed!\n" );
	    exit( EXIT_FAILURE );
	}
    }

    // initialize r_mat = mat
    /*
     * r_mat( i, j ) = mat( i, j )
     * */
    for( index_i = 0; index_i < m_size; index_i++ )
    {
	for( index_j = 0; index_j < n_size; index_j++ )
	{
	    *( *( r_mat + index_i ) + index_j ) = *( *( mat + index_i ) + index_j );
	}
    }

    // householder transformation
    double * * h_mat;    // householder matirx
    double * vec_temp;

    // size h_mat = m_size x m_size
    if( ( h_mat = ( double * * ) malloc( m_size * sizeof( double * ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }
    for( index_i = 0; index_i < m_size; index_i++ )
    {
	if( ( *( h_mat + index_i ) = ( double * ) malloc( m_size * sizeof( double ) ) ) 
		== NULL )
	{
	    fprintf( stderr, "Memory allocation failed!\n" );
	    exit( EXIT_FAILURE );
	}
    }

    // size vec_temp = m_size x 1
    if( ( vec_temp = ( double * ) malloc( m_size * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    for( iter = 0; iter < m_size - 1; iter++ )
    {
	// initialize vec_temp
	/*
	 * vec_temp = 0
	 * */
	for( index_i = 0; index_i < m_size; index_i++ )
	{
	    *( vec_temp + index_i ) = 0;
	}

	// computing norm r_mat( iter : m_size, iter )
	temp_value = gmres_restart_computing_vector_norm_part( r_mat, iter, m_size, iter );

	// computing vec_temp
	/* 
	 * vec_temp( iter : m_size ) = 
	 * r_mat( iter : m_size, iter ) - temp_value * e( iter : m_size )
	 * */
	for( index_i = iter; index_i < m_size; index_i++ )
	{
	    if( index_i == iter )
	    {
		*( vec_temp + index_i ) = *( *( r_mat + index_i ) + iter ) - temp_value;
	    }
	    else
	    {
		*( vec_temp + index_i ) = *( *( r_mat + index_i ) + iter );
	    }
	}

	// computing norm of vec_temp
	temp_value = computing_vector_norm_2( vec_temp, m_size );

	if( temp_value == 0 )
	{
	    continue;
	}
	else
	{
	    // computing householder matrix h_mat
	    /*
	     * h_mat = I - 2 vec_temp vec_temp'
	     * */
	    gmres_restart_computing_householder_matrix( h_mat, vec_temp, m_size );

	    // updating r_mat, vec_lse
	    /*
	     * r_mat = h_mat r_mat
	     * vec_lse = h_mat vec_lse
	     * */
	    gmres_restart_updating_lse( vec_lse, r_mat, h_mat, m_size, n_size );
	}
    }

    // solveing linear equation vec_lse( 1 : n_size ) = r_mat( 1 : n_size, 1 : n_size ) vec_y
    /*
     * r_mat is triangular matrix
     * */
    gmres_restart_updating_solution_lse( vec_y, r_mat, vec_lse, n_size );

    // free memory
#if 0
    for( index_i = 0; index_i < m_size; index_i++ )
    {
	free( *( q_mat + index_i ) );
    }
    free( q_mat );
#endif

    for( index_i = 0; index_i < m_size; index_i++ )
    {
	free( *( r_mat + index_i ) );
    }
    free( r_mat );

    for( index_i = 0; index_i < m_size; index_i++ )
    {
	free( *( h_mat + index_i ) );
    }
    free( h_mat );
}

int gmres_restart_updating_solution_lse( double * solution,
	const double * * mat, const double * rhs, int size )
{
    int index_i, index_j;
    double temp_value;

    for( index_i = size - 1; index_i >= 0; index_i-- )
    {
	temp_value = 0;
	for( index_j = index_i + 1; index_j < size; index_j++ )
	{
	    temp_value += ( *( *( mat + index_i ) + index_j ) ) 
		* ( *( solution + index_j ) );
	}
	*( solution + index_i ) = ( ( *( rhs + index_i ) ) - temp_value ) 
	    / ( *( *( mat + index_i ) + index_i ) );
    }

    return 0;
}

int gmres_restart_updating_lse( double * vec, double * * r_mat,
	const double * * h_mat, int m_size, int n_size )
{
    double * vec_temp;
    double * * r_mat_temp;
    int index_i, index_j, index_k;
    double temp_value;

    if( ( vec_temp = ( double * ) malloc( m_size * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    if( ( r_mat_temp = ( double * * ) malloc( m_size * sizeof( double * ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }
    for( index_i = 0; index_i < m_size; index_i++ )
    {
	if( ( *( r_mat_temp + index_i ) = ( double * ) malloc( n_size * sizeof( double ) ) ) 
		== NULL )
	{
	    fprintf( stderr, "Memory allocation failed!\n" );
	    exit( EXIT_FAILURE );
	}
    }

    // initialize vec_temp
    for( index_i = 0; index_i < m_size; index_i++ )
    {
	*( vec_temp + index_i ) = *( vec + index_i );
    }

    // vec_temp = h_mat * vec
    /*
     * vec_temp( i ) = \sum k=0^m h_mat( i, k ) vec( k )
     * */
    for( index_i = 0; index_i < m_size; index_i++ )
    {
	temp_value = 0;
	for( index_j = 0; index_j < m_size; index_j++ )
	{
	    temp_value += ( *( *( h_mat + index_i ) + index_j ) ) 
		* ( *( vec + index_j ) );
	}
	*( vec_temp + index_i ) = temp_value;
    }

    // vec = vec_temp
    for( index_i = 0; index_i < m_size; index_i++ )
    {
	*( vec + index_i ) = *( vec_temp + index_i );
    }

    // initialize r_mat_temp
    /*
     * r_mat_temp = h_mat * r_mat
     * r_mat_temp( i, j ) = \sum k=0^m h_mat( i, k ) r_mat( k, j )
     * */
    for( index_i = 0; index_i < m_size; index_i++ )
    {
	for( index_j = 0; index_j < n_size; index_j++ )
	{
	    temp_value = 0;
	    for( index_k = 0; index_k < m_size; index_k++ )
	    {
		temp_value += ( *( *( h_mat + index_i ) + index_k ) ) 
		    * ( *( *( r_mat + index_k ) + index_j ) );
	    }
	    *( *( r_mat_temp + index_i ) + index_j ) = temp_value;
	}
    }

    // r_mat = r_mat_temp
    for( index_i = 0; index_i < m_size; index_i++ )
    {
	for( index_j = 0; index_j < n_size; index_j++ )
	{
	    *( *( r_mat + index_i ) + index_j ) = *( *( r_mat_temp + index_i ) + index_j );
	}
    }

    // free memory
    free( vec_temp );
    for( index_i = 0; index_i < m_size; index_i++ )
    {
	free( *( r_mat_temp + index_i ) );
    }
    free( r_mat_temp );

    return 0;
}

int gmres_restart_computing_householder_matrix( double * * mat, const double * vec, int size )
{
    int index_i, index_j;

    // computing mat = -2 vec vec'
    /*
     * mat( i, j ) = -2 vec( i ) vec( j )
     * */
    for( index_i = 0; index_i < size; index_i++ )
    {
	for( index_j = 0; index_j < size; index_j++ )
	{
	    *( *( mat + index_i ) + index_j ) = -2 
		* ( *( vec + index_i ) ) * ( *( vec + index_j ) );
	}
    }

    // computing mat = I + mat
    /*
     * mat( i, i ) = 1 + mat( i, i )
     * */
    for( index_i = 0; index_i < size; index_i++ )
    {
	*( *( mat + index_i ) + index_i ) += 1;
    }

    return 0;
}

double gmres_restart_computing_vector_norm_part( const double * * mat,
	int index_start, int index_end, int column )
{
    double value;
    int index;

    // initialize value
    value = 0;

    for( index = index_start; index < index_end; index++ )
    {
	value += ( *( *( mat + index ) + column ) ) 
	    * ( *( *( mat + index ) + column ) );
    }

    return sqrt( value );
}

int gmres_restart_updating_w_vector( double * w_vector,
	const ARNOLDI_VEC_NODE * v_vector, const double h_ele )
{
    int index;

    for( index = 0; index < ( v_vector->vec_item ).n_size; index++ )
    {
	*( w_vector + index ) -= h_ele * 
	    ( *( ( v_vector->vec_item ).vec_values + index ) );
    }

    return 0;
}

double gmres_restart_inner_product( const double * w_vector, const ARNOLDI_VEC_NODE * v_vector )
{
    double value;
    int index;

    // initialize 
    value = 0;

    for( index = 0; index < ( v_vector->vec_item ).n_size; index++ )
    {
	value += ( *( w_vector + index ) ) * 
	    ( *( ( v_vector->vec_item ).vec_values + index ) );
    }

    return value;
}

int gmres_restart_computing_mvp( double * result_vector,
	const CSR_MAT * mat, const double alpha,
	const ARNOLDI_VEC_NODE * node_pointer )
{
    int index;
    double * temp_vec_1, * temp_vec_2;

    // memory allocation for temporary vector
    if( ( temp_vec_1 = ( double * ) malloc( mat->n_size * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    if( ( temp_vec_2 = ( double * ) malloc( mat->n_size * sizeof( double ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	exit( EXIT_FAILURE );
    }

    // temp_vec_1 = alpha I * v( 1 : n ) - W * v( n + 1 : 2n )
    gmres_computing_residual_vector_part_1( temp_vec_1,
	    mat, alpha, ( node_pointer->vec_item ).vec_values );

    // temp_vec_2 = w * v( 1 : n ) + ( alpha I + T ) * v( n + 1 : 2n )
    gmres_computing_residual_vector_part_2( temp_vec_2,
	    mat, alpha, ( node_pointer->vec_item ).vec_values );

    // updating result_vector
    for( index = 0; index < mat->n_size; index++ )
    {
	*( result_vector + index ) = *( temp_vec_1 + index );
	*( result_vector + mat->n_size + index ) = *( temp_vec_2 + index );
    }

    // free memory
    free( temp_vec_1 );
    free( temp_vec_2 );

    return 0;
}

int AddItemToLinkedList( ARNOLDI_VEC item, ARNOLDI_VEC_LIST * list_pointer )
{
    int value;

    ARNOLDI_VEC_NODE * new_node_pointer;    // new node
    ARNOLDI_VEC_NODE * node_pointer = list_pointer->head_node;    // current node

    // memory allocation for new_node_pointer
    if( ( new_node_pointer = ( ARNOLDI_VEC_NODE * ) malloc( sizeof( ARNOLDI_VEC_NODE ) ) ) 
	    == NULL )
    {
	fprintf( stderr, "Memory allocation failed!\n" );
	value = 0;
	return value;
    }

    // updating size of linked list
    ( list_pointer->list_size )++;

    // copy item to new_node_pointer
    CopyItemToNode( item, new_node_pointer );

    // add node to tail of linked list
    new_node_pointer->next = NULL;
    if( node_pointer == NULL )
    {
	list_pointer->head_node = new_node_pointer;
    }
    else
    {
	while( node_pointer->next != NULL )
	{
	    node_pointer = node_pointer->next;
	}
	node_pointer->next = new_node_pointer;
    }

    value = 1;

    return value;
}

static int CopyItemToNode( ARNOLDI_VEC item, ARNOLDI_VEC_NODE * node_pointer )
{
    // copying struct
    node_pointer->vec_item = item;

    return 0;
}

int InitializeLinkedList( ARNOLDI_VEC_LIST * list_pointer )
{
    list_pointer->head_node = NULL;
    list_pointer->list_size = 0;

#if 0
    puts( "====Linked list initialize successfully!===\n" );
#endif

    return 0;
}

int IsLinkedListFull( const ARNOLDI_VEC_LIST * list_pointer )
{
    int value;

    ARNOLDI_VEC_NODE * node_pointer;

    if( ( node_pointer = ( ARNOLDI_VEC_NODE * ) malloc( sizeof( ARNOLDI_VEC_NODE ) ) ) 
	    == NULL )
    {
	value = 1;
	fprintf( stdout, "Linked list is full!\n\n" );
    }
    else
    {
	value = 0;
	fprintf( stdout, "Linked list is not full!\n\n" );
    }

    return value;
}
