#ifndef GMRES_RESTART_ITERATION_H_
#define GMRES_RESTART_ITERATION_H_

// include file
#include "../include_file/pss_gmres_solver.h"
#include <math.h>

// linked list
/*
 * Arnoldi decomposition Q = [ q1, q2, ... qk ]
 * */
typedef struct arnoldi_vector{
    int n_size;
    double * vec_values;
}ARNOLDI_VEC;

/*
 * node of linked list
 * */
typedef struct arnoldi_vector_node{
    ARNOLDI_VEC vec_item;
    struct arnoldi_vector_node * next;
}ARNOLDI_VEC_NODE;

/*
 * linked list
 * */
typedef struct arnoldi_vector_list{
    ARNOLDI_VEC_NODE * head_node;
    int list_size;
}ARNOLDI_VEC_LIST;

// function prototype
/*
 * linked list function
 * initialize linked list
 * plist->head_node = null, plist->list_size = 0
 * return: EXIT_SUCCESS
 * */
int InitializeLinkedList( ARNOLDI_VEC_LIST * );

/*
 * linked list function
 * determine whether linked list is full or not
 * if full, return 1; else, return 0
 * */
int IsLinkedListFull( const ARNOLDI_VEC_LIST * );

/*
 * linked list function
 * add item to linked list
 * return: success or fail
 * */
int AddItemToLinkedList( ARNOLDI_VEC, ARNOLDI_VEC_LIST * );

/*
 * linked list function
 * empty linked list
 * */
int EmptyLinkedList( ARNOLDI_VEC_LIST * );

/*
 * gmres restart function
 * computing omega_j = A v_j
 * input: gmres_mat, linked list node v_j
 * output: result of matrix-by-vector
 * return: EXIT_SUCCESS
 * */
int gmres_restart_computing_mvp( double *,
	const CSR_MAT *, const double, const ARNOLDI_VEC_NODE * );

/*
 * gmres restart function
 * computing inner product of vectors h_ij = ( w_j, v_i )
 * input: vector w_j, linked list node v_vector
 * return: value of inner product
 * */
double gmres_restart_inner_product( const double *, const ARNOLDI_VEC_NODE * );

/*
 * gmres restart function
 * updating w_vector w_j = w_j - h_ij v_i
 * input: linkes list node v_vector, h_ij
 * output: w_vector
 * return: EXIT_SUCCESS
 * */
int gmres_restart_updating_w_vector( double *, const ARNOLDI_VEC_NODE *, const double );

/*
 * gmres restart function
 * solving least-square equation: y = argmin || \beta e_1 - H_m y_m ||
 * input: up Hessenberg matrix H, size of matrix, \beta
 * output: solution y_m
 * return: EXIT_SUCCESS
 * */
int gmres_restart_lse( double *, const double * *,
	const int, const double );

/*
 * gmres restart function
 * solving least-square equation with householder transformation: || \beat e_1 - H_m y ||
 * H_m = q_m * r_m, q_m' x q_m = identity matrix, r_m is upper triangular matrix
 * size H_m = m_size + 1 x m_size, size q_m = m_size + 1 x m_size + 1, size r_mat = m_size + 1 x m_size
 * || \beta e_1 - H_m y || = || q_m' * \beta e_1 - r_m y ||
 * solving { q_m' * \beta e_1 }( 1 : m_size ) = { r_m y }( 1 : m_size )
 * */
int gmres_restart_householder_lse( double *, const double *, const double * *, const int );

/*
 * gmres restart function
 * norm of mat( row_start : row_end, column )
 * */
double gmres_restart_computing_vector_norm_part( const double * *, int, int, int );

/*
 * gmres restart function
 * h_mat = I - 2 v v'
 * */
int gmres_restart_computing_householder_matrix( double * *, const double *, int );

/*
 * gmres restart function
 * LSE || vec_lse - mat y ||
 * vec_lse = h_mat vec_lse
 * mat = h_mat mat
 * */
int gmres_restart_updating_lse( double *, double * *, const double * *, int, int );

/*
 * gmres restart function
 * Gauss elimination solves linear system: A x = b
 * A is upper triangular matrix
 * */
int gmres_restart_updating_solution_lse( double *, const double * *, const double *, int );

/*
 * gmres restart function
 * x_m = x_0 + [ v_1, v_2, ..., v_m ] y
 *     = x_0 + ( y_1 v_1 + y_2 v_2 + ... + y_m v_m )
 * */
int gmres_restart_updating_solution_linear_system( double *, const double *, const ARNOLDI_VEC_LIST *,
	int, int );

#endif
