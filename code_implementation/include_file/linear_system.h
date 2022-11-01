#ifndef LINEAR_SYSTEM_H_
#define LINEAR_SYSTEM_H_

// struct definition
/*
 * complex value type
 * contains real part and imaginary part
 * */
typedef struct complex_number{
    double re_part;
    double im_part;
}COMPLEX_VALUE;

/*
 * CSR matrix
 * component:
 *     size of matrix
 *     total number of non-zero elements
 *     index of matrix row
 *     count of non-zero elements in every row
 *     index of non-zero elemnts column
 *     value of non-zero elements
 * */
typedef struct csr_sparse_matrix{
    int n_size;
    int n_nz;
    int * n_row;
    int * row_pointers;
    int * column_pointers;
    COMPLEX_VALUE * nz_values;
}CSR_MAT;

/*
 * right-hand side vector
 * complex version
 * */
typedef struct rhs_vec{
    int n_size;
    COMPLEX_VALUE * vec_values;
}RHS_VEC;

// function prototype
/*
 * free memory of sparse of matrix
 * */
int free_memory_matrix( const CSR_MAT * );

/*
 * free memory of rhs vector
 * */
int free_memory_vector( const RHS_VEC * );

/*
 * print check, print formation for vector
 * */
int print_check_vector( const RHS_VEC * );

/*
 * print check, print information for matrix
 * */
int print_check_matrix( const CSR_MAT * );

/*
 * linear system solver
 * input: matrix( A ), rhs vector( b )
 * function: solving A x = b
 * */
void linear_system_solver( const CSR_MAT *, const RHS_VEC * );
#endif
