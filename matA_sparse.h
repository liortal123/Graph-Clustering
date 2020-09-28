/**
 * matA_sparse.h
 *
 * DESCRIPTION:
 * This header is implementing mat A - the adjacency matrix of the given graph and it's methods.
 * Mat A is a sparse matrix implemented with arrays.
 *
 * spmat_allocate_array     - Allocates new sparse matrix implemented with arrays
 * add_row_spmat		    - Adds row i to the sparse matrix
 * building_adjacency_mat   - Builds adjacency matrix from the input file graph
 * create_sub_spmat	        - Creates the sub matrix Ag out of A
 * get_spmat_row		    - Returns an array representation of row i in A
 * mult_adj_double		    - Multiplies Ag and vector of type double
 * mult_adj_int			    - Multiplies Ag and vector of type int
 * free_spmat			    - Frees all memory resources used by mat
 * free_global_spmat_arrays - Frees global arrays used in matA sparse
 */

#ifndef MATA_SPARSE
#define MATA_SPARSE

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>

typedef struct matA_sparse {
	 int *colind,*rowptr, nnz;
	 int n;
} matA_sparse;

/* Forward declaration of struct group */
struct group;

/* Allocates a new array sparse matrix of size n */
matA_sparse* spmat_allocate_array(int n, int nnz);

/* Adds row i to the matrix. Called before any other call,
 * exactly n times in order (i = 0 to n-1) */
void add_row_spmat(matA_sparse *mat, int *row, int row_i);

/* Building and returning the adjacency matrix as sparse matrix with arrays implementation
 * A is the adjacency matrix of the input file graph
 * If two vertices are neighbors the adjacency matrix value will be 1 otherwise 0
 * The function allocates and calculates degree array */
matA_sparse* building_adjacency_mat(int **degree_arr, FILE *fileIn, int n, int M);

/*Creating and allocating the sub matrix Ag out of A*/
matA_sparse* create_sub_spmat(matA_sparse *A, struct group *g);

/* The function gets an index of a row i,
 * and returns an array representation of the row with zeros and non zero values */
int* get_spmat_row(matA_sparse  *mat, int row_i);

/* Multiplies matrix Ag by vector v of type double, into result (result is pre-allocated) */
void mult_adj_double(matA_sparse *mat, const double *v, double *result);

/* Multiplies matrix Ag by vector v of type int, into result (result is pre-allocated) */
void mult_adj_int(matA_sparse *mat, const int *v, double *result);

/* Frees all resources used by mat */
void free_spmat(matA_sparse  *mat);

/*freeing global arrays*/
void free_global_spmat_arrays();

#endif /* SPMAT_H_ */
