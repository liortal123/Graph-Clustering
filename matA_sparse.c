/**
 *
 * matA_sparse.c
 *
 */

#include "matA_sparse.h"
#include "Utility.h"

/* Global array used for inner calculations (allocated once to spare allocations) */
int* A_row = NULL, *Ag_row=NULL;

/* Checking the allocation is successful else ending the program */
void check_allocate_spmat(void *pointer){
	if (pointer==NULL){
		printf("problem with allocation memory");
		exit(1);
	}
}

/* Checking the number of character read/write */
void check_char_num_spmat (int charNum, int expected_char_num, char* message){
	if (charNum!=expected_char_num){
		printf("%s",message);
		exit(1);
	}
}

/* Allocates a new array sparse matrix of size n */
matA_sparse* spmat_allocate_array(int n, int nnz){
	matA_sparse *mat;

	mat = (matA_sparse*)malloc(sizeof(matA_sparse));
	check_allocate_spmat(mat);

	mat->n = n;
	mat->nnz = nnz;
	mat->rowptr = (int*)calloc((n+1),sizeof(int));
	check_allocate_spmat(mat->rowptr);

	if (nnz != 0){
		mat->colind = (int*)malloc(nnz*sizeof(int));
		check_allocate_spmat(mat->colind);
	}

	return mat;
}

/*Allocating memory to sub matrix Ag*/
matA_sparse* allocate_sub_adjacency_mat(matA_sparse *A, group *g){
	int ng, i,j, *A_row_res, nnz = 0;
	matA_sparse *Ag;

	ng = g->size;

	/*find Ag nnz*/
	for (i = 0; i < ng; i++){
		A_row_res = get_spmat_row(A,g->array[i]);
		for (j = 0; j < ng; j++){
			/*for vertex i in group g,
			 * counting only the neighbors that are also vertices in group g*/
			if (A_row_res[g->array[j]] == 1){
				nnz++;
			}
		}
	}

	/*allocation and returning the new Ag as sparse matrix*/
	Ag = spmat_allocate_array(ng,nnz);
	return Ag;
}

/*Creating and allocating the sub matrix Ag out of A*/
matA_sparse* create_sub_spmat(matA_sparse *A, group *g){
	int i, j, ng, *A_row_res;
	matA_sparse* Ag;

	ng = g->size;
	Ag = allocate_sub_adjacency_mat(A, g);

	if (Ag_row == NULL){
		Ag_row = (int*)malloc(A->n*sizeof(int));
		check_allocate_spmat(Ag_row);
	}

	for (i = 0; i < ng; i++){
		/*getting row g[i] of A of size n*/
		A_row_res = get_spmat_row(A,g->array[i]);

		for (j = 0; j < ng; j++){
			/*for vertex i in group g, extracting the sub row of Ag*/
			Ag_row[j] = A_row_res[g->array[j]];
		}

		/*adding the new row to Ag*/
		add_row_spmat(Ag, Ag_row, i);
	}

	return Ag;
}


/* Adds row i the matrix. Called before any other call,
 * exactly n times in order (i = 0 to n-1) */
void add_row_spmat(matA_sparse *mat, int *row, int row_i)
{
	int i, n ,pointer;

	n = mat->n;

	/*pointer symbolize where to add the next nnz to arrays colind and values */
	pointer=mat->rowptr[row_i];
	if (mat->nnz == pointer){
		mat->rowptr[row_i+1] = pointer;
		return;
	}

	/*filling colind array*/
	for (i=0;i<n;i++){
		if (row[i]==0) {continue;}
		mat->colind[pointer] = i;
		pointer++;
	}

	/*update pointer of next row */
	mat->rowptr[row_i+1] = pointer;
}

/* Building and returning the adjacency matrix as sparse matrix with arrays implementation
 * A is the adjacency matrix of the input file graph
 * If two vertices are neighbors the adjacency matrix value will be 1 otherwise 0
 * The function allocates and calculates degree array */
matA_sparse* building_adjacency_mat(int **degree_arr, FILE *fileIn, int n, int M){
	int i, j, readCharNum, *neighbors, *row_arr;
	matA_sparse *A;

	A = spmat_allocate_array(n, M);
	row_arr = (int *)calloc(n,sizeof(int));
	check_allocate_spmat(row_arr);
	neighbors = (int *)malloc(n*sizeof(int));
	check_allocate_spmat(neighbors);
	*degree_arr = (int *)malloc(n*sizeof(int));
	check_allocate_spmat(*degree_arr);

	for (i = 0; i < n; i++){

		readCharNum = fread(&((*degree_arr)[i]), sizeof(int), 1, fileIn);
		check_char_num_spmat(readCharNum,1, "problem with reading input file");

		readCharNum = fread(neighbors, sizeof(int), (*degree_arr)[i], fileIn);
		check_char_num_spmat(readCharNum,(*degree_arr)[i], "problem with reading input file");

		for (j = 0; j < (*degree_arr)[i]; j++){
			/* the adjacency matrix hold 1 for two vertices that are neighbors */
			row_arr[neighbors[j]] = 1;
		}
		add_row_spmat(A,row_arr,i);

		memset(row_arr, 0, n*sizeof(int));
	}

	free(row_arr);
	free(neighbors);
	fclose(fileIn);

	return A;
}

/* The function gets an index of a row i,
 * and returns an array representation of the row with zeros and non zero values */
int* get_spmat_row(matA_sparse  *mat, int row_i){
	int i, start, end, nz_col;

	if (A_row == NULL){
		A_row = (int*)malloc(mat->n*sizeof(int));
		check_allocate_spmat(A_row);
	}

	/*setting A to all zeros*/
	memset(A_row, 0, mat->n*sizeof(int));

	start = mat->rowptr[row_i];
	end = mat->rowptr[row_i+1];

	for (i = start; i < end; i++) {
		nz_col = mat->colind[i]; /*column in the row that is a non zero column */
		A_row[nz_col] = 1;
	}

	return A_row;
}

/* Multiplies matrix Ag by vector v of type double, into result (result is pre-allocated) */
void mult_adj_double(matA_sparse *mat, const double *v, double *result){
	int n,i,j,diff;
	double sum;

	n = mat->n;
	for (i=0;i<n;i++){
		int rowptr = mat->rowptr[i]; /*rowptr symbolize the index in the colind array */
		diff = mat->rowptr[i+1]-rowptr; /*diff symbolize how many nnz elements in row i */
		sum = 0.0;
		for (j=0;j<diff;j++){
			int col_ind = mat->colind[rowptr+j]; /*col_ind symbolize the column of the nnz in the original matrix */

			/*The value in col ind is 1 (Adjacency matrix)
			* so we don't need to multiply v[col_ind]*values[rowptr+j], we just sum*/
			sum += v[col_ind];
		}

		result[i] = sum;
	}
}

/* Multiplies matrix A by vector v of type int, into result (result is pre-allocated) */
void mult_adj_int(matA_sparse *mat, const int *v, double *result){
	int n,i,j,diff;
	int sum;

	n = mat->n;
	for (i=0;i<n;i++){
		int rowptr = mat->rowptr[i]; /*rowptr symbolize the index in the colind array */
		diff = mat->rowptr[i+1]-rowptr; /*diff symbolize how many nnz elements in row i */
		sum = 0.0;
		for (j=0;j<diff;j++){
			int col_ind = mat->colind[rowptr+j]; /*col_ind symbolize the column of the nnz in the original matrix */

			/*The value in col ind is 1 (Adjacency matrix)
			* so we don't need to multiply v[col_ind]*values[rowptr+j], we just sum*/
			sum += v[col_ind];
		}

		result[i] = sum;
	}
}

/* Frees all resources used by A */
void free_spmat(matA_sparse *mat){
	if(mat->nnz != 0){
		free(mat->colind);
	}
	free(mat->rowptr);
	free(mat);
}

/*freeing global arrays*/
void free_global_spmat_arrays(){
	free(A_row);
	free(Ag_row);

}
