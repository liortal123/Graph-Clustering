/**
 * matB.c
 *
 */

#include "matB.h"
#include "Utility.h"

/* Global arrays used for inner calculations (allocated once to spare allocations) */
double *B_row = NULL, *res_mult_B_int = NULL;

/* Checking the allocation is successful else ending the program */
void check_allocate_matB(void *pointer){
	if (pointer==NULL){
		printf("problem with allocation memory");
		exit(1);
	}
}

/* Allocates new mat B
 * We don't save the full matrix B
 * Mat B is composed of matrix A (matA_sparse), degree_array and M*/
matB* matB_allocate(matA_sparse* A, int* degree_array, int M) {
	matB *B;

	B = (matB*)malloc(sizeof(matB));
	check_allocate_matB(B);

	B->Ag = A;
	B->degree_array = degree_array;
	B->M = M;

	return B;
}

/* Calculate B_hat(g) row or B row
 * g is the sub group that we work on (could be the original graph)
 * If g is the full graph then each row sums to zero (therefore B_hat_row calculation is equal to B)
 * The calculation consists of: B[g]=(Aij-(ki*kj)/M) + adding the norma and subtracting sum from the diagonal.
 * if we want to get B without norma we'll pass norma = 0 */
double* get_B_row(const matB* B, const struct group* g, double norma, int row_index){
		int i, *A_row_i, ng, diag=0, g_index, g_index_degree;
		double sum = 0.0, KiKj;

		if (B_row == NULL){
			B_row = (double*)malloc(B->Ag->n*sizeof(double));
			check_allocate_matB(B_row);
		}

		/*getting row i from the sub matrix Ag*/
		ng = g->size;
		g_index = g->array[row_index];
		g_index_degree = B->degree_array[g_index];
		A_row_i = get_spmat_row(B->Ag, row_index);

		for (i = 0; i < ng; i++){
			/* g->array[i] is vertex (column) i in group g */
			KiKj = g_index_degree*B->degree_array[g->array[i]];
			B_row[i] = A_row_i[i] -KiKj/(double)B->M;

			sum += B_row[i]; /* used only when calculating B_hat */

			if (g->array[i] == g_index){ /* adding the norma to the diagonal */
				B_row[i] += norma;
				diag = i;
			}
		}

		/*decreasing the sum for the calculation of B_hat[g]*/
		B_row[diag] -= sum;

		return B_row;
}




/* Multiplying B and vector of type double into result
 * We calculat Bij*v = (Aij-ki*kj/M)*vj = Aij*vj - (ki*kj/M)vj +norma*vi
 * if we look at the second part of the calculation:
 * sigma_j((ki*kj/M)*vj)=ki * sigma(kj*vj)/M
 * We'll define constant = sigma(kjvj)/M and calculate it once*/
void mult_B_double(const matB *B, const group* g, double norma, const double *vector, double *result){
	double constant, sum = 0.0;
	int ng, i;

	ng = g->size;

	/*Multiplying A*vector and inserting it to result*/
	mult_adj_double(B->Ag, vector , result);

	/*Calcultating constant */
	for (i = 0; i < ng; i++){
		sum += B->degree_array[g->array[i]] * vector[i];
	}
	constant = sum / B->M;

	/*Calculation Aij*vi - ki*constant + norma*vi */
	for (i = 0; i < ng; i++){
		result[i] = result[i] - B->degree_array[g->array[i]]*constant + norma*vector[i];
	}

}

/* Multiplying B and vector of type int into result
 * We calculat Bij*v = (Aij-ki*kj/M)*vj = Aij*vj - (ki*kj/M)vj +norma*vi
 * if we look at the second part of the calculation:
 * sigma_j((ki*kj/M)*vj)=ki * sigma(kj*vj)/M
 * We'll define constant = sigma(kjvj)/M and calculate it once*/
double* mult_B_int(const matB *B, const group* g, double norma, const int *vector){
	double constant, sum = 0.0;
	int ng, i;

	if (res_mult_B_int == NULL){
		res_mult_B_int = (double*)malloc(B->Ag->n*sizeof(double));
		check_allocate_matB(res_mult_B_int);
	}

	ng = g->size;

	/*Multiplying A*vector and inserting it to result*/
	mult_adj_int(B->Ag, vector , res_mult_B_int);

	/*Calcultating constant */
	for (i = 0; i < ng; i++){
		sum += B->degree_array[g->array[i]] * vector[i];
	}
	constant = sum / B->M;

	/*Calculation Aij*vi - ki*constant + norma*vi */
	for (i = 0; i < ng; i++){
		res_mult_B_int[i] =res_mult_B_int[i] - B->degree_array[g->array[i]]*constant + norma*vector[i];
	}

	return res_mult_B_int;

}


/* Frees all memory resources used by matB */
void free_matB(matB* B){
	free(res_mult_B_int);
	free(B_row);
	free(B->degree_array);
	free(B);
}
