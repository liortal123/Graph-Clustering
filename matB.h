/**
 * matB.h
 *
 * DESCRIPTION:
 * This header is implementing mat B - the Modulatirty matrix and it's methods.
 * We don't save the entire matrix B in memory.
 * We hold the adjacency matrix A, degree array and M in struct.
 * We calculate B each time according to the formula: Bij = Aij - KiKj/M
 *
 * matB_allocate - Allocates new mat B - struct that holds A(sparse matrix), degree array and M
 * get_B_row	 - Returns an array representation of row i in B(or B_hat)
 * mult_B_double - Multiplies B and vector of type double
 * mult_B_int	 - Multiplies B and vector of type int
 * free_matB 	 - Frees all memory resources used by matB
 */

#ifndef MATB_H_
#define MATB_H_

#include <stdio.h>
#include <stdlib.h>
#include "matA_sparse.h"

/* All objects necessary in order to create the modulatiry matrix */
typedef struct _matB
{
	int M;
	matA_sparse* Ag;
	int* degree_array;
}matB;

/* Forward declaration of struct group */
struct group;

/* Allocates new mat B
 * We don't save the full matrix B
 * Mat B is composed of matrix A (matA_sparse), degree_array and M*/
matB* matB_allocate(matA_sparse* A, int* degree_array, int M);

/* Calculate B_hat(g) row or B row
 * g is the sub group that we work on (could be the original graph)
 * If g is the full graph then each row sums to zero (therefore B_hat_row calculation is equal to B)
 * The calculation consists of: B[g]=(Aij-(ki*kj)/M) + adding the norma and subtracting sum from the diagonal.
 * if we want to get B without norma we'll pass norma = 0 */
double* get_B_row(const matB* B, const struct group* g, double norma, int row_index);

/* Multiplying B and vector of type double into result
 * We calculat Bij*v = (Aij-ki*kj/M)*vj = Aij*vj - (ki*kj/M)vj +norma*vi
 * if we look at the second part of the calculation:
 * sigma_j((ki*kj/M)*vj)=ki * sigma(kj*vj)/M
 * We'll define constant = sigma(kjvj)/M and calculate it once*/
void mult_B_double(const matB *B, const struct group* g, double norma, const double *vector, double *result);

/* Multiplying B and vector of type int into result
 * We calculat Bij*v = (Aij-ki*kj/M)*vj = Aij*vj - (ki*kj/M)vj +norma*vi
 * if we look at the second part of the calculation:
 * sigma_j((ki*kj/M)*vj)=ki * sigma(kj*vj)/M
 * We'll define constant = sigma(kjvj)/M and calculate it once*/
double* mult_B_int(const matB *B, const struct group* g, double norma, const int *vector);

/* Frees all memory resources used by matB */
void free_matB(matB* B);

#endif /* MATB_H_ */
