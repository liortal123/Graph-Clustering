/**
 * cluster.c
 *
 * DESCRIPTION:
 * This project implements graph clustering according to modularity.
 * This is the main flow of the program.
 * It's functions are implemented in Utility files
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "matA_sparse.h"
#include "Utility.h"
#include "matB.h"

int main(int argc, char *argv[])
{
	int n, *degree_arr, M, *s = NULL, index_P, index_O;
	double norma, *eigenVector = NULL;
	matA_sparse *Ag, *A;
	FILE *fileIn, *fileOut;
	group *g = NULL, **P = NULL, **O = NULL;
	matB *B;

	if (argc<3){
		printf("number of input arguments isn't consistent with the program");
		exit(1);
	}

	/* opening the file */
	n = open(argv, &fileIn, &fileOut, &M);

	/* building and allocating the sparse adjacency matrix A */
	A = building_adjacency_mat(&degree_arr, fileIn, n, M);

	/*creating and allocating struct B*/
	B = matB_allocate(A, degree_arr, M);

	/*creating the initial group g with all the vertices of the graph*/
	g = creating_g(n);

	/* allocating arrays O and P, both used as stacks.
	 * P holds groups that aren't ready and O for the final groups*/
	P = (group**)calloc(n,sizeof(group*));
	check_allocate(P);
	O = (group**)calloc(n,sizeof(group*));
	check_allocate(O);

	/* algorithm 3 */
	P[0] = g;
	index_P = 1; /* the number of groups in P is starting from 0 */
	index_O = 0;

	while (index_P > 0){
		/* removing group g from P */
		g = P[--index_P];
		P[index_P] = NULL;

		/*allocation and creating sub matrix Ag*/
		if (g->size != n){
			Ag = create_sub_spmat(A, g);
			B->Ag = Ag;
		}

		/* norma1 calculation */
		norma = norma_calc(B, g);

		/* finding eigen pair */
		eigenVector = finding_eigen_pair(B, g, norma);

		/* calculating s */
		s = calc_s(B, g, eigenVector);

		/* improving s */
		improving_s(B, g, s);

		/* calculating the division according to s */
		index_P = divide_g(P, O, g, s, &index_O, index_P);

		/*freeing Ag*/
		if (B->Ag->n!=n){
			B->Ag = NULL;
			free_spmat(Ag);
		}


	}

	/* writing the algorithm results -the graph clustering- to fileOut */
	write_to_output_file(O, fileOut, index_O);

	/*freeing all allocation*/
	free_all(O, P, B, A);

	return 0;
}
