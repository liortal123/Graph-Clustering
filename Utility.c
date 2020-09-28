/**
 * Utility.c
 *
 */

#include "Utility.h"

#define epsilon 0.00001

/* Global arrays used for inner calculations (allocated once to spare allocations) */
 double *vector = NULL, *nextVector = NULL, *eigenVector_Mult_B = NULL;
 double *improve = NULL, *score = NULL;
 int *s = NULL, *unmoved = NULL, *indices = NULL;

 /* Checking the allocation is successful else ending the program */
 void check_allocate(void *pointer){
 	if (pointer==NULL){
 		printf("problem with allocation memory");
 		exit(1);
 	}
 }

 /* Checking the number of character read/write */
 void check_char_num (int charNum, int expected_char_num, char* message){
 	if (charNum!=expected_char_num){
 		printf("%s",message);
 		exit(1);
 	}
 }

 /* Opening and reading the size of the adjacency matrix and calculating nnz
  * Input file: the first value is n=|V|, the function returns it */
 int open(char *argv[], FILE **fileIn, FILE **fileOut, int *nnz){
 	int readCharNum, n, sumOfChars;

 	*fileIn = fopen(argv[1], "r");
 	*fileOut = fopen(argv[2], "w");
 	if(*fileIn==NULL || *fileOut==NULL){
 		printf("problem with opening input file/ output file");
 		exit(1);
 	}

 	/*counting how many chars in the file */
 	fseek(*fileIn, 0, SEEK_END);
 	sumOfChars = ftell(*fileIn)/sizeof(int);
 	rewind(*fileIn);	/* rewinding the file to the beginning */

 	/*reading the first value n*/
 	readCharNum = fread(&n, sizeof(int), 1, *fileIn);
 	check_char_num(readCharNum,1, "problem with reading input file");

 	/*calculating M = nnz.
 	 * nnz = sumOfChars - n (number of vertices) - 1 (the first value in file is n)*/
 	*nnz = sumOfChars - n - 1;

 	if (*nnz<=0){
 		printf("M=0 zero division problem");
 		exit(1);
 	}

 	return n;
 }

/* Creating the initial group g with all the vertices of the graph */
group *creating_g(int n){
	int i;
	group *g;

	g = (group*)malloc(sizeof(group));
	check_allocate(g);

	g->array = (int*)malloc(n*sizeof(int));
	check_allocate(g->array);
	g->size = n;
	for (i=0; i<n; i++){
		g->array[i] = i;
	}

	return g;
}

/* Calculating norma1 of matrix B.
 * Because B is a symmetric matrix summing over the columns is the same as over the rows.
 * Therefore we calculte the norma according the B rows.
 * norma1 = maximum i of sum(j=1,..,size of group g) abs(B[g]ij) */
double norma_calc(matB *B, group *g){
		int i, j, ng;
		double sum, max, *B_row_res;

		ng = g->size;
		max = 0.0;

		for (i = 0; i < ng; i++){
			/*getting row i of B_hat[g]*/
			B_row_res = get_B_row(B, g, 0, i);

			/*calculation sum abs*/
			sum = 0.0;
			for(j=0;j<ng;j++){
				sum += fabs(B_row_res[j]);
			}

			/*saving the maximum (norma1)*/
			if (sum>max){
				max = sum;
			}
		}
		return max;
}

/* Power iteration to calculate an estimation of the eigen vector */
double* power_iteration(matB *B, group *g, double norma){
	double dotProduct, *pointer;
	int i, diff, ng, iter_count=0;

	if (vector == NULL || nextVector == NULL){
		vector = (double*)malloc(B->Ag->n*sizeof(double));
		check_allocate(vector);
		nextVector = (double*)malloc(B->Ag->n*sizeof(double));
		check_allocate(nextVector);
	}

	ng = g->size;

	/* creating a random vector */
	srand(time(NULL));
	for (i = 0; i < ng; i++){
		vector[i] = rand();
	}

	/* power iteration */
	diff = 0;
	while (diff != ng){
			iter_count++;
			diff = 0;

			/* eigenVector = B*vector */
			mult_B_double(B, g, norma, vector, nextVector);

			/*calculating dot product*/
			dotProduct = (double) 0.0;
			for (i = 0; i < ng; i++){
				dotProduct += nextVector[i]*nextVector[i];
			}
			dotProduct = sqrt(dotProduct);

			if (dotProduct==0){
				printf("zero division problem - vector is all zeros");
				exit(1);
			}

			/*Normalizing eigenVector*/
			for (i = 0; i < ng; i++){
				nextVector[i] = nextVector[i]/(dotProduct);
				if (fabs(nextVector[i]-vector[i])<epsilon){
					diff++;
				}
			}

			/*switching vector and next vector pointers*/
			pointer = vector;
			vector = nextVector;
			nextVector = pointer;

			if (iter_count >= 10000*(ng+1)){
				printf("too many iteration in power iteration");
				exit(1);
			}
		}

	return vector;
}

/* Finding eigen pair.
 * The function returns the eigen vector calculated with power iteration.
 * eigen value = (eigen vector*B*eigen vector)/((eigen vector)*(eigen vector)) */
double* finding_eigen_pair(matB *B, group* g, double norma){
	int i, ng;
	double *eigenVector , eigenValue, numerator = 0, denominator = 0;

	ng = g->size;

	/*finding eigen vector*/
	eigenVector = power_iteration(B, g, norma);

	/*finding eigen value*/
	if (eigenVector_Mult_B == NULL){
		eigenVector_Mult_B = (double*)malloc(B->Ag->n*sizeof(double));
		check_allocate(eigenVector_Mult_B);
	}
	mult_B_double(B, g, norma, eigenVector,eigenVector_Mult_B);

	for (i=0;i<ng;i++){
		/* numerator = eigen vector*B*eigen vector */
		numerator += eigenVector[i]*eigenVector_Mult_B[i];
		/* denominator = (eigen vector)*(eigen vector) */
		denominator += eigenVector[i]*eigenVector[i];
	}
	if (denominator == 0){
		printf("Zero devision problem (eigenVector = 0)");
		exit(1);
	}
	eigenValue = (numerator/denominator) - norma;

	/* if eigenValue is negative, the group is indivisible
	 * therefore we return eigenVector = 1 to symbolize it */
	if (eigenValue<=0){
		for (i = 0; i < ng; i++){
			eigenVector[i] = 1;
		}
	}
	return eigenVector;
}

/* Calculating s(transpose)Bs */
double calc_sBs(matB *B, group* g, int* s){
	double *mult_Bs, sBs_result;
	int i,ng;

	ng = g->size;
	mult_Bs = mult_B_int(B, g, 0, s);

	sBs_result = 0.0;
	for (i = 0; i < ng; i++){
		sBs_result += s[i]*mult_Bs[i];
	}

	return sBs_result;
}

/* Calculating s using eigenVector.
 * Calculating the multiplication of s(transpose)Bs.
 * If s(transpose)Bs is negative the network is not divisible,
 * therefore returning s = 1 (the division s is only for one group). */
int* calc_s(matB *B, group* g, double* eigenVector){
	int ng, i;
	double sBs_value;

	if (s == NULL){
		s = (int*)malloc(B->Ag->n*sizeof(int));
		check_allocate(s);
	}

	ng = g->size;

	/*dividing to the groups according to eigenVector*/
	for (i = 0; i < ng; i++){
		if (eigenVector[i]>0){
			s[i] = 1;
		}
		else {
			s[i] = -1;
		}
	}

	/*checking sBs value, and if it's negative returning s=1 */
	sBs_value = calc_sBs(B, g, s);
	if (sBs_value<=0){
		for (i=0;i<ng;i++){
			s[i] = 1;
		}
	}
	return s;
}

/* Calculating deltaQ.
 * row_index (k) indicates the index that we changed in s (k is the vertex that was moved).
 * deltaQ = 4s[k]*sigma_j(B[g]kj*s[j])-4*B[g]kk */
double deltaQ(matB *B, group *g, int *s, int row_index){
	int i, ng;
	double sum = 0.0, *B_row_res;

	ng = g->size;
	B_row_res = get_B_row(B, g, 0, row_index);

	for(i = 0; i < ng; i++){
		sum += B_row_res[i]*s[i];
	}
	sum = 4*s[row_index]*sum -4*B_row_res[row_index];
	return sum;
}

/* Allocating unmoved, indices, improve and score arrays */
void allocate_improving_s(int n, int **unmoved, int **indices, double **improve, double **score){
	if (*unmoved == NULL || *improve == NULL || *indices == NULL || *score==NULL){
		*unmoved = (int*)malloc(n*sizeof(int));
		check_allocate(*unmoved);
		*indices = (int*)malloc(n*sizeof(int));
		check_allocate(indices);
		*improve = (double*)malloc(n*sizeof(double));
		check_allocate(*improve);
		*score = (double*)malloc(n*sizeof(double));
		check_allocate(*score);
		}
}

/* Improving the division s using modularity maximization - Algorithm 4
 * We used some math calculations to improve the multiplication s(t)Bs (deltaQ). */
void improving_s(matB *B, group *g, int* s){
	int ng, unmoved_size, i, j, k, iter=0;
	int index_max_score=0, max_index, max_index_improvement=0, prev_max_index = 0;
	double delQ, max_score=0, max_improvement=0, *B_row_res=NULL;

	allocate_improving_s(B->Ag->n, &unmoved, &indices, &improve, &score);
	ng = g->size;

	do{
		iter++;
		/*calculating unmoved - originally an array from 0,..,ng-1*/
		for (i = 0; i < ng; i++){
				unmoved[i] = i;
			}
		unmoved_size = ng;

		for (i = 0; i < ng; i++){
			if (i!=0){
				B_row_res = get_B_row(B, g, 0, prev_max_index);
			}

			/* j is the index of k (the item) in unmoved */
			for (j = 0; j < unmoved_size; j++){
				k = unmoved[j];
				s[k] = -s[k];

			    /* computing deltaQ for the move of each unmoved vertex */
				if (i!=0){
				/* we use the previous score[k] calculated in iteration i-1
				 * in order to calculate current score[k] more efficiently */
					score[k] = score[k] + 8*s[k]*s[prev_max_index]*B_row_res[k];
				}

				else {
					score[k] = deltaQ(B, g, s, k);
				}

				s[k] = -s[k];
				if (j==0 || score[k]>max_score){
					index_max_score = j;
					max_score = score[k];
				}
			}

			/* moving the vertex with max score to the other group
			 * (s[max_index]=-s[max_index]) */
			max_index = unmoved[index_max_score];
			s[max_index] = -s[max_index];
			indices[i] = max_index;

			if (i==0){
				improve[i] = max_score;
			}
			else {
				improve[i] = improve[i-1] + max_score;
			}

			if (i==0 || improve[i] > max_improvement){
				max_improvement = improve[i];
				max_index_improvement = i;
			}

			prev_max_index = unmoved[index_max_score];

		   /* in order to delete an item from the array unmoved
			* we swap it with the last item in the array,
			* and update unmoved_size */
			unmoved[index_max_score] = unmoved[unmoved_size-1];
			unmoved_size--;
		}

		/* finding the maximum improvement of s and updating s accordingly */
		for (i = ng-1; i >= max_index_improvement+1; i--){
			j = indices[i];
			s[j] = -s[j];
		}

		if (max_index_improvement == ng-1){
			delQ = 0;
		}
		else {
			delQ = improve[max_index_improvement];
		}

		if (iter>10000*(ng+1)){
			printf("too many iteration in modularity maximization");
			exit(1);
		}
	}
	while (delQ>epsilon);
}

/* Dividing g into two groups according to s (after improvement).
 * Adding the new group to P or O according to it's size  */
int divide_g(group** P, group **O, group* g, int *s, int *index_O, int index_P){
	int i, ng, countNeg, countPos, pointPos, pointNeg;
	group *groupPos, *groupNeg;

	/*dividing to two groups according to s*/
	ng = g->size;
	countPos = 0;
	countNeg = 0;
	for (i=0;i<ng;i++){
		if (s[i]<=0){
			countNeg++;
		}
		else {
			countPos++;
		}
	}

	/*if the division is just one group, the group can be moved to O*/
	if (countNeg==0 || countPos==0){
		O[*index_O] = g;
		(*index_O)++;
		return index_P;
	}

	/* building 2 groups - groupPos, groupNeg */
	groupPos = (group*)malloc(sizeof(group));
	check_allocate(groupPos);
	groupPos->size = countPos;
	groupPos->array = (int*)malloc(countPos*sizeof(int));
	check_allocate(groupPos->array);

	groupNeg = g; /*using g's memory in order to save allocation*/
	groupNeg->size = countNeg;

	/*using two pointers pointPos, and pointNeg*/
	pointPos = 0; pointNeg = 0;

	/*creating the new group's arrays*/
	for (i = 0; i < ng ; i++){
		if (s[i]>0){
			groupPos->array[pointPos++] = g->array[i];
		}
		else {
			groupNeg->array[pointNeg++] = g->array[i];
		}
	}

	/* adding group positive and group negative to O or P according to it size */
	if (countPos == 1){
		O[(*index_O)++] = groupPos;
	}
	else {
		P[index_P++] = groupPos;
	}

	if (countNeg == 1){
		O[(*index_O)++] = groupNeg;
	}
	else {
		P[index_P++] = groupNeg;
	}

	return index_P;
}

/* Writing the division of the graph to file out.
 * First char in output file is the number of groups,
 * followed by number of vertices in the group and then the vertices.
 * the vertices will be sorted in increasing order */
void write_to_output_file(group** O, FILE* fileOut, int index_O){
	int i, writeCharNum, sizeOfg;
	group *g;

	/*writing to fileOut the number of groups*/
	writeCharNum = fwrite(&index_O, sizeof(int), 1, fileOut);
	check_char_num(writeCharNum, 1, "problem with writing to output file");

	for (i=0;i<index_O;i++){
		g = O[i];
		sizeOfg = g->size;

		/*writing to fileOut the number of vertices in group g*/
		writeCharNum = fwrite(&sizeOfg, sizeof(int), 1, fileOut);
		check_char_num(writeCharNum, 1, "problem with writing to output file");

		/*writing to fileOut the vertices in group g*/
		writeCharNum = fwrite(g->array, sizeof(int), sizeOfg, fileOut);
		check_char_num(writeCharNum, sizeOfg, "problem with writing to output file");

		free(g->array);
		free(g);
	}
	fclose(fileOut);
}

/* Frees all dynamic allocations and memory resources */
void free_all(group** O, group** P, matB* B, matA_sparse *A){
	/*freeing cluster's resources*/
	free(O);
	free(P);
	free_matB(B);

	/*freeing the original mat A */
	free_spmat(A);
	free_global_spmat_arrays();

	/*freeing the global arrays used in Utility*/
	free(vector);
	free(nextVector);
	free(eigenVector_Mult_B);
	free(s);
	free(unmoved);
	free(indices);
	free(improve);
	free(score);
}



