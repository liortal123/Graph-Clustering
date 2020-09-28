/**
 * Utility.h
 *
 * DESCRIPTION:
 * This header holds the primary methods of the graph clustering algorithms
 * In addition it holds the implementation of struct group that holds a group of vertices.
 *
 * check_allocate		- Checks if allocation was successful (else ends the program)
 * open					- Opens input file reads size of the adjacency matrix and calculates nnz
 * creating_g			- Creates group g
 * norma_calc			- Calculates norma1 of matrix B
 * finding_eigen_pair	- Finds leading eigenpair
 * calc_s				- Calculates division s
 * improving_s			- Improves the division s
 * divide_g				- Divides group g to two groups according to the division s.
 * write_to_output_file - Writes results to output file
 * free_all				- Frees all dynamic allocations and memory resources
 *
 */

#ifndef UTILITY_H_
#define UTILITY_H_

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <time.h>
#include <math.h>

#include "matA_sparse.h"
#include "matB.h"

typedef struct group
{
	int  size;
	int *array;
} group;

/* checking the allocation is successful else ending the program */
void check_allocate(void *pointer);

/* Opening and reading the size of the adjacency matrix and calculating nnz
 * Input file: the first value is n=|V|, the function returns it */
int open(char *argv[], FILE **fileIn, FILE **fileOut, int *nnz);

/* creating the initializing group g with all the vertices of the graph */
group* creating_g(int n);

/* Calculating norma1 of matrix B.
 * Because B is a symmetric matrix summing over the columns is the same as over the rows.
 * Therefore we calculte the norma according the B rows.
 * norma1 = maximum i of sum(j=1,..,size of group g) abs(B[g]ij) */
double norma_calc(matB *B, group *g);

/* Finding eigen pair.
 * The function returns the eigen vector calculated with power iteration.
 * eigen value = (eigen vector*B*eigen vector)/((eigen vector)*(eigen vector)) */
double* finding_eigen_pair(matB *B, group* g, double norma);

/* Calculating s using eigenVector.
 * Calculating the multiplication of s(transpose)Bs.
 * If s(transpose)Bs is negative the network is not divisible,
 * therefore returning s = 1 (the division s is only for one group). */
int* calc_s(matB *B, group* g, double* eigenVector);

/* Improving the division s using modularity maximization - Algorithm 4
 * We used some math calculations to improve the multiplication s(t)Bs (deltaQ). */
void improving_s(matB *B, group *g, int* s);

/* Dividing g into two groups according to s (after improvement).
 * Adding the new group to P or O according to it's size  */
int divide_g(group** P, group **O, group* g, int *s, int* index_O, int index_P);

/* Writing the division of the graph to file out.
 * First char in output file is the number of groups,
 * followed by number of vertices in the group and then the vertices.
 * the vertices will be sorted in increasing order */
void write_to_output_file(group** O, FILE* fileOut,  int index_O);

/* Frees all dynamic allocations and memory resources */
void free_all(group** O, group** P, matB* B, matA_sparse *A);

#endif /* UTILITY_H_ */
