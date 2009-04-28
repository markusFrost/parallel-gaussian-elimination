/*
 * Gaussian.c
 *
 *  Created on: Apr 27, 2009
 *      Author: yzhang8
 */

#include "Gaussian.h"
#include "Gaussian_parallel.h"
#include "Gaussian_Block_Parallel.h"

#define MATRIX_ENTRY_VALUE_MAX 1000;

int main(int argc, char* argv[]) {
	puts("!!!Start Gaussian Elimination!!!");
	if(argc !=4){
		puts("usage: Gaussian size thread_max thread_threshold");
		return EXIT_FAILURE;
	}
	int i, j;
	size = atoi(argv[1]);
	thread_num = atoi(argv[2]);//the maximum number of thread possible
	thread_row_threshold = atoi(argv[3]);//the minimum number of rows a thread should be assigned to
	matrix_A = (double**) malloc(size * sizeof(double));
	vector_B = (double*) malloc(size * sizeof(double));
	double ** test_Matrix = (double**) malloc(size * sizeof(double));
	double *test_Vector = (double*) malloc(size * sizeof(double));
	vector_x = (double*) malloc(size * sizeof(double));
	for (i = 0; i < size; i++) {
		matrix_A[i] = malloc(size * sizeof(double));
		test_Matrix[i] = malloc(size * sizeof(double));
	}
	srand(time(NULL));
	for (i = 0; i < size; i++) {
		vector_B[i] = ((double) rand()) / RAND_MAX * MATRIX_ENTRY_VALUE_MAX;
		test_Vector[i] = vector_B[i];
		for (j = 0; j < size; j++) {
			matrix_A[i][j] = ((double) rand()) / RAND_MAX * MATRIX_ENTRY_VALUE_MAX;
			test_Matrix[i][j] = matrix_A[i][j];
		}
	}
	print_Gaussian("before Gaussian!");
	////////////
	//gaussian_elimination_parallel();
	//gaussian_elimination_all_parallel();
	gaussian_elimination_block_parallel();
	////////////
	double sumax = 0;
	for (i = size - 1; i >= 0; i--) {//backward substitution
		sumax = 0;
		for (j = i + 1; j < size; j++) {
			sumax = sumax + matrix_A[i][j] * vector_x[j];
		}
		vector_x[i] = (vector_B[i] - sumax) / matrix_A[i][i];
	}
	print_Gaussian("after Gaussian!");
	puts("the result!\n----------------------------------------------");
	double* vector_result = (double*) malloc(size * sizeof(double));
	double l2 = 0, temp;
	for (i = 0; i < size; i++) {
		temp = 0;
		for (j = 0; j < size; j++) {
			temp += test_Matrix[i][j] * vector_x[j];
		}
		temp -= test_Vector[i];
		l2 += temp * temp;
	}
	/*for (i = 0; i < size; i++) {
	 for (j = 0; j < size; j++) {
	 printf("%f ", test_Matrix[i][j]);
	 }
	 printf("   %f", vector_x[i]);
	 printf("   %f", test_Vector[i]);
	 printf(" diff %f\n", temp);
	 }*/
	l2 = sqrt(l2);
	printf("-----------------------------------------------------\nthe residue value is %f\n", l2);
	free(vector_B);
	free(test_Vector);
	free(vector_result);
	free(vector_x);
	for (i = 0; i < size; i++) {
		free(matrix_A[i]);
		free(test_Matrix[i]);
	}
	free(matrix_A);
	free(test_Matrix);
	return EXIT_SUCCESS;
}

inline void print_Gaussian(const char* s) {
	int i, j;
	/*puts(s);
	 for (i = 0; i < size; i++) {
	 for (j = 0; j < size; j++) {
	 printf("%f ", matrix_A[i][j]);
	 }
	 printf("   %f", vector_x[i]);
	 printf("   %f\n", vector_B[i]);
	 }*/
}
