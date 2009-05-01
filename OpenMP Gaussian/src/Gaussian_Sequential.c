/*
 * Gaussian_Sequential.c
 *
 *  Created on: Apr 30, 2009
 *      Author: yzhang8
 */

#ifndef GAUSSIAN_SEQUENTIAL_C_
#define GAUSSIAN_SEQUENTIAL_C_

#include "Gaussian.h"
#include "Gaussian_Sequential.h"

int volatile col;

void gaussian_sequential() {
	int i, j;
	time_t time_sec_start, time_sec_finish;
	time_sec_start = time(NULL);
	double *row_swap = malloc(size * sizeof(double)), b_swap;
	char buffer[40];
	for (col = 0; col < size - 1; col++) {// go through each column
		int pivot_row = col;
		for (i = col + 1; i < size; i++) {//find the largest
			if (matrix_A[pivot_row][col] < matrix_A[i][col]) {
				pivot_row = i;
			}
		}
		if (pivot_row != col) {//exchange rows
			b_swap = vector_B[col];
			vector_B[col] = vector_B[pivot_row];
			vector_B[pivot_row] = b_swap;
			for (j = col; j < size; j++) {
				row_swap[j] = matrix_A[col][j];
				matrix_A[col][j] = matrix_A[pivot_row][j];
				matrix_A[pivot_row][j] = row_swap[j];
			}
		}
		pivot_row = col;//reset to the top
		int row, j;
		double k;
		for (row = col + 1; row < size; row++) {
			k = matrix_A[row][col] / matrix_A[pivot_row][col];
			for (j = col; j < size; j++) {
				matrix_A[row][j] -= matrix_A[pivot_row][j] * k;
			}
			vector_B[row] -= vector_B[pivot_row] * k;
		}
		sprintf(buffer, "\nGaussian after proccess col %d\n ", col);
		print_Gaussian(buffer);
	}
	time_sec_finish = time(NULL);
	printf("\nGaussian elimination used %d sec\n", (int) (time_sec_finish - time_sec_start));
	free(row_swap);
}
#endif /* GAUSSIAN_SEQUENTIAL_C_ */
