/*
 * Gaussian_Sequential_Block.c
 *
 *  Created on: Apr 30, 2009
 *      Author: yzhang8
 */
#ifndef GAUSSIAN_SEQUENTIAL_BLOCK_C_
#define GAUSSIAN_SEQUENTIAL_BLOCK_C_

#include "Gaussian.h"
#include "Gaussian_Sequential.h"

int volatile col;

inline int min(int, int);

void gaussian_sequential_block() {
	int i, j;
	time_t time_sec_start, time_sec_finish;
	time_sec_start = time(NULL);
	double *row_swap = malloc(size * sizeof(double)), b_swap;
	double *row_ratio = malloc(size * sizeof(double));
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
		int row, rr, jj, j, rsize, jsize;
		for (row = col+1; row <size; row++) {
			row_ratio[row] = matrix_A[row][col] / matrix_A[col][col];
		}
		for (row = col+1; row <size; row += block_size) {
			rsize = min(row + block_size, size);
			for (j = col; j < size; j += block_size) {
				jsize = min(j + block_size, size);
				for (rr = row; rr < rsize; rr++) {
					for (jj = j; jj < jsize; jj++) {
						matrix_A[rr][jj] -= matrix_A[col][jj] * row_ratio[rr];
					}
				}
			}
		}
		for (row = col+1; row <size; row++) {
			vector_B[row] -= vector_B[col] * row_ratio[row];
		}
		sprintf(buffer, "\nGaussian after proccess col %d\n ", col);
		print_Gaussian(buffer);
	}
	time_sec_finish = time(NULL);
	printf("\nGaussian elimination used %d sec\n", (int) (time_sec_finish - time_sec_start));
	free(row_swap);
	free(row_ratio);
}

inline int min(int a, int b) {
	return ((a > b) ? b : a);
}

#endif /* GAUSSIAN_SEQUENTIAL_BLOCK_C_ */
