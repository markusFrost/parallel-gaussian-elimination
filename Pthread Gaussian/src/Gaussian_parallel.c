/*
 * Gaussian_parallel_all.c
 *
 *  Created on: Apr 27, 2009
 *      Author: yzhang8
 */
#include "Gaussian.h"
#include "Gaussian_parallel.h"

void gaussian_elimination_all_parallel() {
	int i, j;
	time_t time_sec_start, time_sec_finish;
	time_sec_start = time(NULL);
	pthread_t* p_threads = malloc(thread_num * sizeof(pthread_t));
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	max_rows = malloc(thread_num * sizeof(int));
	thread_row_range* thread_arg = malloc(thread_num * sizeof(thread_row_range));
	int useful_thread_num;//the real number of thread used
	int thread_first_row, thread_row_size;
	double *row_swap = malloc(size * sizeof(double)), b_swap;
	char buffer[40];
	for (col = 0; col < size - 1; col++) {// go through each column
		useful_thread_num = (size - col - 1) / thread_row_threshold + 1;//calculate how many threads should be used
		if (useful_thread_num > thread_num) {
			useful_thread_num = thread_num;
		}
		thread_first_row = col + 1;
		thread_row_size = (size - col - 1) / useful_thread_num;
		for (i = 0; i < useful_thread_num; i++) {
			thread_arg[i].first_row = thread_first_row;
			thread_arg[i].thread_id = i;
			if (i != useful_thread_num - 1) {//not the last thread
				thread_arg[i].last_row = thread_first_row + thread_row_size - 1;
				thread_first_row = thread_arg[i].last_row + 1;
			} else {
				thread_arg[i].last_row = size - 1;
			}
		}
		for (i = 0; i < useful_thread_num; i++) {//run threads
			pthread_create(&p_threads[i], &attr, pick_row, (void*) &thread_arg[i]);
		}
		for (i = 0; i < useful_thread_num; i++) {//wait for all thread to complete finding the largest pivot
			pthread_join(p_threads[i], NULL);
		}
		int pivot_row = max_rows[0];
		for (i = 1; i < useful_thread_num; i++) {//find the largest
			if (matrix_A[pivot_row][col] < matrix_A[max_rows[i]][col]) {
				pivot_row = max_rows[i];
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
		for (i = 0; i < useful_thread_num; i++) {//run threads
			pthread_create(&p_threads[i], &attr, compute_row, (void*) &thread_arg[i]);
		}
		for (i = 0; i < useful_thread_num; i++) {//wait for all thread to complete calculate the rows
			pthread_join(p_threads[i], NULL);
		}
		sprintf(buffer, "\nGaussian after proccess col %d\n ", col);
		print_Gaussian(buffer);
	}
	time_sec_finish = time(NULL);
	printf("\nGaussian elimination used %d sec\n", (time_sec_finish - time_sec_start));
	free(p_threads);
	free(thread_arg);
}

void gaussian_elimination_parallel() {
	print_Gaussian("before Gaussian!");
	time_t time_sec_start, time_sec_finish;
	time_sec_start = time(NULL);
	pthread_t* p_threads = malloc(thread_num * sizeof(pthread_t));
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	thread_row_range* thread_arg = malloc(thread_num * sizeof(thread_row_range));
	int i,j, thread_first_row, thread_row_size,useful_thread_num;
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

		useful_thread_num = (size - col - 1) / thread_row_threshold + 1;//calculate how many threads should be used
		if (useful_thread_num > thread_num) {
			useful_thread_num = thread_num;
		}
		thread_first_row = col + 1;
		thread_row_size = (size - col - 1) / useful_thread_num;
		for (i = 0; i < useful_thread_num; i++) {
			thread_arg[i].first_row = thread_first_row;
			thread_arg[i].thread_id = i;
			if (i != useful_thread_num - 1) {//not the last thread
				thread_arg[i].last_row = thread_first_row + thread_row_size - 1;
				thread_first_row = thread_arg[i].last_row + 1;
			} else {
				thread_arg[i].last_row = size - 1;
			}
		}

		for (i = 0; i < useful_thread_num; i++) {//run threads
			pthread_create(&p_threads[i], &attr, compute_row, (void*) &thread_arg[i]);
		}
		for (i = 0; i < useful_thread_num; i++) {//wait for all thread to complete calculate the rows
			pthread_join(p_threads[i], NULL);
		}
		sprintf(buffer, "\nGaussian after proccess col %d\n ", col);
		print_Gaussian(buffer);
	}
	time_sec_finish = time(NULL);
	printf("\nGaussian elimination used %d sec\n", (time_sec_finish - time_sec_start));
	free(p_threads);
	free(thread_arg);
}

void* compute_row(void* s) {
	thread_row_range* range = (thread_row_range*) s;
	int pivot_row = col;
	int row, j;
	double k;
	for (row = range->first_row; row <= range->last_row; row++) {
		k = matrix_A[row][col] / matrix_A[pivot_row][col];
		for (j = col; j < size; j++) {
			matrix_A[row][j] -= matrix_A[pivot_row][j] * k;
		}
		vector_B[row] -= vector_B[pivot_row] * k;
	}
	pthread_exit(0);
}

void* pick_row(void* s) {
	thread_row_range* range = (thread_row_range*) s;
	int pivot_row = col;
	int row;
	for (row = range->first_row; row <= range->last_row; row++) {
		if (matrix_A[pivot_row][col] < matrix_A[row][col]) {
			pivot_row = row;
		}
	}
	max_rows[range->thread_id] = pivot_row;
	pthread_exit(0);
}

