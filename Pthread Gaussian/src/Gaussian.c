/*
 ============================================================================
 Name        : Gaussian.c
 Author      : Ryan Zhang
 Version     :
 Copyright   : It's patent
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>
#include <math.h>

#define MATRIX_ENTRY_VALUE_MAX 1000;

struct range {
	int first_row;
	int last_row;
	int thread_id;
};
typedef struct range thread_row_range;//rows for each thread;

double** matrix_A;
double* vector_B;
double *vector_x;
int col;//for each thread
int* max_rows;//for pivoting
int size;

inline void print_Gaussian();
void* pick_row(void*);
void* compute_row(void*);

int main(int argc, char* argv[]) {
	puts("!!!Start Gaussian Elimination!!!");
	time_t time_sec_start, time_sec_finish;
	int i, j;
	size = atoi(argv[1]);
	int thread_num = atoi(argv[2]);//the maximum number of thread possible
	int thread_row_threshold = atoi(argv[3]);//the minimum number of rows a thread should be assigned to
	int useful_thread_num;//the real number of thread used
	thread_row_range* thread_arg = malloc(thread_num * sizeof(thread_row_range));
	max_rows = malloc(thread_num * sizeof(int));
	pthread_t* p_threads = malloc(thread_num * sizeof(pthread_t));
	pthread_attr_t attr;
	pthread_attr_init(&attr);
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
		//vector_B[i] = rand() % MATRIX_ENTRY_VALUE_MAX;
		test_Vector[i] = vector_B[i];
		for (j = 0; j < size; j++) {
			matrix_A[i][j] = ((double) rand()) / RAND_MAX * MATRIX_ENTRY_VALUE_MAX;
			//matrix_A[i][j] = rand() % MATRIX_ENTRY_VALUE_MAX ;
			test_Matrix[i][j] = matrix_A[i][j];
		}
	}
	print_Gaussian("before Gaussian!");
	time_sec_start = time(NULL);
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
	free(thread_arg);
	free(p_threads);
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
	free(max_rows);
	return EXIT_SUCCESS;
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
