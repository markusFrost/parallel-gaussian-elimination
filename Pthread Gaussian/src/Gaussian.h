/*
 * Gaussian.h
 *
 *  Created on: Apr 27, 2009
 *      Author: yzhang8
 */

#ifndef GAUSSIAN_H_
#define GAUSSIAN_H_

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>
#include <math.h>


struct range {
	int first_row;
	int last_row;
	int thread_id;
};
typedef struct range thread_row_range;//rows for each thread;

double** matrix_A;
double* vector_B;
double *vector_x;
int size, thread_num, thread_row_threshold;

inline void print_Gaussian();

#endif /* GAUSSIAN_H_ */
