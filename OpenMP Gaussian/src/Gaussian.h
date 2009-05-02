/*
 * Gaussian.h
 *
 *  Created on: Apr 30, 2009
 *      Author: yzhang8
 */

#ifndef GAUSSIAN_H_
#define GAUSSIAN_H_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "omp.h"

double** matrix_A;
double* vector_B;
double *vector_x;
int size, thread_num, block_size;

inline void print_Gaussian();


#endif /* GAUSSIAN_H_ */
