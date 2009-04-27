/*
 * Gausssian_parallel_all.h
 *
 *  Created on: Apr 27, 2009
 *      Author: yzhang8
 */

#ifndef GAUSSSIAN_PARALLEL_ALL_H_
#define GAUSSSIAN_PARALLEL_ALL_H_

int col;
int *max_rows;

void* pick_row(void* s);
void* compute_row(void* s);
void gaussian_elimination_all_parallel();
void gaussian_elimination_parallel();

#endif /* GAUSSSIAN_PARALLEL_ALL_H_ */
