/*
 * linalg.h
 *
 * Linear algebra
 *
 *  Created on: 13 Nov 2011
 *      Author: dan
 */

#ifndef LINALG_H_
#define LINALG_H_

#include "numerics.h"

int Cholesky(SafeArray<scalar> &In, SafeArray<scalar> &Out, int D);
void TriSolve(SafeArray<scalar> &M, SafeArray<scalar> &x, SafeArray<scalar> &Out, int D);

#endif /* LINALG_H_ */
