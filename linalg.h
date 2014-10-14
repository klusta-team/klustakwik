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

integer Cholesky(SafeArray<scalar> &In, SafeArray<scalar> &Out, integer D);
integer MaskedCholesky(SafeArray<scalar> &In, SafeArray<scalar> &Out, integer D, vector<integer> &Masked, vector<integer> &Unmasked);
void TriSolve(SafeArray<scalar> &M, SafeArray<scalar> &x, SafeArray<scalar> &Out, integer D);
void MaskedTriSolve(SafeArray<scalar> &M, SafeArray<scalar> &x,
	SafeArray<scalar> &Out, integer D,
	vector<integer> &Masked, vector<integer> &Unmasked);

#endif /* LINALG_H_ */
