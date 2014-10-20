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

// Matrix which has the form of a block matrix and a nonzero diagonal
// Unmasked gives the indices corresponding to the block matrix
// Masked gives the indices corresponding to the diagonal
class BlockPlusDiagonalMatrix
{
public:
	vector<scalar> Block;
	vector<scalar> Diagonal;
	vector<integer> &Unmasked;
	vector<integer> &Masked;
	integer NumUnmasked, NumMasked;
	BlockPlusDiagonalMatrix(vector<integer> &_Masked, vector<integer> &_Unmasked);
	void compare(scalar *Flat);
};

integer Cholesky(SafeArray<scalar> &In, SafeArray<scalar> &Out, integer D);
integer MaskedCholesky(SafeArray<scalar> &In, SafeArray<scalar> &Out, integer D, vector<integer> &Masked, vector<integer> &Unmasked);
integer BPDCholesky(BlockPlusDiagonalMatrix &In, BlockPlusDiagonalMatrix &Out);
void TriSolve(SafeArray<scalar> &M, SafeArray<scalar> &x, SafeArray<scalar> &Out, integer D);
void MaskedTriSolve(SafeArray<scalar> &M, SafeArray<scalar> &x,
	SafeArray<scalar> &Out, integer D,
	vector<integer> &Masked, vector<integer> &Unmasked);

#endif /* LINALG_H_ */
