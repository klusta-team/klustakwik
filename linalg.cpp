/*
 * linalg.cpp
 *
 * Linear algebra routines
 *
 *  Created on: 13 Nov 2011
 *      Author: dan
 */

#include "linalg.h"
#include<math.h>
#include<vector>
#include "numerics.h"

using namespace std;

// Cholesky Decomposition
// In provides upper triangle of input matrix (In[i*D + j] >0 if j>=i);
// which is the top half of a symmetric matrix
// Out provides lower triange of output matrix (Out[i*D + j] >0 if j<=i);
// such that Out' * Out = In.
// D is number of dimensions
//
// returns 0 if OK, returns 1 if matrix is not positive definite
integer Cholesky(SafeArray<scalar> &In, SafeArray<scalar> &Out, integer D)
{
    integer i, j, k;
    scalar sum;

    // empty output array
    for (i=0; i<D*D; i++) Out[i] = 0;

    // main bit
    for (i=0; i<D; i++) {
        for (j=i; j<D; j++) {    // j>=i
            sum = In[i*D + j];

            for (k=i-1; k>=0; k--) sum -= Out[i*D + k] * Out[j*D + k]; // i,j >= k
            if (i==j) {
                if (sum <=0) return(1); // Cholesky decomposition has failed
                Out[i*D + i] = (scalar)sqrt(sum);
            }
            else {
                Out[j*D + i] = sum/Out[i*D + i];
            }
        }
    }

    return 0; // for sucess
}

integer MaskedCholesky(SafeArray<scalar> &In, SafeArray<scalar> &Out, integer D, vector<integer> &Masked, vector<integer> &Unmasked)
{
	integer i, j, k;
	integer ii, jj, kk;
	scalar sum;

	integer NumUnmasked = (integer)Unmasked.size();

	// empty output array
	for (i = 0; i<D*D; i++) Out[i] = 0;

	// main bit for unmasked features
	//for (i = 0; i<D; i++)
	for (ii = 0; ii < NumUnmasked; ii++)
	{
		i = Unmasked[ii];
		//for (j = i; j<D; j++) 
		for (jj = ii; jj < NumUnmasked; jj++)
		{    // j>=i
			j = Unmasked[jj];
			sum = In[i*D + j];

			//for (k = i - 1; k >= 0; k--)
			for (kk = ii - 1; kk >= 0; kk--)
			{
				k = Unmasked[kk];
				sum -= Out[i*D + k] * Out[j*D + k]; // i,j >= k
			}
			if (i == j) {
				if (sum <= 0) return(1); // Cholesky decomposition has failed
				Out[i*D + i] = (scalar)sqrt(sum);
			}
			else {
				Out[j*D + i] = sum / Out[i*D + i];
			}
		}
	}
	// main bit for masked features
	for (ii = 0; ii < (integer)Masked.size(); ii++)
	{
		i = Masked[ii];
		scalar sum = In[i*D + i];
		if (sum <= 0)
			return 1; // Cholesky failed
		Out[i*D + i] = (scalar)sqrt(sum);
	}

	return 0; // for sucess
}


// Solve a set of linear equations M*Out = x.
// Where M is lower triangular (M[i*D + j] >0 if j>=i);
// D is number of dimensions
void TriSolve(SafeArray<scalar> &M, SafeArray<scalar> &x,
                SafeArray<scalar> &Out, integer D)
{
    for(integer i=0; i<D; i++)
    {
        scalar *MiD = &M[i*D];
        scalar sum = x[i];
        for(integer j=0; j<i; j++) // j<i
            //sum += M[i*D + j] * Out[j];
            sum += MiD[j] * Out[j];
        //Out[i] = - sum / M[i*D + i];
        Out[i] = - sum / MiD[i];
    }
}

//void MaskedTriSolve(SafeArray<scalar> &M, SafeArray<scalar> &x,
//                	SafeArray<scalar> &Out, integer D,
//					vector<integer> &Masked, vector<integer> &Unmasked)
//{
//	integer NumUnmasked = (integer)Unmasked.size();
//	integer NumMasked = (integer)Masked.size();
//	for (integer ii = 0; ii < NumUnmasked; ii++)
//	{
//		integer i = Unmasked[ii];
//		scalar sum = x[i];
//		for (integer jj = 0; jj < ii; jj++) // j<i
//		{
//			integer j = Unmasked[jj];
//			sum += M[i*D + j] * Out[j];
//		}
//		Out[i] = - sum / M[i*D + i];
//	}
//	for (integer ii = 0; ii < NumMasked; ii++)
//	{
//		integer i = Masked[ii];
//		Out[i] = -x[i] / M[i*D + i];
//	}
//}

// fast version with pointers and restricts
void FastMaskedTriSolve(scalar * __restrict M, scalar * __restrict x,
	scalar * __restrict Out, integer D,
	integer * __restrict Masked, integer * __restrict Unmasked,
	integer NumMasked, integer NumUnmasked)
{
	for (integer ii = 0; ii < NumUnmasked; ii++)
	{
		const integer i = Unmasked[ii];
		scalar sum = x[i];
		scalar * __restrict MiD = M + i*D;
		for (integer jj = 0; jj < ii; jj++) // j<i
		{
			const integer j = Unmasked[jj];
			sum += MiD[j] * Out[j];
		}
		Out[i] = -sum / MiD[i];
	}
	for (integer ii = 0; ii < NumMasked; ii++)
	{
		const integer i = Masked[ii];
		Out[i] = -x[i] / M[i*D + i];
	}
}

void MaskedTriSolve(SafeArray<scalar> &M, SafeArray<scalar> &x,
	SafeArray<scalar> &Out, integer D,
	vector<integer> &Masked, vector<integer> &Unmasked)
{
	FastMaskedTriSolve(&(M[0]), &(x[0]), &(Out[0]), D, &(Masked[0]), &(Unmasked[0]), Masked.size(), Unmasked.size());
}
