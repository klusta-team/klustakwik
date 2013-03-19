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

using namespace std;

// Cholesky Decomposition
// In provides upper triangle of input matrix (In[i*D + j] >0 if j>=i);
// which is the top half of a symmetric matrix
// Out provides lower triange of output matrix (Out[i*D + j] >0 if j<=i);
// such that Out' * Out = In.
// D is number of dimensions
//
// returns 0 if OK, returns 1 if matrix is not positive definite
int Cholesky(SafeArray<scalar> &In, SafeArray<scalar> &Out, int D)
{
	int i, j, k;
	scalar sum;

	// empty output array
	for (i=0; i<D*D; i++) Out[i] = 0;

	// main bit
	for (i=0; i<D; i++) {
		for (j=i; j<D; j++) {	// j>=i
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

// Solve a set of linear equations M*Out = x.
// Where M is lower triangular (M[i*D + j] >0 if j>=i);
// D is number of dimensions
void TriSolve(SafeArray<scalar> &M, SafeArray<scalar> &x,
				SafeArray<scalar> &Out, int D)
{
	for(int i=0; i<D; i++)
	{
		scalar *MiD = &M[i*D];
		scalar sum = x[i];
		for(int j=0; j<i; j++) // j<i
			//sum += M[i*D + j] * Out[j];
			sum += MiD[j] * Out[j];
		//Out[i] = - sum / M[i*D + i];
		Out[i] = - sum / MiD[i];
	}
}
