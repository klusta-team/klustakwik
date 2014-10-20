/*
 * util.cpp
 *
 * Various utility functions.
 *
 *  Created on: 11 Nov 2011
 *      Author: dan
 */

// Disable some Visual Studio warnings
#define _CRT_SECURE_NO_WARNINGS

#include "util.h"
#include<stdlib.h>
#include<stdarg.h>
#include "numerics.h"

const scalar HugeScore = (scalar)1e32;

/* integer random number between min and max*/
integer irand(integer min, integer max)
{
    return (rand() % (max - min + 1) + min);
}

FILE *fopen_safe(char *fname, char *mode) {
    FILE *fp;

    fp = fopen(fname, mode);
    if (!fp) {
        fprintf(stderr, "Could not open file %s\n", fname);
        abort();
    }

    return fp;
}

// Print a matrix
void MatPrint(FILE *fp, scalar *Mat, integer nRows, integer nCols) {
    integer i, j;

    for (i=0; i<nRows; i++) {
        for (j=0; j<nCols; j++) {
            fprintf(fp, "%.5g ", Mat[i*nCols + j]);
            Output("%.5g ", Mat[i*nCols + j]);
        }
        fprintf(fp, "\n");
        Output("\n");
    }
}

void CompareVectors(scalar *A, scalar *B, integer N)
{
	scalar meanerr = 0.0;
	scalar maxerr = 0.0;
	integer nerr = 0;
	integer ntotal = 0;
	for (integer i = 0; i < N; i++)
	{
		scalar err = fabs(A[i] - B[i]);
		ntotal++;
		if (err > 0)
		{
			nerr++;
			meanerr += err;
			if (err > maxerr)
				maxerr = err;
		}
	}
	if (nerr)
	{
		meanerr /= nerr;
		cout << "Comparison error n=" << nerr << " (" << (100.0*nerr) / ntotal << "%), mean=" << meanerr << ", max=" << maxerr << endl;
	}
	else
		cout << "No comparison error." << endl;
}