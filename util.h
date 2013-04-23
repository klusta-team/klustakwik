 /*
 * util.h
 *
 *  Created on: 11 Nov 2011
 *      Author: dan
 */

#ifndef UTIL_H_
#define UTIL_H_

#include<stdio.h>
#include "numerics.h"
#include "log.h"

extern const scalar HugeScore;

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

int irand(int min, int max);
FILE *fopen_safe(char *fname, char *mode);
void MatPrint(FILE *fp, scalar *Mat, int nRows, int nCols);

#endif /* UTIL_H_ */
