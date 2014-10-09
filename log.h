/*
 * log.h
 *
 * Used for logging output to stdout and to a file.
 *
 *  Created on: 11 Nov 2011
 *      Author: dan
 */

#ifndef LOG_H_
#define LOG_H_

#include <stdio.h>

extern FILE *logfp;

void Output(char *fmt, ...);
void Error(char *fmt, ...);
void FlushLog();

#endif /* LOG_H_ */
