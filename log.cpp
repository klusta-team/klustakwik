/*
 * log.cpp
 *
 * Used for logging output to stdout and to a file.
 *
 *  Created on: 11 Nov 2011
 *      Author: dan
 */

#include "log.h"
#include "parameters.h"
#include<stdio.h>
#include <stdarg.h>

FILE *logfp;

// Write to screen and log file
void Output(char *fmt, ...) {
    va_list arg;
    char str[STRLEN];

    if (!Screen && !Log) return;
    va_start(arg, fmt);
    vsnprintf(str,STRLEN,fmt,arg);
    va_end(arg);

    if (Screen) printf("%s", str);
    if (Log) fprintf(logfp, "%s", str);
}

// Print an error message and abort
void Error(char *fmt, ...) {
    va_list arg;
    char str[STRLEN];

    if (!Screen && !Log) return;
    va_start(arg, fmt);
    vsnprintf(str,STRLEN,fmt,arg);
    va_end(arg);

    if (Screen) fprintf(stderr, "%s", str);
    if (Log) fprintf(logfp, "%s", str);
}
