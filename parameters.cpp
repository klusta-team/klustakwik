/*
 * parameters.cpp
 *
 * Definition and reading of the parameters. Include "parameters.h" to
 * use these values.
 *
 *  Created on: 11 Nov 2011
 *      Author: dan
 */

// Disable some Visual Studio warnings
#define _CRT_SECURE_NO_WARNINGS

#include "parameters.h"
#include "math.h"
#include "log.h"
#include "util.h"
#include <stdlib.h>
#include <string.h>
#include "numerics.h"

// See parameters.h for an explanation
PARAMETERS_TABLE(DEFINE_PARAMETERS)

char HelpString[] = "\n\
\n\
Masked KlustaKwik\n\
\n\
Uses the CEM algorithm with masks to do automatic clustering.\n\n\
";

integer num_used_arguments = 0;

void params_error()
{
    fprintf(stderr, "Usage: MaskedKlustaKwik FileBase ElecNo [Arguments]\n\n");
    fprintf(stderr, "Arguments (with default values): \n\n");
    print_params(stderr);
    exit(1);
}

void SetupParams(integer argc, char **argv) {
    char fname[STRLEN];

    init_params(argc, argv);

    // PARAMETER DEFINITIONS GO HERE
    // This line reads the PARAMETERS_TABLE from parameters.h and reads each
    // of the named parameters from the command line file.
    PARAMETERS_TABLE(INPUT_PARAMETERS)

    if (argc<3)
    {
        fprintf(stderr, "Not enough command line arguments.\n\n");
        params_error();
    }

    strcpy(FileBase, argv[1]);
    ElecNo = atoi(argv[2]);

    if(2*num_used_arguments!=argc-3)
    {
        fprintf(stderr, "%d %d\n", num_used_arguments, argc);
        fprintf(stderr, "Unrecognised command line arguments.\n\n");
        params_error();
    }

    //if(UseDistributional)
    //{
    //    UseMaskedMStep = false;
    //    UseMaskedEStep = false;
    //}

    if (Screen) print_params(stdout);

    // open log file, if required
    if (Log) {
        sprintf(fname, "%s.klg.%d", FileBase, ElecNo);
        logfp = fopen_safe(fname, "w");
        print_params(logfp);
    }
}

typedef struct entry_t {
    integer t;
    char *name;
    void *addr;
    struct entry_t *next;
} entry;


entry *top, *bottom;
integer argc;
char **argv;
extern char HelpString[];

char help = 0;

/* returns 1 if the parameter was found and changed, else zero. */
integer change_param(char *name, char *value)
{
    entry *e;
    integer changed = 0;

    for(e=bottom; e; e = e->next) if (!strcmp(name, e->name)) {
        switch (e->t) {
        case FLOAT:
            *((scalar *) e->addr) = atof(value); break;
        case INT:
            *((integer *) e->addr) = atoi(value); break;
        case BOOLEAN:
            if (*value == '0')
                *((char *) e->addr) = 0;
            else
                *((char *) e->addr) = 1;
            break;
        case STRING:
            strncpy((char *)e->addr, value, STRLEN); break;
        }
        changed = 1;
        break;
    }
    num_used_arguments += changed;
    return changed;
}

void init_params(integer ac, char **av)
{
    argc = ac;
    argv = av;
}

void search_command_line(char *name)
{
    integer i;

    for(i=0; i<argc-1; i++)
        if (argv[i][0] == '-' && !strcmp(argv[i]+1, name))
            change_param(argv[i] + 1, argv[i+1]);
    if (argv[argc-1][0] == '-' && !strcmp(argv[argc-1]+1, name))
        change_param(argv[argc-1] + 1, "");
}

void add_param(integer t, char *name, void *addr)
{
    entry *e;
    if (top == NULL) {
        bottom = top = e = (entry *) malloc(sizeof(entry));
    } else {
        e = (entry *) malloc(sizeof(entry));
        top->next = e;
        top = e;
    }
    if (e == NULL) {printf("parameter manager out of memory!\n"); exit(1);}
    top->t = t;
    top->name = name;
    top->addr = addr;
    top->next = NULL;
    search_command_line(name);
}

void print_params(FILE *fp)
{
    entry *e;

    add_param(BOOLEAN, "help", &help);
    if (help)
    {
        fprintf(fp, HelpString);
        PARAMETERS_TABLE(PRINT_PARAMETERS);
        exit(0);
    }
    else
    {
        for(e=bottom; e; e = e->next) {
            fprintf(fp, "%s\t", e->name);
            switch (e->t) {
            case FLOAT:
                fprintf(fp, SCALARFMT "\n", *(scalar *)(e->addr)); break;
            case INT:
                fprintf(fp, "%d\n", *(integer *)(e->addr)); break;
            case BOOLEAN:
                fprintf(fp, "%d\n", *(char *)(e->addr)); break;
            case STRING:
                fprintf(fp, "%s\n", (char *)(e->addr)); break;
            }
        }
    }
}
