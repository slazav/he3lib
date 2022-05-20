#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "he3tab.h"

/* This is a command line interface for the he3lib.
                                       slazav, 2014 */

/* Print short error message and exit */
void
error(const char * err){
  printf("Error: %s\n", err);
  exit(1);
}

/* Parse range setting in form <v1>:<step>:<v2>
   Return number of steps. */
int
parse_range(const char * str,
      double * v1, double * v2, double *st){
  const char *del=":";
  char *tok;

  tok = strtok((char*)str,  del);
  if (tok) *v1 = atof(tok);
  else error("bad parameter");

  tok = strtok(NULL, del);
  if (tok) *v2 = atof(tok);
  else { *v2=*v1; *st=0; return 1;}

  tok = strtok(NULL, del);
  if (tok){ /* we have step! */
    *st = *v2;
    *v2 = atof(tok);
  }
  else {
    *st=(*v2-*v1)/(20-1);
  }

  if (*v2==*v1){ *st=0; return 1;}

  if (*v1<*v2) *st = fabs(*st);
  else *st = -fabs(*st);

  if (*st==0) error("zero step");

  return (int)ceil(fabs((*v2-*v1)/(*st)))+1;
}

void
call_func(struct tab_t * F, char *argv[], int is_cmpl){
  const int MAXARG=8; /* max number of arguments */
  int N=1;            /* table length */
  int i,j;
  /* arguments: first and last value, step, current value */
  double A1[MAXARG], A2[MAXARG], dA[MAXARG], A[MAXARG];

  /* parse arguments */
  for (j=0; j<F->narg; j++){
    int n1 = parse_range(argv[j], &A1[j], &A2[j], &dA[j]);
    if (N!=1 && n1!=1 && n1!=N) error("wrong table size");
    if (n1!=1) N=n1;
    A[j] = atof(argv[j]);
  }

  /* print table header */
  if (N>1){
    printf("# %s\n", F->comm);
    printf("# %s ", F->args);
    if (is_cmpl) printf("%s\n", F->name);
    else         printf("Re(%s) Im(%s)\n", F->name, F->name);
  }

  /* call the function */
  for (i=0;i<N;i++){
    for (j=0; j<F->narg; j++){
      A[j] = A1[j] + i*dA[j];
      if (N>1) printf("%e ", A[j]);
    }
    if (is_cmpl) {
      // Calling fortran function which returns complex values
      // can depend on the compiler. Sometimes values returned through
      // an extra argument f(&ret, &A...). See -ff2c option if gfortran.
      complex (*ff)() = F->func;
      complex ret = {0,0};
      switch(F->narg){
        case 0: ret=ff(); break;
        case 1: ret=ff(&A[0]); break;
        case 2: ret=ff(&A[0],&A[1]); break;
        case 3: ret=ff(&A[0],&A[1],&A[2]); break;
        case 4: ret=ff(&A[0],&A[1],&A[2],&A[3]); break;
        case 5: ret=ff(&A[0],&A[1],&A[2],&A[3],&A[4]); break;
        case 6: ret=ff(&A[0],&A[1],&A[2],&A[3],&A[4],&A[5]); break;
        case 7: ret=ff(&A[0],&A[1],&A[2],&A[3],&A[4],&A[5],&A[6]); break;
        case 8: ret=ff(&A[0],&A[1],&A[2],&A[3],&A[4],&A[5],&A[6],&A[7]); break;
        default:
          error("Error: unsupported function prototype\n");
      }
      printf("%e %e\n", ret.r, ret.i);
    }
    else {
      double (*ff)() = F->func;
      double ret = 0;
      switch(F->narg){
        case 0: ret=ff(); break;
        case 1: ret=ff(&A[0]); break;
        case 2: ret=ff(&A[0],&A[1]); break;
        case 3: ret=ff(&A[0],&A[1],&A[2]); break;
        case 4: ret=ff(&A[0],&A[1],&A[2],&A[3]); break;
        case 5: ret=ff(&A[0],&A[1],&A[2],&A[3],&A[4]); break;
        case 6: ret=ff(&A[0],&A[1],&A[2],&A[3],&A[4],&A[5]); break;
        case 7: ret=ff(&A[0],&A[1],&A[2],&A[3],&A[4],&A[5],&A[6]); break;
        case 8: ret=ff(&A[0],&A[1],&A[2],&A[3],&A[4],&A[5],&A[6],&A[7]); break;
        default:
          error("Error: unsupported function prototype\n");
      }
      printf("%e\n", ret);
    }
  }
}


int
main(int argc, char *argv[]){
  char *cmd, *c;
  int i,j;

  /* no arguments -- print command list */
  if (argc<2){
    printf("### he3lib command line interface, slazav, 03.2014 ###\n");
    printf("### Command list:\n");
    for (i=0; func_tab[i].name!=NULL ; i++){
      struct tab_t * F = &func_tab[i];
      if (F->narg==0)
        printf("%-22s %s\n", F->name, F->comm);
      else
        printf("%s(%s)\t%s\n", F->name, F->args, F->comm);
    }
    exit(1);
  }

  /* first argument -- command name */
  cmd=argv[1];

  /* look for the command in a function list */
  for (i=0; func_tab[i].name!=NULL ; i++){
    struct tab_t *F = &func_tab[i];
    if (strcmp(F->name, cmd)!=0) continue;
    /* check number of arguments */
    if (argc-2!=F->narg){
      printf("%s\n", F->comm);
      printf("Not enough argumnets (%d expected)\n", F->narg);
      printf("Usage: %s %s\n", F->name, F->args);
      exit(1);
    }
    switch (F->type) {
    case CONST:
      printf("%e\n", *(double *)F->func);
      break;
    case DFUNC:
      call_func(F, argv+2, 0);
      break;
    case CFUNC:
      call_func(F, argv+2, 1);
      break;
    }
    exit(0);
  }
  printf("Unknown command: %s\n", cmd);
  exit(1);
}
