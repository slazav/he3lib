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

/* Parse range setting in form <v1>:<v2>:<step>
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
  *st = tok? atof(tok):(*v2-*v1)/(20-1);

  if (*v1<*v2) *st = fabs(*st);
  else *st = -fabs(*st);

  if (*st==0) error("zero step");

  return (int)ceil(abs((*v2-*v1)/(*st)))+1;
}

void
call_func(struct tab_t * F, char *argv[]){
  const int MAXARG=5; /* max number of arguments */
  int N=1;            /* table length */
  int i,j;
  /* arguments: first and last value, step, current value */
  double A1[MAXARG], A2[MAXARG], dA[MAXARG], A[MAXARG];
  double (*ff)() = F->func;

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
    printf("# %s %s\n", F->args, F->name);
  }

  /* call the function */
  for (i=0;i<N;i++){
    for (j=0; j<F->narg; j++){
      A[j] = A1[j] + i*dA[j];
      if (N>1) printf("%e ", A[j]);
    }
    switch(F->narg){
      case 0: printf("%e\n", ff()); break;
      case 1: printf("%e\n", ff(&A[0])); break;
      case 2: printf("%e\n", ff(&A[0],&A[1])); break;
      case 3: printf("%e\n", ff(&A[0],&A[1],&A[2])); break;
      case 4: printf("%e\n", ff(&A[0],&A[1],&A[2],&A[3])); break;
      case 5: printf("%e\n", ff(&A[0],&A[1],&A[2],&A[3],&A[4])); break;
      default:
        error("Error: unsupported function prototype\n");
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
    /* constants */
    for (i=0; const_tab[i].name!=NULL ; i++){
      struct tab_t *F = &const_tab[i];
      printf("%-22s %s\n", F->name, F->comm);
    }
    /* functions */
    for (i=0; func_tab[i].name!=NULL ; i++){
      struct tab_t * F = &func_tab[i];
      printf("%s(%s)\t%s\n", F->name, F->args, F->comm);
    }
    exit(1);
  }

  /* first argument -- command name */
  cmd=argv[1];
  /* look for the command in a constant list */
  for (i=0; const_tab[i].name!=NULL ; i++){
    struct tab_t *F = &const_tab[i];
    if (strcmp(F->name, cmd)==0){
      printf("%e\n", *(double *)F->func);
      exit(0);
    }
  }
  /* look for the command in a function list */
  for (i=0; func_tab[i].name!=NULL ; i++){
    struct tab_t *F = &func_tab[i];
    if (strcmp(F->name, cmd)==0){
      /* check number of arguments */
      if (argc-2!=F->narg){
        printf("# %s\n", F->comm);
        printf("Usage: %s %s\n", F->name, F->args, F->comm);
        exit(1);
      }
      call_func(F, argv+2);
      exit(0);
    }
  }
  printf("Unknown command: %s\n", cmd);
  exit(1);
}




#if 0
/* Constants */
#if NARGIN == 0
  double *out;

  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  out = (double*)mxGetPr(plhs[0]);
  *out = FUNC;
  return;
#endif

/* Functions */
#if NARGIN > 0
  int m[NARGIN], n[NARGIN], s[NARGIN], i, maxm, maxn;
  double *in[NARGIN];
  double *out;

  if (nrhs != NARGIN)
    mexErrMsgTxt("wrong number of arguments");

  if (nrhs > 5)
    mexErrMsgTxt("functions with > 5 arguments are not supported in mexfunc.c");


  /* Get input arguments, calculate maximal size */
  maxm = maxn = 1;
  for (i=0; i<NARGIN; i++){
    if (!mxIsDouble(prhs[i]))
      mexErrMsgTxt("non-numeric argument");
    m[i] = mxGetM(prhs[i]);
    n[i] = mxGetN(prhs[i]);
    if (maxm<m[i]) maxm=m[i];
    if (maxn<n[i]) maxn=n[i];
    in[i] = mxGetPr(prhs[i]);
  }

  /* check that aguments are either constants or matrices of the same size */
  for (i=0; i<NARGIN; i++){
    s[i] = (m[i] == 1 && n[i] == 1); /* scalar/vector */
    if (!s[i] && (m[i] != maxm && n[i] != maxn))
      mexErrMsgTxt("dimensions of arguments must agree");
  }

  /* allocate space for output data */
  plhs[0] = mxCreateDoubleMatrix(maxm, maxn, mxREAL);
  if (plhs[0] == NULL)
    mexErrMsgTxt("can't allocate memory");
  out = mxGetPr(plhs[0]);

  /* calculate values */
  for (i=0; i<maxm*maxn; i++){
#if NARGIN == 1
    out[i] = FUNC(in[0]+(s[0]?0:i));
#endif
#if NARGIN == 2
    out[i] = FUNC(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i));
#endif
#if NARGIN == 3
    out[i] = FUNC(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                  in[2]+(s[2]?0:i));
#endif
#if NARGIN == 4
    out[i] = FUNC(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                  in[2]+(s[2]?0:i), in[3]+(s[3]?0:i));
#endif
#if NARGIN == 5
    out[i] = FUNC(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                  in[2]+(s[2]?0:i), in[3]+(s[3]?0:i),
                  in[4]+(s[4]?0:i));
#endif
  }
  return;
#endif
}
#endif