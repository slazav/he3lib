#include <math.h>
#include "he3tab.h"

/* This is a command line interface for the he3lib.
                                       slazav, 2014 */

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
      struct tab_t F=const_tab[i];
      printf("%-22s %s\n", F.name, F.comm);
    }
    /* functions */
    for (i=0; func_tab[i].name!=NULL ; i++){
      struct tab_t F=func_tab[i];
      printf("%s(%s)\t%s\n", F.name, F.args, F.comm);
    }
    exit(1);
  }

  /* first argument -- command name */
  cmd=argv[1];
  /* look for the command in a constant list */
  for (i=0; const_tab[i].name!=NULL ; i++){
    struct tab_t F=const_tab[i];
    if (strcmp(const_tab[i].name, cmd)==0){
      printf("%e\n", *(double *)const_tab[i].func);
      exit(0);
    }
  }
  /* look for the command in a function list */
  for (i=0; func_tab[i].name!=NULL ; i++){
    struct tab_t F=func_tab[i];
    double (*ff)() = F.func;
    if (strcmp(F.name, cmd)==0){
      double A[16]; /* large enough to keep all args*/
      /* check number of arguments */
      if (argc-2!=F.narg){
        printf("%s %s -- %s\n", F.name, F.args, F.comm);
        exit(1);
      }
      /* convert arguments to real numbers */
      for (j=0; j<F.narg; j++){
        A[j] = atof(argv[2+j]);
      }

      /* call the function */
      switch(F.narg){
        case 0: printf("%e\n", ff()); break;
        case 1: printf("%e\n", ff(&A[0])); break;
        case 2: printf("%e\n", ff(&A[0],&A[1])); break;
        case 3: printf("%e\n", ff(&A[0],&A[1],&A[2])); break;
        case 4: printf("%e\n", ff(&A[0],&A[1],&A[2],&A[3])); break;
        case 5: printf("%e\n", ff(&A[0],&A[1],&A[2],&A[3],&A[4])); break;
        default:
          printf("Error: unsupported function prototype\n");
          exit(1);
      }
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