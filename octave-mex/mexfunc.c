#include "mex.h"
#include "math.h"
#include <string.h>

/* This is a universal matlab interface for
 * constants and functions of 1..5 aguments.
 * All arguments and returned values are double.
 * Vectors or matrices can be used as arguments,
 * vector and matrix arguments can be mixed
 * with constants:  res = func(1:10, 1, 2:11, 2)
 * will return a vector with length 10
 *
 * FUNC and NARGIN must be defined during compilation!
 *
 *                               slazav, 2013 */

/* I can either get function prototype from the h-file or create
   it here. I want to support h-file in a working condition
   so I use it. */
# if 1
#include "../he3.h"
#else
#if NARGIN == 0
extern double FUNC;
#endif
#if NARGIN == 1
double FUNC(double *a1);
#endif
#if NARGIN == 2
double FUNC(double *a1, double *a2);
#endif
#if NARGIN == 3
double FUNC(double *a1, double *a2, double *a3);
#endif
#if NARGIN == 4
double FUNC(double *a1, double *a2, double *a3,
            double *a4);
#endif
#if NARGIN == 5
double FUNC(double *a1, double *a2, double *a3,
            double *a4, double *a5);
#endif
#if NARGIN == 6
double FUNC(double *a1, double *a2, double *a3,
            double *a4, double *a5, double *a6);
#endif
#if NARGIN == 7
double FUNC(double *a1, double *a2, double *a3,
            double *a4, double *a5, double *a6,
            double *a7);
#endif
#if NARGIN == 8
double FUNC(double *a1, double *a2, double *a3,
            double *a4, double *a5, double *a6,
            double *a7, double *a8);
#endif
#endif

#include "../he3tab.h" /* for help messages */
/* stringification of the function name, see
  https://gcc.gnu.org/onlinedocs/cpp/Stringification.html*/
#define str1(x) str2(x)
#define str2(x) #x
#define sFUNC str1(FUNC)

void
mexFunction(int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[]){

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

  if (nrhs != NARGIN){
    /* look for the command in a function list */
    for (i=0; func_tab[i].name!=NULL ; i++){
      if ( strlen(func_tab[i].name) == strlen(sFUNC)-1 &&
           strncmp(func_tab[i].name, sFUNC, strlen(sFUNC)-1)==0) break;
    }
    struct tab_t *F = &func_tab[i];
    if (F->name){
      mexPrintf("# Function: %s: %s\n", F->name, F->comm);
      mexPrintf("# Usage: %s %s\n", F->name, F->args, F->comm);
    }
    mexErrMsgTxt("wrong number of arguments");
  }


  if (nrhs > 8)
    mexErrMsgTxt("functions with > 8 arguments are not supported in mexfunc.c");


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
#if NARGIN == 6
    out[i] = FUNC(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                  in[2]+(s[2]?0:i), in[3]+(s[3]?0:i),
                  in[4]+(s[4]?0:i), in[5]+(s[5]?0:i));
#endif
#if NARGIN == 7
    out[i] = FUNC(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                  in[2]+(s[2]?0:i), in[3]+(s[3]?0:i),
                  in[4]+(s[4]?0:i), in[5]+(s[5]?0:i),
                  in[6]+(s[6]?0:i));
#endif
#if NARGIN == 8
    out[i] = FUNC(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                  in[2]+(s[2]?0:i), in[3]+(s[3]?0:i),
                  in[4]+(s[4]?0:i), in[5]+(s[5]?0:i),
                  in[6]+(s[6]?0:i), in[7]+(s[7]?0:i));
#endif
  }
  return;
#endif
}