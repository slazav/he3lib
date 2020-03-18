#include "mex.h"
#include "math.h"
#include <string.h>

/* This is a universal matlab/octave mex interface for
 * all he3lib constants and functions.
 * All arguments and returned values are double.
 * Vectors or matrices can be used as arguments,
 * vector and matrix arguments can be mixed
 * with constants:  res = func(1:10, 1, 2:11, 2)
 * will return a vector with length 10
 *
 *                               slazav, 2013-2020 */

//extern "C" {
#include "he3tab.h"
//}

void
mexFunction(int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[]){

  struct tab_t *F = NULL;
  char fn[1024];

  if (nlhs < 1) return;
  // first argument - a string is a function name
  if (nrhs < 1)
    mexErrMsgTxt("name of function expected");
  if (!mxIsString(prhs[0]))
    mexErrMsgTxt("name of function in the first argument should be a string");
  if (mxGetString(prhs[0], fn, sizeof(fn))!=0)
    mexErrMsgTxt("can't convert string value");

  /* look for the command in a constant list */
  for (int i=0; const_tab[i].name!=NULL ; i++){
    if (strcmp(const_tab[i].name, fn)!=0) continue;
    double *out;
    F = &const_tab[i];
    if (nrhs-1 != 0){
      mexPrintf("%s: ", fn);
      mexErrMsgTxt("zero arguments expected");
    }
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    out = (double*)mxGetPr(plhs[0]);
    *out = *(double*)F->func;
    return;
  }
  /* look for the command in a function list */
  for (int i=0; func_tab[i].name!=NULL ; i++){
    if (strcmp(func_tab[i].name, fn)!=0) continue;
    F = &func_tab[i];
    break;
  }

  if (F == NULL){
    mexPrintf("%s: ", fn);
    mexErrMsgTxt("unknown function name");
  }

  if (nrhs-1 != F->narg){
    mexPrintf("%s: %d argument(s) expected (%s)\n", fn, F->narg, F->args);
    mexErrMsgTxt("wrong number of arguments");
  }

  /* Functions */
  int m[F->narg], n[F->narg], s[F->narg], maxm, maxn;
  double *in[F->narg];
  double *out;

  /* Get input arguments, calculate maximal size */
  maxm = maxn = 1;
  for (int i=0; i<F->narg; i++){
    if (!mxIsDouble(prhs[i+1]))
      mexErrMsgTxt("non-numeric argument");
    m[i] = mxGetM(prhs[i+1]);
    n[i] = mxGetN(prhs[i+1]);
    if (maxm<m[i]) maxm=m[i];
    if (maxn<n[i]) maxn=n[i];
    in[i] = mxGetPr(prhs[i+1]);
  }

  /* check that aguments are either constants or matrices of the same size */
  for (int i=0; i<F->narg; i++){
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
  double (*ff)() = F->func;
  for (int i=0; i<maxm*maxn; i++){
    switch(F->narg){
      case 1: out[i] = ((fun1_t)ff)(in[0]+(s[0]?0:i));
              break;
      case 2: out[i] = ((fun2_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i));
              break;
      case 3: out[i] = ((fun3_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                                    in[2]+(s[2]?0:i));
              break;
      case 4: out[i] = ((fun4_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                                    in[2]+(s[2]?0:i), in[3]+(s[3]?0:i));
              break;
      case 5: out[i] = ((fun5_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                                    in[2]+(s[2]?0:i), in[3]+(s[3]?0:i),
                                    in[4]+(s[4]?0:i));
              break;
      case 6: out[i] = ((fun6_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                                    in[2]+(s[2]?0:i), in[3]+(s[3]?0:i),
                                    in[4]+(s[4]?0:i), in[5]+(s[5]?0:i));
              break;
      case 7: out[i] = ((fun7_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                                    in[2]+(s[2]?0:i), in[3]+(s[3]?0:i),
                                    in[4]+(s[4]?0:i), in[5]+(s[5]?0:i),
                                    in[6]+(s[6]?0:i));
              break;
      case 8: out[i] = ((fun8_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                                    in[2]+(s[2]?0:i), in[3]+(s[3]?0:i),
                                    in[4]+(s[4]?0:i), in[5]+(s[5]?0:i),
                                    in[6]+(s[6]?0:i), in[7]+(s[7]?0:i));
              break;
      default:
        mexErrMsgTxt("Error: unsupported function prototype\n");
    }
  }
  return;
}
