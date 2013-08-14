#include "mex.h"
#include "math.h"
#include "../he3.h"

/* FUNC and NARGIN must be defined! */
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

#if NARGIN > 0

  int m[NARGIN], n[NARGIN], s[NARGIN], i, maxm, maxn;
  double *in[NARGIN];
  double *out;

  if (nrhs != NARGIN)
    mexErrMsgTxt("wrong number of arguments");

  if (nrhs > 5)
    mexErrMsgTxt("functions with > 5 arguments are not supported in mexfunc.c");

  maxm = nrhs? 0:1;
  maxn = nrhs? 0:1;

  for (i=0; i<NARGIN; i++){
    if (!mxIsDouble(prhs[i]))
      mexErrMsgTxt("non-numeric argument");
    m[i] = mxGetM(prhs[i]);
    n[i] = mxGetN(prhs[i]);
    if (maxm<m[i]) maxm=m[i];
    if (maxn<n[i]) maxn=n[i];
    in[i] = mxGetPr(prhs[i]);
  }
  for (i=0; i<NARGIN; i++){
    s[i] = (m[i] == 1 && n[i] == 1); /* scalar/vector */
    if (!s[i] && (m[i] != maxm && n[i] != maxn))
      mexErrMsgTxt("dimensions of arguments must agree");
  }

  plhs[0] = mxCreateDoubleMatrix(maxm, maxn, mxREAL);
  if (plhs[0] == NULL)
    mexErrMsgTxt("can't allocate memory");
  out = mxGetPr(plhs[0]);

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
