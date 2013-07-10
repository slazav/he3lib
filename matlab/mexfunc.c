#include "mex.h"
#include "math.h"
#include "../he3.h"


void
mexFunction(int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[]){

/* Constants */
#if NARGIN == 0
  double *out;

  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  out = (double*)mxGetPr(plhs[0]);
  *out = he3_const_.FUNC;
  return;
#endif

/* Functions of one argument */
#if NARGIN == 1
  int m, n, i;
  double *in, *out;

  if (nrhs != 1)
    mexErrMsgTxt("1 argument required"); 

  if (!mxIsDouble(prhs[0]))
    mexErrMsgTxt("argument is non-numeric");

  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);

  in = mxGetPr(prhs[0]);
  plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
  if (plhs[0] == NULL)
    mexErrMsgTxt("can't allocate memory");
  out = mxGetPr(plhs[0]);

  for (i=0; i<m*n; i++){
    out[i] = FUNC(in+i);
    if (out[i]<0) out[i]=NAN;
  }
  return;
#endif

/* Functions of two arguments */
#if NARGIN == 2
  int m, n, m2, n2, i;
  double *in1, *in2, *out;
  int mode;

  if (nrhs != 2)
    mexErrMsgTxt("2 arguments required");

  if (!mxIsDouble(prhs[0]))
    mexErrMsgTxt("argument 1 is non-numeric");

  if (!mxIsDouble(prhs[1]))
    mexErrMsgTxt("argument 2 is non-numeric");

  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);

  m2 = mxGetM(prhs[1]);
  n2 = mxGetN(prhs[1]);

  if (m == m2 && n == n2){
    mode=0;
  }
  else if (m==1 && n==1 && m2>0 && n2>0){ /* 1st arg - constant, 2nd - vector*/
    m = m2; n = n2; mode = 1;
  }
  else if (m2==1 && n2==1 && m>0 && n>0){ /* 2st arg - constant, 1nd - vector*/
    mode = 2;
  }
  else
    mexErrMsgTxt("dimensions of arguments must agree");

  in1 = mxGetPr(prhs[0]);
  in2 = mxGetPr(prhs[1]);
  plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
  if (plhs[0] == NULL)
    mexErrMsgTxt("can't allocate memory");
  out = mxGetPr(plhs[0]);

  for (i=0; i<m*n; i++){
    out[i] = FUNC(in1+(mode==1?0:i), in2+(mode==2?0:i));
    if (out[i]<0) out[i]=NAN;
  }
  return;
#endif
}
