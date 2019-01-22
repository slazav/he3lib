#include <octave/oct.h>
#include <string.h>

/* This is a universal octave interface for
 * constants and functions with 1..5 aguments.
 * All arguments and returned values are double.
 * Vectors or matrices can be used as arguments,
 * vector and matrix arguments can be mixed
 * with constants:  res = func(1:10, 1, 2:11, 2)
 * will return a vector with length 10
 *
 * FUNC, FUNC_, DESCR and NARGIN must be defined during compilation!
 *
 *                               slazav, 2013-2018 */

/* I can either get function prototype from the h-file or create
   it here. I want to support h-file in a working condition
   so I use it. */

# if 1
extern "C" {
#include "../he3.h"
}
#else
#if NARGIN == 0
extern double FUNC_;
#endif
#if NARGIN == 1
double FUNC_(double *a1);
#endif
#if NARGIN == 2
double FUNC_(double *a1, double *a2);
#endif
#if NARGIN == 3
double FUNC_(double *a1, double *a2, double *a3);
#endif
#if NARGIN == 4
double FUNC_(double *a1, double *a2, double *a3,
            double *a4);
#endif
#if NARGIN == 5
double FUNC_(double *a1, double *a2, double *a3,
            double *a4, double *a5);
#endif
#if NARGIN == 6
double FUNC_(double *a1, double *a2, double *a3,
            double *a4, double *a5, double *a6);
#endif
#if NARGIN == 7
double FUNC_(double *a1, double *a2, double *a3,
            double *a4, double *a5, double *a6,
            double *a7);
#endif
#if NARGIN == 8
double FUNC_(double *a1, double *a2, double *a3,
            double *a4, double *a5, double *a6,
            double *a7, double *a8);
#endif
#endif

DEFUN_DLD(FUNC, args, nargout, "DESCR") {
  if (nargout==0) {}
  /*check number of arguments */
  if (args.length() != NARGIN)
    error("wrong number of arguments, %d expected", NARGIN);

  if (args.length() > 8)
    error("functions with > 8 arguments are not supported in octfunc.c");

/* Constants */
#if NARGIN == 0
  Matrix out(1,1);
  *out.fortran_vec()=FUNC_;
  return octave_value(out);
#endif

/* Functions */
#if NARGIN > 0
  /* Get input arguments, calculate maximal size */
  NDArray ina[NARGIN];
  double  *in[NARGIN];
  bool s[NARGIN];
  dim_vector dv0(1,1);
  int numel = 1;

  for (int i=0; i<NARGIN; i++){
    if (!args(i).isnumeric())
      error("numeric matrix or scalar expected in argument %d", i);

    ina[i] = args(i).array_value(); // array
    in[i] = ina[i].fortran_vec();   // pointer to its raw data
    s[i] = ina[i].numel()==1; // has 1 element?

    if (!s[i]) {
      if (numel!=1 && ina[i].ndims()!=dv0.ndims())
        error("wrong dimensions of argument %d", i);
      if (numel!=1 && ina[i].numel()!=numel)
        error("wrong number of elements in argument %d", i);
      dv0 = ina[i].dims();
      numel = ina[i].numel();
    }
  }

  /* allocate space for output data */
  NDArray out(dv0);
  double *o = out.fortran_vec();

  /* calculate values */
  for (int i=0; i<numel; i++){
    OCTAVE_QUIT; // octave can break here;
#if NARGIN == 1
    o[i] = FUNC_(in[0]+(s[0]?0:i));
#endif
#if NARGIN == 2
    o[i] = FUNC_(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i));
#endif
#if NARGIN == 3
    o[i] = FUNC_(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                in[2]+(s[2]?0:i));
#endif
#if NARGIN == 4
    o[i] = FUNC_(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                in[2]+(s[2]?0:i), in[3]+(s[3]?0:i));
#endif
#if NARGIN == 5
    o[i] = FUNC_(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                in[2]+(s[2]?0:i), in[3]+(s[3]?0:i),
                in[4]+(s[4]?0:i));
#endif
#if NARGIN == 6
    o[i] = FUNC_(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                in[2]+(s[2]?0:i), in[3]+(s[3]?0:i),
                in[4]+(s[4]?0:i), in[5]+(s[5]?0:i));
#endif
#if NARGIN == 7
    o[i] = FUNC_(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                in[2]+(s[2]?0:i), in[3]+(s[3]?0:i),
                in[4]+(s[4]?0:i), in[5]+(s[5]?0:i),
                in[6]+(s[6]?0:i));
#endif
#if NARGIN == 8
    o[i] = FUNC_(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                in[2]+(s[2]?0:i), in[3]+(s[3]?0:i),
                in[4]+(s[4]?0:i), in[5]+(s[5]?0:i),
                in[6]+(s[6]?0:i), in[7]+(s[7]?0:i));
#endif
  }
  return octave_value(out);
#endif
}
