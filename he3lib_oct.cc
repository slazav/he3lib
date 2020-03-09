#include <octave/oct.h>
#include <string.h>
#include <iostream>

/* This is a universal octave interface for
 * constants and functions with 1..5 aguments.
 * All arguments and returned values are double.
 * Vectors or matrices can be used as arguments,
 * vector and matrix arguments can be mixed
 * with constants:  res = func(1:10, 1, 2:11, 2)
 * will return a vector with length 10
 *
 *                               slazav, 2013-2020 */

extern "C" {
#include "he3.h"
}

DEFUN_DLD(he3lib, args, nargout, "he3 library") {

  int max_na = 8;

  std::string fn; // function name
  std::string fd; // function description
  int na = -1; // number of arguments
  double cnst;
  typedef double (*fun1_t)(double*);
  typedef double (*fun2_t)(double*, double*);
  typedef double (*fun3_t)(double*, double*, double*);
  typedef double (*fun4_t)(double*, double*, double*, double*);
  typedef double (*fun5_t)(double*, double*, double*, double*, double*);
  typedef double (*fun6_t)(double*, double*, double*, double*, double*, double*);
  typedef double (*fun7_t)(double*, double*, double*, double*, double*, double*, double*);
  typedef double (*fun8_t)(double*, double*, double*, double*, double*, double*, double*, double*);

  void *func;

  // first argument - a string is a function name
  if (args.length() < 1)
    error("name of function expected");

  if (!args(0).is_string())
    error("name of function in the first argument should be a string");

  fn = args(0).string_value();

  // read function table, find wich function we need,
  // number of arguments, description
  if (0) {}
  #define STR(s) #s
  #define CONST(name, args, descr)\
     else if (fn + '_' == STR(name)){ cnst = name; na = args; fd = descr;}
  #define FUNC(name, args, descr)\
     else if (fn + '_' == STR(name)){ func = (void*)name; na = args; fd = descr;}
  #include "he3funcs.h"
  else
    error("unknown function name: %s", fn.c_str());

  if (args.length() != 1 + na)
    error("%s -- %d arguments expected", fd.c_str(), na);

  if (args.length() > max_na+1)
    error("functions with > 8 arguments are not supported in octfunc.c");

  /* Constants */
  if (na == 0) {
    Matrix out(1,1);
    *out.fortran_vec()=cnst;
    return octave_value(out);
  }

  /* Functions */
  /* Get input arguments, calculate maximal size */
  NDArray ina[na];
  double  *in[na];
  bool s[na];
  dim_vector dv0(1,1);
  int numel = 1;

  for (int i=0; i<na; i++){
    if (!args(i+1).isnumeric())
      error("%s -- numeric matrix or scalar expected in argument %d", fd.c_str(), i);

    ina[i] = args(i+1).array_value(); // array
    in[i] = ina[i].fortran_vec();     // pointer to its raw data
    s[i] = ina[i].numel()==1;         // has 1 element?

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
    if (na==1) o[i] =
      ((fun1_t)func)(in[0]+(s[0]?0:i));
    if (na==2) o[i] =
      ((fun2_t)func)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i));
    if (na==3) o[i] =
      ((fun3_t)func)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                     in[2]+(s[2]?0:i));
    if (na==4) o[i] =
      ((fun4_t)func)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                     in[2]+(s[2]?0:i), in[3]+(s[3]?0:i));
    if (na==5) o[i] =
      ((fun5_t)func)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                     in[2]+(s[2]?0:i), in[3]+(s[3]?0:i),
                     in[4]+(s[4]?0:i));
    if (na==6) o[i] =
      ((fun6_t)func)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                     in[2]+(s[2]?0:i), in[3]+(s[3]?0:i),
                     in[4]+(s[4]?0:i), in[5]+(s[5]?0:i));
    if (na==7) o[i] =
      ((fun7_t)func)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                     in[2]+(s[2]?0:i), in[3]+(s[3]?0:i),
                     in[4]+(s[4]?0:i), in[5]+(s[5]?0:i),
                     in[6]+(s[6]?0:i));
    if (na==8) o[i] =
      ((fun8_t)func)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                     in[2]+(s[2]?0:i), in[3]+(s[3]?0:i),
                     in[4]+(s[4]?0:i), in[5]+(s[5]?0:i),
                     in[6]+(s[6]?0:i), in[7]+(s[7]?0:i));
  }
  return octave_value(out);
}
