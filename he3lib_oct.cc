#include <octave/oct.h>
#include <string.h>
#include <iostream>

/* This is a universal octave interface for
 * all he3lib constants and functions.
 * All arguments and returned values are double.
 * Vectors or matrices can be used as arguments,
 * vector and matrix arguments can be mixed
 * with constants:  res = func(1:10, 1, 2:11, 2)
 * will return a vector with length 10
 *
 *                               slazav, 2013-2020 */

extern "C" {
#include "he3tab.h"
}

DEFUN_DLD(he3lib, args, nargout, "he3 library") {

  int max_na = 8;

  double cnst;
  typedef double (*fun1_t)(double*);
  typedef double (*fun2_t)(double*, double*);
  typedef double (*fun3_t)(double*, double*, double*);
  typedef double (*fun4_t)(double*, double*, double*, double*);
  typedef double (*fun5_t)(double*, double*, double*, double*, double*);
  typedef double (*fun6_t)(double*, double*, double*, double*, double*, double*);
  typedef double (*fun7_t)(double*, double*, double*, double*, double*, double*, double*);
  typedef double (*fun8_t)(double*, double*, double*, double*, double*, double*, double*, double*);

  // first argument - a string is a function name
  if (args.length() < 1)
    error("name of function expected");
  if (!args(0).is_string())
    error("name of function in the first argument should be a string");
  std::string fn = args(0).string_value();

  /* look for the command in a constant list */
  struct tab_t *F = NULL;
  for (int i=0; const_tab[i].name!=NULL ; i++){
    if (strcmp(const_tab[i].name, fn.c_str())==0){
      F = &const_tab[i];
      Matrix out(1,1);
      *out.fortran_vec()=*(double *)F->func;
      return octave_value(out);
    }
  }

  /* look for the command in a function list */
  for (int i=0; func_tab[i].name!=NULL ; i++){
    if (strcmp(func_tab[i].name, fn.c_str())==0){
      F = &func_tab[i];
      break;
    }
  }

  if (F == NULL)
    error("%s: unknown function name", fn.c_str());
  if (args.length() != 1 + F->narg)
    error("%s: %d argument(s) expected (%s)", fn.c_str(), F->narg, F->args);


  /* Get input arguments, calculate maximal size */
  NDArray ina[F->narg];
  double  *in[F->narg];
  bool s[F->narg];
  dim_vector dv0(1,1);
  int numel = 1;
  for (int i=0; i<F->narg; i++){
    if (!args(i+1).isnumeric())
      error("%s: numeric matrix or scalar expected in argument %d", fn.c_str(), i);

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
  double (*ff)() = F->func;
  for (int i=0; i<numel; i++){
    OCTAVE_QUIT; // octave can break here;
    switch(F->narg){
      case 1: o[i] = ((fun1_t)ff)(in[0]+(s[0]?0:i));
              break;
      case 2: o[i] = ((fun2_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i));
              break;
      case 3: o[i] = ((fun3_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                                  in[2]+(s[2]?0:i));
              break;
      case 4: o[i] = ((fun4_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                                  in[2]+(s[2]?0:i), in[3]+(s[3]?0:i));
              break;
      case 5: o[i] = ((fun5_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                                  in[2]+(s[2]?0:i), in[3]+(s[3]?0:i),
                                  in[4]+(s[4]?0:i));
              break;
      case 6: o[i] = ((fun6_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                                  in[2]+(s[2]?0:i), in[3]+(s[3]?0:i),
                                  in[4]+(s[4]?0:i), in[5]+(s[5]?0:i));
              break;
      case 7: o[i] = ((fun7_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                                  in[2]+(s[2]?0:i), in[3]+(s[3]?0:i),
                                  in[4]+(s[4]?0:i), in[5]+(s[5]?0:i),
                                  in[6]+(s[6]?0:i));
              break;
      case 8: o[i] = ((fun8_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                                  in[2]+(s[2]?0:i), in[3]+(s[3]?0:i),
                                  in[4]+(s[4]?0:i), in[5]+(s[5]?0:i),
                                  in[6]+(s[6]?0:i), in[7]+(s[7]?0:i));
              break;
      default:
        error("Error: unsupported function prototype\n");
    }
  }
  return octave_value(out);
}
