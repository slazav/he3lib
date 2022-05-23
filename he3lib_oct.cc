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

  // first argument - a string is a function name
  if (args.length() < 1)
    error("name of function expected");
  if (!args(0).is_string())
    error("name of function in the first argument should be a string");
  std::string fn = args(0).string_value();

  /* look for the command in a function list */
  struct tab_t *F = NULL;
  for (int i=0; func_tab[i].name!=NULL ; i++){
    if (strcmp(func_tab[i].name, fn.c_str())!=0) continue;
    F = &func_tab[i];
    break;
  }

  if (F == NULL)
    error("%s: unknown function name", fn.c_str());
  if (args.length() != 1 + F->narg)
    error("%s: %d argument(s) expected (%s)", fn.c_str(), F->narg, F->args);


  /* constant */
  if (F->type == CONST){
    Matrix out(1,1);
    *out.fortran_vec()=*(double *)F->func;
    return octave_value(out);
  }

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

  if (F->type == DFUNC) {
    /* allocate space for output data */
    NDArray out(dv0);
    double *o = out.fortran_vec();

    /* calculate values */
    double (*ff)() = F->func;
    for (int i=0; i<numel; i++){
      OCTAVE_QUIT; // octave can break here;
      switch(F->narg){
        case 1: o[i] = ((dfun1_t)ff)(in[0]+(s[0]?0:i));
                break;
        case 2: o[i] = ((dfun2_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i));
                break;
        case 3: o[i] = ((dfun3_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                                    in[2]+(s[2]?0:i));
                break;
        case 4: o[i] = ((dfun4_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                                    in[2]+(s[2]?0:i), in[3]+(s[3]?0:i));
                break;
        case 5: o[i] = ((dfun5_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                                    in[2]+(s[2]?0:i), in[3]+(s[3]?0:i),
                                    in[4]+(s[4]?0:i));
                break;
        case 6: o[i] = ((dfun6_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                                    in[2]+(s[2]?0:i), in[3]+(s[3]?0:i),
                                    in[4]+(s[4]?0:i), in[5]+(s[5]?0:i));
                break;
        case 7: o[i] = ((dfun7_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                                    in[2]+(s[2]?0:i), in[3]+(s[3]?0:i),
                                    in[4]+(s[4]?0:i), in[5]+(s[5]?0:i),
                                    in[6]+(s[6]?0:i));
                break;
        case 8: o[i] = ((dfun8_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
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

  if (F->type == CFUNC) {
    /* allocate space for output data */
    ComplexNDArray out(dv0);
    complex *o = (complex*)out.fortran_vec();

    /* calculate values */
    complex (*ff)() = F->func;
    for (int i=0; i<numel; i++){
      OCTAVE_QUIT; // octave can break here;
      switch(F->narg){
        case 1: o[i] = ((cfun1_t)ff)(in[0]+(s[0]?0:i));
                break;
        case 2: o[i] = ((cfun2_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i));
                break;
        case 3: o[i] = ((cfun3_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                                    in[2]+(s[2]?0:i));
                break;
        case 4: o[i] = ((cfun4_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                                    in[2]+(s[2]?0:i), in[3]+(s[3]?0:i));
                break;
        case 5: o[i] = ((cfun5_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                                    in[2]+(s[2]?0:i), in[3]+(s[3]?0:i),
                                    in[4]+(s[4]?0:i));
                break;
        case 6: o[i] = ((cfun6_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                                    in[2]+(s[2]?0:i), in[3]+(s[3]?0:i),
                                    in[4]+(s[4]?0:i), in[5]+(s[5]?0:i));
                break;
        case 7: o[i] = ((cfun7_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
                                    in[2]+(s[2]?0:i), in[3]+(s[3]?0:i),
                                    in[4]+(s[4]?0:i), in[5]+(s[5]?0:i),
                                    in[6]+(s[6]?0:i));
                break;
        case 8: o[i] = ((cfun8_t)ff)(in[0]+(s[0]?0:i), in[1]+(s[1]?0:i),
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

  error("Error: unsupported function type\n");
}
