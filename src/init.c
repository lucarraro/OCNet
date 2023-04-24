/* Reinhard Furrer,  fall 2019, based on spam's version and 'Writing R extensions'. */

//#include <R.h>
#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R_ext/Lapack.h>
//#include <Rinternals.h>

#include "invperm.c"

/* to get all functions:
   nm -g ./OCNet.so | grep " T "
*/



/* .Call calls */
extern SEXP _OCNet_paths_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _OCNet_WSC(SEXP, SEXP, SEXP, SEXP);
extern SEXP _OCNet_NN_OCN(SEXP, SEXP, SEXP, SEXP);
extern SEXP _OCNet_NN_river(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _OCNet_NN_FD(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"inv_permutation", (DL_FUNC) &inv_permutation, 1},
  {"_OCNet_paths_cpp", (DL_FUNC) &_OCNet_paths_cpp, 5},
  {"_OCNet_WSC", (DL_FUNC) &_OCNet_WSC, 4},
  {"_OCNet_NN_OCN", (DL_FUNC) &_OCNet_NN_OCN, 4},
  {"_OCNet_NN_river", (DL_FUNC) &_OCNet_NN_river, 6},
  {"_OCNet_NN_FD", (DL_FUNC) &_OCNet_NN_FD, 5},
  {NULL, NULL, 0}
};


/* .Fortran calls */

extern void F77_NAME(allinone    )( void *, void *, void *, void *, void *, void *,void *, void *,  void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(permandsolve)( void *, void *, void *, void *, void *, void *,void *, void *,  void *, void *, void *, void *, void *);


static const R_FortranMethodDef FortranEntries[] = {
    {"permandsolve",        (DL_FUNC) &F77_NAME(permandsolve ),13},
    {"allinone",            (DL_FUNC) &F77_NAME(allinone     ),21},
    {NULL, NULL, 0}
};

void R_init_OCNet(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

