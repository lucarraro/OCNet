/* Reinhard Furrer, fall 2019 
   stripped-down version of 'inv_permutation' from Matrix.
*/  


#include <R.h>          /* includes Rconfig.h */
#include <Rdefines.h>   /* includes Rinternals.h  */
 
/* 
  Inverse Permutation
   inv_permutation <- function(p) { p[p] <- seq_along(p) ; p }
*/
static R_INLINE
SEXP inv_permutation(SEXP p_)
{
    int np = 1; 
    if(!isInteger(p_)) {p_ = PROTECT(coerceVector(p_, INTSXP)); np++; }
    int *p = INTEGER(p_), n = LENGTH(p_);
    SEXP val = PROTECT(allocVector(INTSXP, n));
    int *v = INTEGER(val);
    v--; 
    // alternative, yet MM prefers the version below:
    // for(int i=0; i < n; ) v[p[i]] = ++i;
    for(int i=0; i < n; ) {
	int j = p[i]; v[j] = ++i;
    }
    UNPROTECT(np);
    return val;
}
