#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "../../zinm.h"


// Register method 'R_call_mle_nm()'.
SEXP R_call_mle_nm (SEXP, SEXP, SEXP);

R_CallMethodDef callMethods[] = {
   {"R_call_mle_nm", (DL_FUNC) &R_call_mle_nm, 1},
   {NULL, NULL, 0},
};

SEXP
R_call_mle_nm
(
   SEXP x_,
   SEXP d_,
   SEXP n_
)
{

   const uint32_t *x = (uint32_t *) INTEGER(x_);
   const uint32_t  d = (uint32_t) INTEGER(d_)[0];
   const uint32_t  n = (uint32_t) INTEGER(n_)[0];

   zinm_par_t *par = new_zinm_par(d);
   if (par == NULL) {
      fprintf(stderr, "error in function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      goto clean_and_return;
   }

   if (!mle_nm(x, d, n, par)) {
      fprintf(stderr, "error in function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      goto clean_and_return;
   }

   // Success. Copy the (useful) values of 'par' to 'ret'.
   SEXP RET;
   PROTECT(RET = allocVector(REALSXP, d+2));
   REAL(RET)[0] = par->alpha;
   for (size_t i = 0 ; i < d+2 ; i++) {
      REAL(RET)[i+1] = par->p[i];
   }

   UNPROTECT(1);

clean_and_return:
   if (par != NULL) free(par);
   return RET;

}
