#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "../../zinm.h"

void
R_call_mle_nm
(
   int *x,
   int *d,
   int *n,
   double *ret
)
{

   zinm_par_t *par = new_zinm_par(*(size_t *)d);
   if (par == NULL) {
      fprintf(stderr, "error in function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      goto clean_and_return;
   }

   if (!mle_nm((uint32_t *)x, *(uint32_t *)d, *(uint32_t *)n, par)) {
      fprintf(stderr, "error in function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      goto clean_and_return;
   }

   // Success. Copy the (useful) values of 'par' to 'ret'.
   ret[0] = par->alpha;
   for (size_t i = 0 ; i < (*d)+1 ; i++) {
      ret[i+1] = par->p[i];
   }

clean_and_return:
   if (par != NULL) free(par);
   return;

}
