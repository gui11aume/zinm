#include <assert.h>
#include <float.h>	/* for FLT_MAX and FLT_MIN*/
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct _gen_beta_param
    {   double a,b;		/* the parameters of the beta distribution */
        /* The remaining fields are precomputed values that can
	 * save some computation when many numbers are drawn from the
	 * same distribution.
	 */
	double min_ab, max_ab;  /* min(a,b) and max(a,b) */
    	double sum_ab;		/* a+b */
	double param[3];	/* various precomputed parameters,
				 * different for each generation method
				 */
    } gen_beta_param;

/* Initialize the variables of the already allocated generator in gen .
 * This must be called before gen can be used.
 */
void gen_beta_initialize(gen_beta_param *gen, double a, double b);


/* Generate one number from the distribution specified by gen */
double gen_gauss (void);
double gen_beta(double a, double b);
double gen_gamma (double a, double b);
