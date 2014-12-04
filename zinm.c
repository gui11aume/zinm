#include "zinm.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


double
eval_nb_f
(
         double a,
   const tab_t *tab
)
{

   double retval;
   double prev;
   // Convenience variables.
   const unsigned int *val = tab->val;
   const unsigned int *num = tab->num;

   size_t nobs = num[0];
   double mean = num[0] * val[0];

   prev = digamma(a + val[0]);
   retval = num[0] * prev;
   // Iterate over the occurrences and compute the new value
   // of digamma either by the recurrence relation, or by
   // a new call to 'digamma()', whichever is faster.
   for (size_t i = 1 ; i < tab->size ; i++) {
      nobs += num[i];
      mean += num[i] * val[i];
      prev = (val[i] - val[i-1] == 1) ?
         prev + 1.0 / (a-1 + val[i]) :
         digamma(a + val[i]);
      retval += num[i] * prev;
   }

   mean /= nobs;
   retval += nobs*(log(a) - digamma(a) - log(a + mean));

   return retval;

}


double
eval_nb_dfda
(
   double a,
   const tab_t *tab
)
{

   double retval;
   double prev;
   // Convenience variables.
   const unsigned int *val = tab->val;
   const unsigned int *num = tab->num;

   size_t nobs = num[0];
   double mean = num[0] * val[0];

   prev = trigamma(a + val[0]);
   retval = num[0] * prev;
   // Iterate over the occurrences and compute the new value
   // of trigamma either by the recurrence relation, or by
   // a new call to 'trigamma()', whichever is faster.
   for (size_t i = 1 ; i < tab->size ; i++) {
      nobs += num[i];
      mean += num[i] * val[i];
      prev = (val[i] - val[i-1] == 1) ?
         prev - 1.0 / sq(a-1 +val[i]) :
         trigamma(a + val[i]);
      retval += num[i] * prev;
   }

   mean /= nobs;
   retval += nobs*(mean/(a*(a+mean)) - trigamma(a));

   return retval;

}


double
eval_zinm_f
(
   double a,
   double p,
   unsigned int nz,
   double m
)
{
   return a*nz / (p*(1-pow(p,a))) - m/(1-p);
}


double
eval_zinm_g
(
         double a,
         double p,
   const tab_t *tab
)
{

   // Convenience variables.
   const unsigned int *val = tab->val;
   const unsigned int *num = tab->num;

   unsigned int nz = 0;
   double retval = 0.0;
   double prev = digamma(a + val[0]);

   if (val[0] > 0) {
      retval += num[0] * prev;
      nz += num[0];
   }

   // Iterate over the occurrences and compute the new value
   // of digamma either by the recurrence relation, or by
   // a new call to 'digamma()', whichever is faster.
   const size_t imin = val[0] == 0 ? 1 : 0;
   for (size_t i = imin ; i < tab->size ; i++) {
      nz += num[i];
      prev = (val[i] - val[i-1] == 1) ?
         prev + 1.0 / (a-1 + val[i]) :
         digamma(a + val[i]);
      retval += num[i] * prev;
   }

   retval += nz*(log(p) / (1-pow(p,a)) - digamma(a));
   return retval;

}


double
eval_zinm_dfda
(
   double a,
   double p,
   unsigned int nz
)
{
   const double ppa = pow(p,a);
   return nz*(1-ppa+a*ppa*log(p)) / (p*sq(1-ppa));
}

double
eval_zinm_dfdp
(
   double a,
   double p,
   unsigned int nz,
   double m
)
{
   const double ppa = pow(p,a);
   return -(nz*a*(1-(a+1)*ppa) / sq(p*(1-ppa)) + m/sq(1-p));
}


double
eval_zinm_dgda
(
         double a,
         double p,
   const tab_t *tab
)
{
   
   // Convenience variables.
   const unsigned int *val = tab->val;
   const unsigned int *num = tab->num;
   const double ppa = pow(p,a);

   unsigned int nz = 0;
   double retval = 0.0;
   double prev = trigamma(a + val[0]);

   if (val[0] > 0) {
      retval += num[0] * prev;
      nz += num[0];
   }

   // Iterate over the occurrences and compute the new value
   // of digamma either by the recurrence relation, or by
   // a new call to 'trigamma()', whichever is faster.
   const size_t imin = val[0] == 0 ? 1 : 0;
   for (size_t i = imin ; i < tab->size ; i++) {
      nz += num[i];
      prev = (val[i] - val[i-1] == 1) ?
         prev - 1.0 / sq(a-1 + val[i]) :
         trigamma(a + val[i]);
      retval += num[i] * prev;
   }

   retval += nz*(sq(log(p))*ppa / sq(1-ppa) - trigamma(a));
   return retval;

}


double
ll_zinm
(
         double a,
         double p,
         double pi,
   const tab_t *tab
)
{

   // Convenience variables.
   const unsigned int *val = tab->val;
   const unsigned int *num = tab->num;
   const unsigned int z0 = val[0] == 0 ? num[0] : 0;
   const double logp_ = log(1-p);

   unsigned int nobs = z0;
   double retval = z0*log(pi*pow(p,a) + 1-pi);

   const size_t imin = val[0] == 0 ? 1 : 0;
   for (size_t i = imin ; i < tab->size ; i++) {
      retval += num[i] * (lgamma(a+val[i]) + val[i]*logp_);
      nobs += num[i];
   }
   // Ignore constant factorial terms.
   retval += (nobs-z0) * (a*log(p) - lgamma(a) + log(pi));

   return retval;

}


zinm_par_t *
mle_zinm
(
   size_t *x,
   size_t dim,
   size_t nobs
)
{

   tab_t *tab = tabulate(x, dim, nobs);
   double *means = malloc(dim * sizeof(double));
   if (means == NULL) {
      fprintf(stderr, "memory error: %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }

   // Compute the means in all dimensions and the grand mean.
   compute_means(x, dim, nobs, means);
   double mean = 0.0;
   for (size_t i = 0 ; i < dim ; i++) mean += means[i];

   // Extract the number of all-zero observaions.
   const unsigned int z0 = tab->val[0] == 0 ? tab->num[0] : 0;

   double deficit[11] = {0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1};
   double init_a[12] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1};
   double init_p[12] = {0,0,0,0,0,0,0,0,0,0,0,.5};

   // Deplete some 0s from the observations, compute alpha
   // and p0 with standard negative binomial estimation
   // and keep the values to be used as initial conditions.
   for (size_t i = 0 ; i < 11 ; i++) {
      double newmean = mean;
      if (tab->val[0] == 0) {
         tab->num[0] = z0 * (1-deficit[i]);
         newmean /= (1.0 - z0*deficit[i]/nobs);
      }
      double alpha = nb_est_alpha(tab);
      init_a[i] = alpha;
      init_p[i] = alpha / (alpha + newmean);
   }

   // Rest 'tab'.
   if (tab->val[0] == 0) tab->num[0] = z0;

   // Try initial conditions. Number 12 is a safety in case
   // all the rest failed during the first phase.
   double max_loglik = -1.0/0.0; // -inf
   zinm_par_t *par = new_zinm_par(dim);
   if (par == NULL) {
      fprintf(stderr, "memory error: %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }
   for (size_t i = 0 ; i < 12 ; i++) {

      if (init_a[i] < 0) continue;  // Skip failures.

      double a = init_a[i];
      double p = init_p[i];

      double grad;
      unsigned int iter = 0;

      double f = eval_zinm_f(a, p, nobs-z0, mean);
      double g = eval_zinm_g(a, p, tab);

      // Newton-Raphson iterations.
      while ((grad = f*f+g*g) > ZINM_TOL && iter++ < ZINM_MAXITER) {

         double dfda, dfdp, dgda, dgdp;
         dfda = dgdp = eval_zinm_dfda(a, p, nobs-z0);
         dfdp = eval_zinm_dfdp(a, p, nobs-z0, mean);
         dgda = eval_zinm_dgda(a, p, tab);

         double denom = dfdp*dgda - dfda*dgdp;
         double da = (f*dgdp - g*dfdp) / denom;
         double dp = (g*dfda - f*dgda) / denom;
         f = eval_zinm_f(a, p, nobs-z0, mean);
         g = eval_zinm_g(a, p, tab);

         for (int j = 0 ; j < ZINM_MAXITER && f*f+g*g > grad ; j++) {
            da /= 2;
            dp /= 2;
            f = eval_zinm_f(a, p, nobs-z0, mean);
            g = eval_zinm_g(a, p, tab);
         }

         a = a+da;
         p = p+dp;

      }

      double pi = (nobs-z0) / nobs / (1-pow(p,a));
      double loglik = ll_zinm(a, p, pi, tab);
      if (loglik > max_loglik) {
         max_loglik = loglik;
         par->alpha = a;
         par->pi = pi;
         par->p[0] = p;
      }
            
   }

   free(means);
   return par;

}

double
nb_est_alpha
(
   tab_t *tab
)
{

   // Find upper and lower bouds for a(lpha).
   double a = 1.0;
   double a_lo;
   double a_hi;
   if (eval_nb_f(a, tab) < 0) {
      a /= 2;
      while (eval_nb_f(a, tab) < 0) a /= 2;
      a_lo = a;
      a_hi = a*2;
   }
   else {
      a *= 2;
      while (eval_nb_f(a, tab) > 0) a *= 2;
      a_lo = a/2;
      a_hi = a;
   }

   // Input is pathological.
   if (a_lo > 128) return -1.0;

   double new_a = (a_lo + a_hi) / 2;
   for (int i = 0 ; i < ZINM_MAXITER ; i++) {
      a = (new_a < a_lo || new_a > a_hi) ?
         (a_lo + a_hi) / 2 :
         new_a;
      double f = eval_nb_f(a, tab);
      if (f < 0) a_hi = a; else a_lo = a;
      if ((a_hi - a_lo) < ZINM_TOL) break;
      double dfda = eval_nb_dfda(a, tab);
      new_a = a - f / dfda;
   }

   return a;

}

zinm_par_t *
mle_nm
(
   size_t *x,
   size_t dim,
   size_t nobs
)
{

   // Estimate alpha.
   tab_t *tab = tabulate(x, dim, nobs);
   double alpha = nb_est_alpha(tab);
   free(tab);

   // Return NULL if failed.
   if (alpha < 0) return NULL;

   // Compute the means in all dimensions.
   double *means = malloc(dim * sizeof(double));
   if (means == NULL) {
      fprintf(stderr, "memory error: %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }

   compute_means(x, dim, nobs, means);

   zinm_par_t *par = new_zinm_par(dim);
   if (par == NULL) {
      free(means);
      fprintf(stderr, "memory error: %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }

   double mean = 0;
   for (size_t i = 0 ; i < dim ; i++) mean += means[i];
   par->alpha = alpha;
   par->p[0] = alpha / (alpha + mean);
   for (size_t i = 1 ; i < dim+1 ; i++) {
      par->p[i] = par->p[0] / alpha * means[i-1];
   }

   free(means);

   return par;

}

// helper functions //

tab_t *
tabulate
(
   size_t * x,
   size_t   dim,
   size_t   nobs
)
{

   histo_t *histo = new_histo();
   if (histo == NULL) return NULL;
   for (size_t i = 0 ; i < nobs ; i++) {
      size_t sum = 0;
      for (size_t j = 0 ; j < dim ; j++) sum += x[i+j*nobs];
      if (histo_push(&histo, sum)) {
         free(histo);
         return NULL;
      }
   }

   tab_t *tab = compress_histo(histo);
   free(histo);
   return tab;

}


void
compute_means
(
   size_t * x,
   size_t   dim,
   size_t   nobs,
   double * means
)
{

   memset(means, 0, dim*sizeof(double));
   for (size_t i = 0 ; i < nobs ; i++) {
   for (size_t j = 0 ; j < dim ; j++) {
      means[j] += x[i+dim*j];
   }
   }

   for (size_t j = 0 ; j < dim ; j++) {
      means[j] /= nobs;
   }

   return;

}


zinm_par_t *
new_zinm_par
(
   size_t r
)
{

   zinm_par_t *new = calloc(1, sizeof(zinm_par_t) + (r+1)*sizeof(double));
   if (new == NULL) {
      fprintf(stderr, "memory error: %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }
   new->r = r;
   new->pi = 1.0;

   return new;

}

histo_t *
new_histo
(void)
{

   size_t initsize = sizeof(histo_t) + HISTO_INIT_SIZE * sizeof(int);
   histo_t *new = calloc(1, initsize);
   if (new == NULL) {
      fprintf(stderr, "memory error: %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }
   new->size = HISTO_INIT_SIZE;

   return new;

}


int
histo_push
(
   histo_t **histo_addr,
   size_t val
)
{

   // Convenience variable.
   histo_t *histo = *histo_addr;
   if (val >= histo->size) {
      size_t newsize = 2*val * (sizeof(int));
      histo_t *new = realloc(histo, sizeof(histo_t) + newsize);
      if (new == NULL) {
         fprintf(stderr, "memory error: %s:%d\n", __FILE__, __LINE__);
         return 1;
      }
      *histo_addr = histo = new;
      size_t added_size = (2*val - histo->size) * sizeof(int);
      memset(histo->num + histo->size, 0, added_size);
      histo->size = 2*val;
   }

   histo->num[val]++;
   return 0;

}

tab_t *
compress_histo
(
   histo_t *histo
)
{

   size_t size = 0;
   for (size_t i = 0 ; i < histo->size ; i++) {
      size += (histo->num[i] > 0);
   }
   tab_t *new = malloc(sizeof(tab_t) + 2*size * sizeof(int));
   if (new == NULL) {
      fprintf(stderr, "memory error: %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }
   new->size = size;
   new->num = new->val + size;

   size_t j = 0;
   for (size_t i = 0 ; i < histo->size ; i++) {
      if (histo->num[i] > 0) {
         new->val[j] = i;
         new->num[j] = histo->num[i];
         j++;
      }
   }

   return new;

}

// The code below was copied from the following link
// http://pmtksupport.googlecode.com/svn/trunk/lightspeed2.3/util.c
// Written by Tom Minka (unless otherwise noted).

/* The digamma function is the derivative of gammaln.

   Reference:
    J Bernardo,
    Psi ( Digamma ) Function,
    Algorithm AS 103,
    Applied Statistics,
    Volume 25, Number 3, pages 315-317, 1976.

    From http://www.psc.edu/~burkardt/src/dirichlet/dirichlet.f
    (with modifications for negative numbers and extra precision)
*/
double
digamma
(
   double x
)
{
  double neginf = -INFINITY;
  static const double c = 12,
    digamma1 = -0.57721566490153286,
    trigamma1 = 1.6449340668482264365, /* pi^2/6 */
    s = 1e-6,
    s3 = 1./12,
    s4 = 1./120,
    s5 = 1./252,
    s6 = 1./240,
    s7 = 1./132;
    //s8 = 691./32760,
    //s9 = 1./12,
    //s10 = 3617./8160;
  double result;
  /* Illegal arguments */
  if((x == neginf) || isnan(x)) {
    return NAN;
  }
  /* Singularities */
  if((x <= 0) && (floor(x) == x)) {
    return neginf;
  }
  /* Negative values */
  /* Use the reflection formula (Jeffrey 11.1.6):
   * digamma(-x) = digamma(x+1) + pi*cot(pi*x)
   *
   * This is related to the identity
   * digamma(-x) = digamma(x+1) - digamma(z) + digamma(1-z)
   * where z is the fractional part of x
   * For example:
   * digamma(-3.1) = 1/3.1 + 1/2.1 + 1/1.1 + 1/0.1 + digamma(1-0.1)
   *               = digamma(4.1) - digamma(0.1) + digamma(1-0.1)
   * Then we use
   * digamma(1-z) - digamma(z) = pi*cot(pi*z)
   */
  if(x < 0) {
    return digamma(1-x) + M_PI / tan(-M_PI*x);
  }
  /* Use Taylor series if argument <= S */
  if(x <= s) return digamma1 - 1/x + trigamma1*x;
  /* Reduce to digamma(X + N) where (X + N) >= C */
  result = 0;
  while(x < c) {
    result -= 1/x;
    x++;
  }
  /* Use de Moivre's expansion if argument >= C */
  /* This expansion can be computed in Maple via asympt(Psi(x),x) */
  if(x >= c) {
    double r = 1/x;
    result += log(x) - 0.5*r;
    r *= r;
    result -= r * (s3 - r * (s4 - r * (s5 - r * (s6 - r * s7))));
  }
  return result;
}

/* The trigamma function is the derivative of the digamma function.

   Reference:

    B Schneider,
    Trigamma Function,
    Algorithm AS 121,
    Applied Statistics, 
    Volume 27, Number 1, page 97-99, 1978.

    From http://www.psc.edu/~burkardt/src/dirichlet/dirichlet.f
    (with modification for negative arguments and extra precision)
*/
double trigamma(double x)
{
  double neginf = -INFINITY,
    small = 1e-4,
    large = 8,
    trigamma1 = 1.6449340668482264365, /* pi^2/6 = Zeta(2) */
    tetragamma1 = -2.404113806319188570799476,  /* -2 Zeta(3) */
    b2 =  1./6,  /* B_2 */
    b4 = -1./30, /* B_4 */
    b6 =  1./42, /* B_6 */
    b8 = -1./30, /* B_8 */
    b10 = 5./66; /* B_10 */
  double result;
  /* Illegal arguments */
  if((x == neginf) || isnan(x)) {
    return NAN;
  }
  /* Singularities */
  if((x <= 0) && (floor(x) == x)) {
    return neginf;
  }
  /* Negative values */
  /* Use the derivative of the digamma reflection formula:
   * -trigamma(-x) = trigamma(x+1) - (pi*csc(pi*x))^2
   */
  if(x < 0) {
    result = M_PI/sin(-M_PI*x);
    return -trigamma(1-x) + result*result;
  }
  /* Use Taylor series if argument <= small */
  if(x <= small) {
    return 1/(x*x) + trigamma1 + tetragamma1*x;
  }
  result = 0;
  /* Reduce to trigamma(x+n) where ( X + N ) >= B */
  while(x < large) {
    result += 1/(x*x);
    x++;
  }
  /* Apply asymptotic formula when X >= B */
  /* This expansion can be computed in Maple via asympt(Psi(1,x),x) */
  if(x >= large) {
    double r = 1/(x*x);
    result += 0.5*r + (1 + r*(b2 + r*(b4 + r*(b6 + r*(b8 + r*b10)))))/x;
  }
  return result;
}
