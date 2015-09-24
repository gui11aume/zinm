#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "zinm.h"

// Constants.
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define HASH_SIZE 105943
#define ZINM_MAXITER 32
#define ZINM_TOL 1e-6

#define sq(x) ((x)*(x))

// Type definitions.
struct tab_t;
struct kv_t;

typedef struct tab_t tab_t;
typedef struct kv_t kv_t;

struct tab_t {
   size_t     size;
   uint32_t * num;
   uint32_t    val[];
};

struct kv_t {
   uint32_t  key;
   uint32_t val;
};


// Declaration of local functions.
void    compute_means(const uint32_t *, uint32_t, uint32_t, double *);
double  eval_nb_dfda(double, const tab_t *);
double  eval_nb_f(double, const tab_t *);
double  eval_zinm_dfda(double, double, uint32_t);
double  eval_zinm_dfdp(double, double, uint32_t, double);
double  eval_zinm_dgda(double, double, const tab_t *);
double  eval_zinm_f(double, double, uint32_t, double);
double  eval_zinm_g(double, double, const tab_t *);
int     update(kv_t *, uint32_t);
double  ll_zinm(double, double, double, const tab_t *);
double  nb_est_alpha(tab_t *);
kv_t  * new_counter(void);
tab_t * tabulate(const uint32_t *, uint32_t, uint32_t);

// Declaration of mathematical functions.
double digamma(double);
double trigamma(double);



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
   const uint32_t *val = tab->val;
   const uint32_t *num = tab->num;

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
   const uint32_t *val = tab->val;
   const uint32_t *num = tab->num;

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
   uint32_t nz,
   double sum
)
{
   return a*nz / (p*(1-pow(p,a))) - sum/(1-p);
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
   const uint32_t *val = tab->val;
   const uint32_t *num = tab->num;

   uint32_t nz = 0;
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
   uint32_t nz
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
   uint32_t nz,
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
   const uint32_t *val = tab->val;
   const uint32_t *num = tab->num;
   const double ppa = pow(p,a);

   uint32_t nz = 0;
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
   const uint32_t *val = tab->val;
   const uint32_t *num = tab->num;
   const uint32_t z0 = val[0] == 0 ? num[0] : 0;
   const double logp_ = log(1-p);

   uint32_t nobs = z0;
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
   uint32_t *x,
   uint32_t dim,
   uint32_t nobs,
   zinm_par_t *par
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
   double sum = mean * nobs;

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

   // Reset 'tab'.
   if (tab->val[0] == 0) tab->num[0] = z0;

   // Try initial conditions. Number 12 is a safety in case
   // all the rest failed during the first phase.
   double max_loglik = -1.0/0.0; // -inf

   for (size_t i = 0 ; i < 12 ; i++) {

      if (init_a[i] < 0) continue;  // Skip failures.

      double a = init_a[i];
      double p = init_p[i];

      double grad;
      unsigned int iter = 0;

      double f = eval_zinm_f(a, p, nobs-z0, sum);
      double g = eval_zinm_g(a, p, tab);

      // Newton-Raphson iterations.
      while ((grad = f*f+g*g) > sq(ZINM_TOL) && iter++ < ZINM_MAXITER) {

         double dfda, dfdp, dgda, dgdp;
         dfda = dgdp = eval_zinm_dfda(a, p, nobs-z0);
         dfdp = eval_zinm_dfdp(a, p, nobs-z0, sum);
         dgda = eval_zinm_dgda(a, p, tab);

         double denom = dfdp*dgda - dfda*dgdp;
         double da = (f*dgdp - g*dfdp) / denom;
         double dp = (g*dfda - f*dgda) / denom;
         f = eval_zinm_f(a+da, p+dp, nobs-z0, sum);
         g = eval_zinm_g(a+da, p+dp, tab);

         for (int j = 0 ; j < ZINM_MAXITER && f*f+g*g > grad ; j++) {
            da /= 2;
            dp /= 2;
            f = eval_zinm_f(a+da, p+dp, nobs-z0, sum);
            g = eval_zinm_g(a+da, p+dp, tab);
         }

         a = a+da;
         p = p+dp;

      }

      double pi = (nobs-z0) / (1-pow(p,a)) / nobs;
      if (pi > 1) pi = 1.0;
      if (pi < 0) pi = 0.0;
      double loglik = ll_zinm(a, p, pi, tab);
      if (loglik > max_loglik) {
         max_loglik = loglik;
         par->alpha = a;
         par->pi = pi;
         par->p[0] = p;
      }
            
   }

   free(tab);
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

int
mle_nm
(
   const uint32_t *x,
         uint32_t dim,
         uint32_t nobs,
         zinm_par_t *par
)
// SYNOPSIS:
//   Estimate the parameters of a negative multinomial of negative
//   binomial distribution with the Newton-Raphson method. The
//   negative binomial is the special case that 'dim' is equal to 1.
//
// PARAMETERS:
//   x: data (counts)
//   dim: the number of dimensions of the data
//   nobs: the size of the sample (number of observations)
//   par: a parameter set used as return value
//
// DETAILS:
//   The data must be passed as a succession of vectors of dimension
//   'dim'. This is to facilitate reading from a file where counts
//   are typically given in columns whereas the file is read by row.
//   For instance if 'dim' is 3 and the first observation is the
//   vector (1,2,4), the first elements of 'x' must be 1, 2, 4, ...
//
//   The caller is responsible for the memory allocation of the
//   parameters. A typical call would thus look as shown below.
//
//       zinm_par_t *par = new_zinm_par(dim);
//       if (par == NULL) { /* out of memory */ }
//       if (!mle_nm(x, dim, nobs, par)) { /* something wrong /* }
//
//   Upon success, 'par->alpha' contains the alpha parameter and the
//   other parameters are in 'par->p[0]', ... , 'par->p[dim+1]'.
//
//   See the test cases for simple examples of usage.
//
// RETURN:
//   Upon success '(int) 1', and upon failure '(int) 0'.
//
// SIDE EFFECTS:
//   None (but calls 'malloc()' to allocate temporary variables).
{

   tab_t *tab = NULL;
   double *means = NULL;

   int status = 0;

   // Estimate alpha.
   tab = tabulate(x, dim, nobs);
   if (tab == NULL) {
      fprintf(stderr, "error in function '%s()': %s:%d\n",
            __func__, __FILE__, __LINE__);
      goto clean_and_return;
   }

   double alpha = nb_est_alpha(tab);

   if (alpha < 0) {
      fprintf(stderr, "domain error in function '%s(): %s:%d'\n",
            __func__, __FILE__, __LINE__);
      goto clean_and_return;
   }

   // Compute the means in all dimensions.
   means = malloc(dim * sizeof(double));
   if (means == NULL) {
      fprintf(stderr, "memory error in function '%s()': %s:%d\n",
            __func__, __FILE__, __LINE__);
      goto clean_and_return;
   }

   compute_means(x, dim, nobs, means);

   double sum = 0;
   for (size_t i = 0 ; i < dim ; i++) sum += means[i];

   // Write the estimates to 'par'.
   par->alpha = alpha;
   par->p[0] = alpha / (alpha + sum);
   for (size_t i = 1 ; i < dim+1 ; i++) {
      par->p[i] = par->p[0] / alpha * means[i-1];
   }

   // Success.
   status = 1;

clean_and_return:
   if (means != NULL) free(means);
   if (tab != NULL) free(tab);

   return status;

}

// helper functions //

tab_t *
tabulate
(
   const uint32_t * x,
         uint32_t   dim,
         uint32_t   nobs
)
{

   tab_t *tab = NULL; // The return value.
   kv_t *counter = NULL;

   // Count the "row sums" of 'x'.
   counter = new_counter();

   if (counter == NULL) {
      fprintf(stderr, "error in function '%s()': %s:%d\n",
            __func__, __FILE__, __LINE__);
      goto clean_and_return;
   }

   for (size_t i = 0 ; i < nobs ; i++) {
      size_t sum = 0;
      for (size_t j = 0 ; j < dim ; j++) sum += x[i+j*nobs];
      if (!update(counter, sum)) {
         fprintf(stderr, "counting error in function '%s()': %s:%d\n",
               __func__, __FILE__, __LINE__);
         goto clean_and_return;
      }
   }

   // Convert counts to 'tab_t'.
   size_t sz = 0;
   for (size_t i = 0 ; i < HASH_SIZE ; i++) sz += (counter[i].val > 0);

   tab = malloc(sizeof(tab_t) + 2*sz * sizeof(int));
   if (tab == NULL) {
      fprintf(stderr, "memory error in function '%s()': %s:%d\n",
            __func__, __FILE__, __LINE__);
      goto clean_and_return;
   }

   tab->size = sz;
   tab->num = tab->val + sz;

   size_t pos = 0;
   for (size_t i = 0 ; i < HASH_SIZE ; i++) {
      if (counter[i].val > 0) {
         tab->val[pos] = counter[i].key;
         tab->num[pos] = counter[i].val;
         pos++;
      }
   }


clean_and_return:
   if (counter != NULL) free(counter);
   return tab;

}


void
compute_means
(
   const uint32_t *x,
         uint32_t dim,
         uint32_t nobs,
         double *means
)
{

   memset(means, 0, dim*sizeof(double));
   for (size_t i = 0 ; i < nobs ; i++) {
   for (size_t j = 0 ; j < dim ; j++) {
      means[j] += x[j+dim*i];
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
   int r
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




// Implementation of a counter of (small) 32 bit integers.
// The implementation is a hash table without hashing. Small numbers
// are stored at the beginning of the array, large numbers could
// be anywhere. Collisions (which happen between a small and a large
// number or between two large numbers) are resolved by quadratic
// probing. The size of the table must be a prime for the probing
// to scan every spot.
//
// A counter is just an array of entries (an 8 byte struct with
// the integer and its count). The counter can hold 'HASH_SIZE'
// entries. Once full, it cannot be incremented.
//
// This structure is a simple (dense) array for small integers and
// a sparse array for large inteters. The main advantage is that
// it will be fast and simple for small values, and that it will
// not choke on large values.

kv_t *
new_counter
(void)
// Total size of the table: 830 kB.
{

   kv_t *new = calloc(HASH_SIZE, sizeof(kv_t));

   if (new == NULL) {
      fprintf(stderr, "memory error in function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      return NULL;
   }

   return new;

}


int
update
(
    kv_t *table,
    uint32_t itg
)
// Return 1 if the operation succeeds, and 0 otherwise
// (i.e. when the table is full).
{

   for (size_t i = 0 ; i < HASH_SIZE ; i++) {

      // Quadratic probing.
      kv_t *e = table + ((itg + i*i) % HASH_SIZE);
      if (e->key == itg || e->val == 0) {
         // Found the integer or an empty spot.
         e->key = itg;
         e->val++;
         return 1;
      }

   }

   // The whole table has been scanned and there
   // is nowhere to store this integer. Return the
   // code for failure.
   
   return 0;

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
