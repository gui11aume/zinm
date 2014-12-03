#include "zinm.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


double
dlda
(
         double a,
         double m,
         size_t n,
   const tab_t *tab
)
{

   double retval;
   double prev;
   // Convenience variables.
   const unsigned int *val = tab->val;
   const unsigned int *num = tab->num;

   prev = digamma(a + val[0]);
   retval = n*(log(a) - digamma(a) - log(a + m)) + num[0] * prev;
   // Iterate over the occurrences and compute the new value
   // of digamma either by the recurrence relation, or by
   // a new call to 'digamma()', whichever is faster.
   for (size_t i = 1 ; i < tab->size ; i++) {
      prev = (val[i] - val[i-1] == 1) ?
         prev + 1.0 / (a-1 + val[i]) :
         digamma(a + val[i]);
      retval += num[i] * prev;
   }

   return retval;

}


double
d2lda2
(
   double a,
   double m,
   size_t n,
   const tab_t *tab
)
{

   double retval;
   double prev;
   // Convenience variables.
   const unsigned int *val = tab->val;
   const unsigned int *num = tab->num;

   prev = trigamma(a + val[0]);
   retval = n*(m/(a*(a+m)) - trigamma(a)) + num[0] * prev;
   // Iterate over the occurrences and compute the new value
   // of trigamma either by the recurrence relation, or by
   // a new call to 'trigamma()', whichever is faster.
   for (size_t i = 1 ; i < tab->size ; i++) {
      prev = (val[i] - val[i-1] == 1) ?
         prev - 1.0 / ((a-1 +val[i])*(a-1 +val[i])) :
         trigamma(a + val[i]);
      retval += num[i] * prev;
   }

   return retval;

}


nm_par_t *
mle_nm
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

   // Compute the means in all dimensions.
   compute_means(x, dim, nobs, means);

   // Compute the mean of marginal sums.
   double a = 1.0;
   double mean = 0.0;
   for (size_t i = 0 ; i < dim ; i++) mean += means[i];

   // Find upper and lower bouds for a(lpha).
   double a_lo;
   double a_hi;
   if (dlda(a, mean, nobs, tab) < 0) {
      a /= 2;
      while (dlda(a, mean, nobs, tab) < 0) a /= 2;
      a_lo = a;
      a_hi = a*2;
   }
   else {
      a *= 2;
      while (dlda(a, mean, nobs, tab) > 0) a *= 2;
      a_lo = a/2;
      a_hi = a;
   }

   // Input is pathological.
   if (a_lo > 128) {
      free(means);
      return NULL;
   }

   double new_a = (a_lo + a_hi) / 2;
   for (int i = 0 ; i < ZINM_NR_MAXITER ; i++) {
      a = (new_a < a_lo || new_a > a_hi) ?
         (a_lo + a_hi) / 2 :
         new_a;
      double eval = dlda(a, mean, nobs, tab);
      if (eval < 0) a_hi = a; else a_lo = a;
      if ((a_hi - a_lo) < ZINM_NR_TOLERANCE) break;
      new_a = a - eval / d2lda2(a, mean, nobs, tab);
   }

   free(tab);

   nm_par_t *par = new_nm_par(dim);
   if (par == NULL) {
      free(means);
      fprintf(stderr, "memory error: %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }
   par->alpha = a;
   par->p[0] = a / (a + mean);
   for (size_t i = 1 ; i < dim+1 ; i++) {
      par->p[i] = par->p[0] / a * means[i-1];
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


nm_par_t *
new_nm_par
(
   size_t r
)
{

   nm_par_t *new = calloc(1, sizeof(nm_par_t) + (r+1)*sizeof(double));
   if (new == NULL) {
      fprintf(stderr, "memory error: %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }
   new->r = r;

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
