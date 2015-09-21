#include "pso.h"

#define ITER 10000000
#define RAND (1+(0.5*(drand48()-.5)))

// Global lock for mutex.
pthread_mutex_t lock;

// Basic parameters.
typedef struct {
   double l;      // log-likelihood
   double a;      // alpha
   double p;
   double pi;
} par_t;

// Arguments of `particle`.
typedef struct {
         par_t *opt;
   const tab_t *tab;
         int *signal;
} parg_t;


// --                  memory functions                  -- //

par_t *
params_new
(void)
{
   par_t *new = calloc(1, sizeof(par_t));
   if (new == NULL) {
      fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }

   new->l = -1.0/0.0;

   return new;

}


par_t *
params_set
(
   par_t *par,
   double a,
   double p,
   double pi
)
{

   par->a = a;
   par->p = p;
   par->pi = pi;

   return par;

}


// --                  main functions                  -- //

double
loglik
(
         par_t *par,
   const tab_t *tab
)
{

   // Convenience variables.
   const unsigned int *val = tab->val;
   const unsigned int *num = tab->num;
   const unsigned int z0 = val[0] == 0 ? num[0] : 0;
   const double a = par->a;
   const double p = par->p;
   const double pi = par->pi;
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

   return par->l = retval;

}


par_t *
params_change
(
   par_t * restrict new,
   par_t * restrict old,
   int which
)
{

   switch (which) {
   case 0:
      new->a = gen_gamma(old->a, 1);
      if (new->a < 0.01) new->a = 0.01;
      break;
   case 1:
      //new->p = gen_beta(old->p, 1-old->p);
      //if (new->p > .99) new->p = .99;
      //if (new->p < .01) new->p = .01;
      new->p = drand48();
      break;
   case 2:
      //new->pi = gen_beta(old->pi, 1-old->pi);
      //if (new->pi > .99) new->pi = .99;
      //if (new->pi < .01) new->pi = .01;
      new->pi = drand48();
   }

   return new;

}


void *
particle
(void *arg)
{

   // Unpack arguments.
         int *signal = ((parg_t *)arg)->signal;
         par_t *opt = ((parg_t *)arg)->opt;
   const tab_t *tab = ((parg_t *)arg)->tab;

   par_t *par = params_new();
   par_t *try = params_new();
   
   if ((par == NULL) || (try == NULL)) return NULL;

   // Copy best conditions.
   memcpy(par, opt, sizeof(par_t));
   memcpy(try, opt, sizeof(par_t));

   // Start simulated annealing.

   double new = par->l;
   double total = 0.0;
   double T = 0.0;
   for (int iter = 0 ; iter < ITER ; iter++) {
      double old = new;
      if (*signal > 0) {
         // Copy best conditions.
         pthread_mutex_lock(&lock);
            memcpy(par, opt, sizeof(par_t));
            (*signal)--;
         pthread_mutex_unlock(&lock);
      }

      params_change(try, par, iter % 3);

      // Keep the global best hit. 
      new = loglik(try, tab);
      total += abs(new-old);
      pthread_mutex_lock(&lock);
      if (new > opt->l) {
         fprintf(stderr, "%f\n", new);
         *signal = 6;
         memcpy(opt, try, sizeof(par_t));
      }
      pthread_mutex_unlock(&lock);

      // Auto-adjust temperature every 25 cycles.
      if (iter % 25 == 0) {
         T = total / 25 / 0.69;
         if (T < 0.001) T = 0.001;
         total = 0.0;
      }
      else {
         T /= .95;
      }

      // Metropolis move.
      double delta = (new-old)/T;
      if ((delta > 0) || (drand48() < exp(delta))) {
         memcpy(par, try, sizeof(par_t));
      }

   }

   free(par);
   free(try);

   return NULL;

}


double
pso
(
         double a,
         double p,
         double pi,
   const tab_t *tab
)
{

   //outf = fopen("tracks.txt", "w");
   int n_threads = 12;

   par_t *opt = params_new();
   if (opt == NULL) {
      fprintf(stderr, "memory error: cannot allocate 'opt'");
      return 1.0/0.0;
   }

   params_set(opt, a, p, pi);

   // Initialize random generator.
   srand48(time(NULL));

   // Initialize mutex.
   int err = pthread_mutex_init(&lock, NULL);
   if (err) {
      fprintf(stderr, "error initializing mutex (%d)\n", err);
      return 1.0/0.0;
   }

   // Initialize arguments.
   int signal = 0;
   parg_t arg = {
      .opt = opt,
      .tab = tab,
      .signal = &signal,
   };

   // Allocate 'tid'.
   fprintf(stderr, "starting threads...\n");
   pthread_t *tid = (pthread_t *) malloc(n_threads * sizeof(pthread_t));
   for (int i = 0 ; i < n_threads ; i++) {
      err = pthread_create(tid+i, NULL, &particle, &arg);
      if (err) {
         fprintf(stderr, "error creating thread (%d)\n", err);
         return 1.0/0.0;
      }
   }

   // Wait for threads to return.
   for (int i = 0 ; i < n_threads ; i++) {
      pthread_join(tid[i], NULL);
   }

   pthread_mutex_destroy(&lock);
   fprintf(stderr, "%f (%f, %f, %f)\n",
         opt->l, opt->a, opt->p, opt->pi);
   double retval = opt->l;
   free(opt);
   free(tid);

   return retval;

}
