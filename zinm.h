#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef _PUBLIC_ZINM_HEADER
#define _PUBLIC_ZINM_HEADER

#define HISTO_INIT_SIZE 128
#define ZINM_MAXITER 32
#define ZINM_TOL 1e-6

#define sq(x) ((x)*(x))

struct histo_t;
struct zinm_part_t;
struct tab_t;

typedef struct histo_t histo_t;
typedef struct tab_t tab_t;
typedef struct zinm_par_t zinm_par_t;

struct histo_t {
   size_t size;
   unsigned int num[];
};

struct zinm_par_t {
   size_t r;
   double pi;
   double alpha;
   double p[];
};

struct tab_t {
   size_t size;
   unsigned int *num;
   unsigned int val[];
};

tab_t      * compress_histo(histo_t *);
void         compute_means(size_t *, size_t, size_t, double *);
double       digamma(double);
double       eval_nb_dfda(double, const tab_t *);
double       eval_nb_f(double, const tab_t *);
double 	     eval_zinm_dfda(double, double, unsigned int);
double 	     eval_zinm_dfdp(double, double, unsigned int, double);
double 	     eval_zinm_dgda(double, double, const tab_t *);
double 	     eval_zinm_f(double, double, unsigned int, double);
double 	     eval_zinm_g(double, double, const tab_t *);
int          histo_push(histo_t **, size_t);
double	     ll_zinm(double, double, double, const tab_t *);
zinm_par_t * mle_nm(size_t *, size_t, size_t);
zinm_par_t * mle_zinm(size_t *, size_t, size_t);
double       nb_est_alpha(tab_t *);
histo_t    * new_histo(void);
zinm_par_t * new_zinm_par(size_t);
double       trigamma(double);
tab_t      * tabulate(size_t *, size_t, size_t);

#endif
