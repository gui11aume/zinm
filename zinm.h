#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef _ZINM_HEADER
#define _ZINM_HEADER

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

// Exported functions.
double	    ll_zinm(double, double, double, const tab_t *);
int          mle_nm(size_t *, size_t, size_t, zinm_par_t *par);
zinm_par_t * mle_zinm(size_t *, size_t, size_t, zinm_par_t *par);
double       nb_est_alpha(tab_t *);
histo_t    * new_histo(void);
zinm_par_t * new_zinm_par(size_t);
tab_t      * tabulate(size_t *, size_t, size_t);

#endif
