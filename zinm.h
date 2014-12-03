#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef _PUBLIC_ZINM_HEADER
#define _PUBLIC_ZINM_HEADER

#define HISTO_INIT_SIZE 128
#define ZINM_NR_MAXITER 32
#define ZINM_NR_TOLERANCE 1e-6

struct histo_t;
struct tab_t;
typedef struct histo_t histo_t;
typedef struct tab_t tab_t;

struct histo_t {
   size_t size;
   unsigned int num[];
};

struct tab_t {
   size_t size;
   unsigned int *num;
   unsigned int val[];
};

tab_t   * compress_histo(histo_t *);
double    compute_mean(tab_t *);
double    d2lda2(double, double, size_t, const tab_t *);
double    dlda(double, double, size_t, const tab_t *);
double    digamma(double);
int       histo_push(histo_t **, size_t);
double	  mle_nm(size_t *, size_t, size_t);
histo_t * new_histo(void);
double    trigamma(double);
tab_t   * tabulate(size_t *, size_t, size_t);

#endif
