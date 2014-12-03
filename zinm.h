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
struct nm_part_t;
struct tab_t;

typedef struct histo_t histo_t;
typedef struct nm_par_t nm_par_t;
typedef struct tab_t tab_t;

struct histo_t {
   size_t size;
   unsigned int num[];
};

struct nm_par_t {
   size_t r;
   double alpha;
   double p[];
};

struct tab_t {
   size_t size;
   unsigned int *num;
   unsigned int val[];
};

tab_t    * compress_histo(histo_t *);
void       compute_means(size_t *, size_t, size_t, double *);
double     d2lda2(double, double, size_t, const tab_t *);
double     dlda(double, double, size_t, const tab_t *);
double     digamma(double);
int        histo_push(histo_t **, size_t);
nm_par_t * mle_nm(size_t *, size_t, size_t);
histo_t  * new_histo(void);
nm_par_t * new_nm_par(size_t);
double     trigamma(double);
tab_t    * tabulate(size_t *, size_t, size_t);

#endif
