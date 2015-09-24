#ifndef _ZINM_HEADER
#define _ZINM_HEADER

struct zinm_part_t;
typedef struct zinm_par_t zinm_par_t;

struct zinm_par_t {
   size_t r;
   double pi;
   double alpha;
   double p[];
};

// Exported functions.
int          mle_nm(const uint32_t *, uint32_t, uint32_t, zinm_par_t *par);
zinm_par_t * mle_zinm(uint32_t *, uint32_t, uint32_t, zinm_par_t *par);
zinm_par_t * new_zinm_par(int);

#endif
