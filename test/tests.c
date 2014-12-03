#include "zinm.h"
#include "unittest.h"

void
test_new_histo
(void)
{

   histo_t *histo;
   histo = new_histo();
   test_assert_critical(histo != NULL);
   for (size_t i = 0 ; i < HISTO_INIT_SIZE ; i++) {
      test_assert(histo->num[i] == 0);
   }
   free(histo);

   redirect_stderr();
   set_alloc_failure_rate_to(1.0);
   histo = new_histo();
   reset_alloc();
   unredirect_stderr();
   test_assert(histo == NULL);
   test_assert(strcmp(caught_in_stderr(), "") != 0);

   return;

}

void
test_histo_push
(void)
{

   histo_t *histo;
   histo = new_histo();
   test_assert_critical(histo != NULL);
   for (size_t i = 0 ; i < (1 << 16) ; i++) {
      test_assert(histo_push(&histo, i) == 0);
   }
   test_assert(histo->size == (1 << 16));
   for (size_t i = 0 ; i < (1 << 16) ; i++) {
      test_assert(histo->num[i] == 1);
   }
   free(histo);

   histo = new_histo();
   test_assert_critical(histo != NULL);
   for (size_t i = (1 << 16) ; i > 0 ; i--) {
      test_assert(histo_push(&histo, i-1) == 0);
   }
   test_assert(histo->size == 2*((1 << 16)-1));
   for (size_t i = 0 ; i < (1 << 16) ; i++) {
      test_assert(histo->num[i] == 1);
   }
   free(histo);

   histo = new_histo();
   test_assert_critical(histo != NULL);
   redirect_stderr();
   set_alloc_failure_rate_to(1.0);
   test_assert(histo_push(&histo, HISTO_INIT_SIZE) == 1);
   reset_alloc();
   unredirect_stderr();
   test_assert(strcmp(caught_in_stderr(), "") != 0);
   free(histo);

   return;

}

void
test_compress_histo
(void)
{

   histo_t *histo = new_histo();
   test_assert_critical(histo != NULL);
   for (size_t i = 0 ; i < 4096 ; i += 2) {
      test_assert(histo_push(&histo, i) == 0);
   }

   tab_t *tab;
   tab = compress_histo(histo);
   test_assert_critical(tab != NULL);
   test_assert(tab->size == 2048);
   for (size_t i = 0 ; i < tab->size ; i++) {
      test_assert(tab->val[i] == 2*i);
      test_assert(tab->num[i] == 1);
   }

   free(tab);

   redirect_stderr();
   set_alloc_failure_rate_to(1.0);
   tab = compress_histo(histo);
   reset_alloc();
   unredirect_stderr();
   test_assert(tab == NULL);
   test_assert(strcmp(caught_in_stderr(), "") != 0);

   free(histo);
   return;

}

void
test_tabulate
(void)
{

   size_t x[] = {1,2,3,4,5,6,7,8,9,10,11,12};
   tab_t *tab;
   tab = tabulate(x, 3, 4);
   test_assert_critical(tab != NULL);
   test_assert(tab->size == 4);
   test_assert(tab->val[0] == 15);
   test_assert(tab->val[1] == 18);
   test_assert(tab->val[2] == 21);
   test_assert(tab->val[3] == 24);
   for (int i = 0 ; i < 4 ; i++) {
      test_assert(tab->num[i] == 1);
   }
   free(tab);

   tab = tabulate(x, 4, 3);
   test_assert_critical(tab != NULL);
   test_assert(tab->size == 3);
   test_assert(tab->val[0] == 22);
   test_assert(tab->val[1] == 26);
   test_assert(tab->val[2] == 30);
   for (int i = 0 ; i < 3 ; i++) {
      test_assert(tab->num[i] == 1);
   }
   free(tab);

   size_t y[] = {4096,2048,1024,512};
   tab = tabulate(y, 1, 4);
   test_assert_critical(tab != NULL);
   test_assert(tab->size == 4);
   test_assert(tab->val[0] == 512);
   test_assert(tab->val[1] == 1024);
   test_assert(tab->val[2] == 2048);
   test_assert(tab->val[3] == 4096);
   for (int i = 0 ; i < 4 ; i++) {
      test_assert(tab->num[i] == 1);
   }
   free(tab);

   tab = tabulate(y, 2, 2);
   test_assert_critical(tab != NULL);
   test_assert(tab->size == 2);
   test_assert(tab->val[0] == 2560);
   test_assert(tab->val[1] == 5120);
   for (int i = 0 ; i < 2 ; i++) {
      test_assert(tab->num[i] == 1);
   }
   free(tab);

   return;

}


void
test_compute_means
(void)
{

   double mean;

   size_t x1[1] = {2595};
   compute_means(x1, 1, 1, &mean);
   test_assert(mean == 2595.0);

   size_t x2[7] = {0,0,0,1,1,2,5};
   compute_means(x2, 1, 7, &mean);
   test_assert(fabs(mean-1.28571428) < 1e-6);

   size_t x3[5] = {0,89,231,55,309};
   compute_means(x3, 1, 5, &mean);
   test_assert(fabs(mean-136.8) < 1e-6);

   // 0:112, 1:94, 2:28, 3:12, 4:3, 7:1
   size_t x4[250] = {0,0,0,3,0,0,1,1,1,1,1,2,0,2,0,0,1,0,0,0,1,1,0,1,
      1,0,1,1,0,2,1,0,2,1,1,0,2,1,1,1,1,1,0,0,2,0,2,1,1,1,2,1,0,0,
      1,0,1,0,0,1,0,0,3,2,0,0,0,0,0,2,1,1,1,0,0,1,0,0,1,0,0,1,0,1,
      0,1,2,1,2,1,0,0,0,2,0,0,0,1,2,1,0,1,1,1,2,0,0,0,0,0,2,1,3,0,
      2,3,0,0,0,1,0,0,1,1,0,0,0,0,1,0,0,1,0,3,1,0,0,0,1,1,0,0,0,0,
      0,1,0,0,0,0,1,2,1,0,2,4,0,1,0,1,0,1,0,1,1,1,0,1,1,1,1,0,0,1,
      0,1,1,3,1,1,1,1,0,0,0,0,3,1,3,0,1,1,0,0,0,1,1,0,1,2,4,2,0,0,
      4,0,2,1,0,0,2,1,2,1,7,1,2,3,0,0,1,1,0,3,1,1,1,3,1,1,0,0,0,0,
      1,2,2,0,1,1,0,1,1,1,0,2,3,0,0,0};
   compute_means(x4, 1, 250, &mean);
   test_assert(fabs(mean-0.82) < 1e-6);

   return;

}


void
test_dlda
(void)
{

   // 0:14, 1:5, 2:4, 3:1, 5:1
   size_t x[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                     1,1,1,1,1,2,2,2,2,3,5 };
   double mean;
   compute_means(x, 1, 25, &mean);
   tab_t *tab = tabulate(x, 1, 25);
   test_assert(fabs(dlda(1.0, mean, 25, tab)+0.12747262) < 1e-6);
   test_assert(fabs(dlda(1.1, mean, 25, tab)+0.24215981) < 1e-6);
   test_assert(fabs(dlda(1.2, mean, 25, tab)+0.31636395) < 1e-6);
   test_assert(fabs(dlda(1.3, mean, 25, tab)+0.36350700) < 1e-6);
   test_assert(fabs(dlda(1.4, mean, 25, tab)+0.39225466) < 1e-6);
   test_assert(fabs(dlda(1.5, mean, 25, tab)+0.40834322) < 1e-6);
   test_assert(fabs(dlda(2.0, mean, 25, tab)+0.39975512) < 1e-6);

   free(tab);

   return;

}


void
test_d2lda2
(void)
{

   // 0:14, 1:5, 2:4, 3:1, 5:1
   size_t x[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                    1,1,1,1,1,2,2,2,2,3,5 };
   tab_t *tab = tabulate(x, 1, 25);
   double mean;
   compute_means(x, 1, 25, &mean);
   test_assert(fabs(d2lda2(1.0, mean, 25, tab)+1.41167874) < 1e-6);
   test_assert(fabs(d2lda2(1.1, mean, 25, tab)+0.91683021) < 1e-6);
   test_assert(fabs(d2lda2(1.2, mean, 25, tab)+0.58911102) < 1e-6);
   test_assert(fabs(d2lda2(1.3, mean, 25, tab)+0.36790287) < 1e-6);
   test_assert(fabs(d2lda2(1.4, mean, 25, tab)+0.21643981) < 1e-6);
   test_assert(fabs(d2lda2(1.5, mean, 25, tab)+0.11168877) < 1e-6);
   test_assert(fabs(d2lda2(2.0, mean, 25, tab)-0.08773865) < 1e-6);

   free(tab);

   return;

}


void
test_mle_nm
(void)
{

   // These test cases have been verified with R.
   nm_par_t *par;
   
   // 0:14, 1:5, 2:4, 3:1, 5:1
   size_t x1[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  1,1,1,1,1,2,2,2,2,3,5 };
   par = mle_nm(x1, 1, 25);
   test_assert_critical(par != NULL);
   test_assert(fabs(par->alpha-0.9237) < 1e-3);
   test_assert(fabs(par->p[0]-0.5237) < 1e-3);
   test_assert(fabs(par->p[1]-0.4763) < 1e-3);
   free(par);

   // 0:27, 1:12, 2:8, 3:1, 4:1, 5:1
   size_t x2[50] = {3,0,1,2,0,0,1,0,0,0,0,1,1,0,0,1,2,2,0,0,0,1,2,
      0, 0,0,0,0,4,0,0,0,1,5,1,0,1,2,1,2,2,2,0,0,0,1,0,1,0,0};
   par = mle_nm(x2, 1, 50);
   test_assert_critical(par != NULL);
   test_assert(fabs(par->alpha-1.3436) < 1e-3);
   test_assert(fabs(par->p[0]-0.6267) < 1e-3);
   test_assert(fabs(par->p[1]-0.3732) < 1e-3);
   free(par);

   // 0:12, 1:7, 2:13, 3:4, 4:6, 5:2, 6:1, 7:3, 8:1, 9:1
   size_t x3[50] = {4,5,2,1,2,4,2,2,0,4,2,1,3,6,0,0,7,3,0,8,4,2,0,
      0,2,3,2,3,7,9,2,4,0,4,2,0,0,2,5,1,1,2,1,0,0,0,1,2,1,7};
   par = mle_nm(x3, 1, 50);
   test_assert_critical(par != NULL);
   test_assert(fabs(par->alpha-1.7969) < 1e-3);
   test_assert(fabs(par->p[0]-0.4221) < 1e-3);
   test_assert(fabs(par->p[1]-0.5779) < 1e-3);
   free(par);

   // 0:39, 1:8, 2:2, 3:1
   size_t x4[50] = {1,0,0,0,0,0,3,1,0,1,0,0,0,0,0,0,0,1,0,0,0,2,0,
      2,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0};
   par = mle_nm(x4, 1, 50);
   test_assert_critical(par != NULL);
   test_assert(fabs(par->alpha-0.7073) < 1e-3);
   test_assert(fabs(par->p[0]-0.7021) < 1e-3);
   test_assert(fabs(par->p[1]-0.2978) < 1e-3);
   free(par);

   // 0:59, 1:83, 2:99, 3:67, 4:67, 5:49, 6:27, 7:22, 8:11, 9:6
   // 10:6, 11:3, 12:2, 13:3
   size_t x5[500] = {1,0,0,1,1,2,1,0,5,7,1,3,3,1,6,0,2,5,7,0,5,2,1,
       10,5,3,4,5,7,0,8,6,3,0,2,1,1,0,2,3,7,2,3,2,2,1,0,4,4,2,4,2,
       0,6,3,2,5,2,1,4,3,4,2,2,5,3,2,0,2,8,1,3,1,7,5,1,4,1,1,0,2,
       2,4,1,1,1,4,1,3,4,4,10,5,2,0,7,1,6,1,3,6,4,0,2,4,1,12,2,5,
       6,5,4,1,11,0,1,3,2,4,2,0,2,3,4,0,2,9,9,7,4,2,1,3,3,3,4,2,9,
       2,4,3,2,2,4,2,5,3,0,1,3,2,0,3,3,4,1,3,3,5,7,3,3,2,1,5,5,4,
       6,1,1,1,2,9,5,1,2,4,0,2,1,0,3,2,4,3,1,4,2,1,4,1,6,0,6,5,3,
       5,2,0,1,2,1,0,5,3,2,7,6,4,3,2,5,7,5,5,1,1,3,10,2,0,5,0,1,2,
       0,5,1,2,3,6,4,0,3,1,2,2,4,3,0,3,2,5,4,10,1,2,4,4,2,13,4,3,
       1,5,4,8,5,6,2,3,4,3,1,5,5,1,8,2,0,5,7,3,2,2,4,2,3,1,5,3,7,
       13,1,4,7,5,5,0,3,0,4,2,3,1,2,4,2,8,1,2,5,6,1,1,0,7,2,2,3,5,
       12,2,2,2,0,3,3,4,0,2,5,1,10,0,7,6,5,0,11,2,3,7,3,5,4,2,1,2,
       4,0,2,2,2,0,6,2,3,4,2,3,7,3,5,2,5,0,4,4,6,3,1,2,7,3,0,2,5,
       7,2,2,0,0,0,6,3,0,1,1,5,5,2,6,2,4,6,0,1,2,3,2,2,2,3,4,1,1,
       4,0,2,0,1,3,4,1,2,2,3,1,4,4,3,4,4,1,5,2,13,4,10,5,6,1,0,5,
       0,0,5,6,0,1,8,5,1,3,1,8,1,8,1,6,7,2,8,2,2,3,3,0,4,2,1,9,6,
       0,6,7,1,8,2,2,1,11,3,0,4,2,5,1,6,8,3,4,7,0,4,2,4,1,1,1,6,0,
       4,4,6,2,1,3,1,0,4,9,3,1,4,2,2,0,1};
   par = mle_nm(x5, 1, 500);
   test_assert_critical(par != NULL);
   test_assert(fabs(par->alpha-3.0057) < 1e-3);
   test_assert(fabs(par->p[0]-0.4854) < 1e-3);
   test_assert(fabs(par->p[1]-0.5145) < 1e-3);
   free(par);

   return;

}


int
main(
   int argc,
   char **argv
)
{

   // Register test cases //
   const static test_case_t test_cases[] = {
      {"new_histo", test_new_histo},
      {"histo_push", test_histo_push},
      {"compress_histo", test_compress_histo},
      {"tabulate", test_tabulate},
      {"compute_means", test_compute_means},
      {"dlda", test_dlda},
      {"d2lda2", test_d2lda2},
      {"mle_nm", test_mle_nm},
      {NULL, NULL}
   };

   return run_unittest(argc, argv, test_cases);

}
