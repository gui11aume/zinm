#include "unittest.h"
#include "zinm.c"


void
test_tabulate
(void)
{

   uint32_t x[] = {1,2,3,4,5,6,7,8,9,10,11,12};
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

   uint32_t y[] = {4096,2048,1024,512};
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

   uint32_t z[] = {1,2,23784983};
   tab = tabulate(z, 1, 3);
   test_assert_critical(tab != NULL);
   test_assert(tab->size == 3);
   test_assert(tab->val[0] == 1);
   test_assert(tab->val[1] == 2);
   test_assert(tab->val[2] == 23784983);

   free(tab);
   return;

}


void
test_compute_means
(void)
{

   double mean;

   uint32_t x1[1] = {2595};
   compute_means(x1, 1, 1, &mean);
   test_assert(mean == 2595.0);

   uint32_t x2[7] = {0,0,0,1,1,2,5};
   compute_means(x2, 1, 7, &mean);
   test_assert(fabs(mean-1.28571428) < 1e-6);

   uint32_t x3[5] = {0,89,231,55,309};
   compute_means(x3, 1, 5, &mean);
   test_assert(fabs(mean-136.8) < 1e-6);

   // 0:112, 1:94, 2:28, 3:12, 4:3, 7:1
   uint32_t x4[250] = {0,0,0,3,0,0,1,1,1,1,1,2,0,2,0,0,1,0,0,0,1,1,0,1,
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

   // Compute means in several dimensions.
   double means2[2];

   uint32_t x5[2] = {1,2595};
   compute_means(x5, 2, 1, means2);
   test_assert(means2[0] == 1.0);
   test_assert(means2[1] == 2595.0);

   uint32_t x6[4] = {1,2,11,12};
   compute_means(x6, 2, 2, means2);
   test_assert(means2[0] == 6.0);
   test_assert(means2[1] == 7.0);

   return;

}


void
test_eval_nb_f
(void)
{

   // 0:14, 1:5, 2:4, 3:1, 5:1
   uint32_t x[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                     1,1,1,1,1,2,2,2,2,3,5 };
   tab_t *tab = tabulate(x, 1, 25);
   test_assert(fabs(eval_nb_f(1.0, tab)+0.12747262) < 1e-6);
   test_assert(fabs(eval_nb_f(1.1, tab)+0.24215981) < 1e-6);
   test_assert(fabs(eval_nb_f(1.2, tab)+0.31636395) < 1e-6);
   test_assert(fabs(eval_nb_f(1.3, tab)+0.36350700) < 1e-6);
   test_assert(fabs(eval_nb_f(1.4, tab)+0.39225466) < 1e-6);
   test_assert(fabs(eval_nb_f(1.5, tab)+0.40834322) < 1e-6);
   test_assert(fabs(eval_nb_f(2.0, tab)+0.39975512) < 1e-6);

   free(tab);

   return;

}


void
test_eval_nb_dfda
(void)
{

   // 0:14, 1:5, 2:4, 3:1, 5:1
   uint32_t x[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                    1,1,1,1,1,2,2,2,2,3,5 };
   tab_t *tab = tabulate(x, 1, 25);
   double mean;
   compute_means(x, 1, 25, &mean);
   test_assert(fabs(eval_nb_dfda(1.0, tab)+1.41167874) < 1e-6);
   test_assert(fabs(eval_nb_dfda(1.2, tab)+0.58911102) < 1e-6);
   test_assert(fabs(eval_nb_dfda(1.3, tab)+0.36790287) < 1e-6);
   test_assert(fabs(eval_nb_dfda(1.4, tab)+0.21643981) < 1e-6);
   test_assert(fabs(eval_nb_dfda(1.5, tab)+0.11168877) < 1e-6);
   test_assert(fabs(eval_nb_dfda(2.0, tab)-0.08773865) < 1e-6);

   free(tab);

   return;

}


void
test_eval_zinm_f
(void)
{

   test_assert(fabs(eval_zinm_f(1.0, .5, 1, 2.0)) < 1e-6);
   test_assert(fabs(eval_zinm_f(1.0, .5, 2, 4.0)) < 1e-6);
   test_assert(fabs(eval_zinm_f(1.3, .7, 141, 8.3)-678.0838351) < 1e-6);

   return;

}


void
test_eval_zinm_g
(void)
{

   tab_t *tab;
   // 0:14, 1:5, 2:4, 3:1, 5:1
   uint32_t x1[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  1,1,1,1,1,2,2,2,2,3,5 };
   tab = tabulate(x1, 1, 25);
   test_assert_critical(tab != NULL);
   test_assert(fabs(eval_zinm_g(1, .5, tab)+0.1325713) < 1e-6);
   free(tab);

   return;

}


void
test_eval_zinm_dfda
(void)
{

   test_assert(fabs(eval_zinm_dfda(1, .5, 1)-1.2274112) < 1e-6);
   test_assert(fabs(eval_zinm_dfda(1, .5, 9)-11.0467015) < 1e-6);
   test_assert(fabs(eval_zinm_dfda(2, .5, 9)-12.9096451) < 1e-6);
   test_assert(fabs(eval_zinm_dfda(2, .3, 9)-25.1159846) < 1e-6);
   test_assert(fabs(eval_zinm_dfda(2.4, .3, 9)-26.3620864) < 1e-6);

   return;

}


void
test_eval_zinm_dfdp
(void)
{

   test_assert(eval_zinm_dfdp(1, .5, 0, 0) == 0.0);
   test_assert(eval_zinm_dfdp(1, .5, 1, 0) == 0.0);
   test_assert(fabs(eval_zinm_dfdp(2, .5, 1, 0)+3.5555555) < 1e-6);
   test_assert(fabs(eval_zinm_dfdp(2, .3, 1, 0)+19.5896899) < 1e-6);
   test_assert(fabs(eval_zinm_dfdp(2, .3, 9, 0)+176.3072092) < 1e-6);
   test_assert(fabs(eval_zinm_dfdp(2, .3, 9, 1.7)+179.7765970) < 1e-6);

   return;

}


void
test_eval_zinm_dgda
(void)
{

   tab_t *tab;
   // 0:14, 1:5, 2:4, 3:1, 5:1
   uint32_t x1[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  1,1,1,1,1,2,2,2,2,3,5 };
   tab = tabulate(x1, 1, 25);
   test_assert_critical(tab != NULL);
   test_assert(fabs(eval_zinm_dgda(1, .5, tab)+2.2547559) < 1e-6);
   test_assert(fabs(eval_zinm_dgda(2, .5, tab)+1.2605630) < 1e-6);
   test_assert(fabs(eval_zinm_dgda(2, .3, tab)+1.8764955) < 1e-6);
   free(tab);

   return;

}


void
test_ll_zinm
(void)
{

   tab_t *tab;
   // 0:14, 1:5, 2:4, 3:1, 5:1
   uint32_t x1[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  1,1,1,1,1,2,2,2,2,3,5 };
   tab = tabulate(x1, 1, 25);
   test_assert_critical(tab != NULL);
   test_assert(fabs(ll_zinm(1, .5, 1, tab)+22.5329303) < 1e-6);
   test_assert(fabs(ll_zinm(1, .5, .7, tab)+22.7832550) < 1e-6);
   test_assert(fabs(ll_zinm(2, .5, .7, tab)+23.7608409) < 1e-6);
   test_assert(fabs(ll_zinm(2, .3, .7, tab)+31.6978553) < 1e-6);
   free(tab);

   return;

}


void
test_nb_est_alpha
(void)
{

   tab_t *tab;
   // 0:14, 1:5, 2:4, 3:1, 5:1
   uint32_t x1[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  1,1,1,1,1,2,2,2,2,3,5 };
   tab = tabulate(x1, 1, 25);
   test_assert_critical(tab != NULL);
   test_assert(fabs(nb_est_alpha(tab)-0.9237) < 1e-3);
   free(tab);

   // 0:27, 1:12, 2:8, 3:1, 4:1, 5:1
   uint32_t x2[50] = {3,0,1,2,0,0,1,0,0,0,0,1,1,0,0,1,2,2,0,0,0,1,2,
      0, 0,0,0,0,4,0,0,0,1,5,1,0,1,2,1,2,2,2,0,0,0,1,0,1,0,0};
   tab = tabulate(x2, 1, 50);
   test_assert_critical(tab != NULL);
   test_assert(fabs(nb_est_alpha(tab)-1.3436) < 1e-3);
   free(tab);

   // 0:12, 1:7, 2:13, 3:4, 4:6, 5:2, 6:1, 7:3, 8:1, 9:1
   uint32_t x3[50] = {4,5,2,1,2,4,2,2,0,4,2,1,3,6,0,0,7,3,0,8,4,2,0,
      0,2,3,2,3,7,9,2,4,0,4,2,0,0,2,5,1,1,2,1,0,0,0,1,2,1,7};
   tab = tabulate(x3, 1, 50);
   test_assert_critical(tab != NULL);
   test_assert(fabs(nb_est_alpha(tab)-1.7969) < 1e-3);
   free(tab);

   // 0:39, 1:8, 2:2, 3:1
   uint32_t x4[50] = {1,0,0,0,0,0,3,1,0,1,0,0,0,0,0,0,0,1,0,0,0,2,0,
      2,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0};
   tab = tabulate(x4, 1, 50);
   test_assert_critical(tab != NULL);
   test_assert(fabs(nb_est_alpha(tab)-0.7073) < 1e-3);
   free(tab);

   // 0:59, 1:83, 2:99, 3:67, 4:67, 5:49, 6:27, 7:22, 8:11, 9:6
   // 10:6, 11:3, 12:2, 13:3
   uint32_t x5[500] = {1,0,0,1,1,2,1,0,5,7,1,3,3,1,6,0,2,5,7,0,5,2,1,
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
   tab = tabulate(x5, 1, 500);
   test_assert_critical(tab != NULL);
   test_assert(fabs(nb_est_alpha(tab)-3.0057) < 1e-3);

   free(tab);

}


void
test_mle_nm
(void)
{

   // These test cases have been verified with R.
   zinm_par_t *par = new_zinm_par(1);
   if (par == NULL) {
      fprintf(stderr, "test error line %d\n", __LINE__);
      return;
   }
   
   // 0:14, 1:5, 2:4, 3:1, 5:1
   uint32_t x1[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  1,1,1,1,1,2,2,2,2,3,5 };
   test_assert(mle_nm(x1, 1, 25, par));
   test_assert(fabs(par->alpha-0.9237) < 1e-3);
   test_assert(fabs(par->p[0]-0.5237) < 1e-3);
   test_assert(fabs(par->p[1]-0.4763) < 1e-3);

   // 0:27, 1:12, 2:8, 3:1, 4:1, 5:1
   uint32_t x2[50] = {3,0,1,2,0,0,1,0,0,0,0,1,1,0,0,1,2,2,0,0,0,1,2,
      0, 0,0,0,0,4,0,0,0,1,5,1,0,1,2,1,2,2,2,0,0,0,1,0,1,0,0};
   test_assert(mle_nm(x2, 1, 50, par));
   test_assert(fabs(par->alpha-1.3436) < 1e-3);
   test_assert(fabs(par->p[0]-0.6267) < 1e-3);
   test_assert(fabs(par->p[1]-0.3732) < 1e-3);

   // 0:12, 1:7, 2:13, 3:4, 4:6, 5:2, 6:1, 7:3, 8:1, 9:1
   uint32_t x3[50] = {4,5,2,1,2,4,2,2,0,4,2,1,3,6,0,0,7,3,0,8,4,2,0,
      0,2,3,2,3,7,9,2,4,0,4,2,0,0,2,5,1,1,2,1,0,0,0,1,2,1,7};
   test_assert(mle_nm(x3, 1, 50, par));
   test_assert(fabs(par->alpha-1.7969) < 1e-3);
   test_assert(fabs(par->p[0]-0.4221) < 1e-3);
   test_assert(fabs(par->p[1]-0.5779) < 1e-3);

   // 0:39, 1:8, 2:2, 3:1
   uint32_t x4[50] = {1,0,0,0,0,0,3,1,0,1,0,0,0,0,0,0,0,1,0,0,0,2,0,
      2,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0};
   test_assert(mle_nm(x4, 1, 50, par));
   test_assert(fabs(par->alpha-0.7073) < 1e-3);
   test_assert(fabs(par->p[0]-0.7021) < 1e-3);
   test_assert(fabs(par->p[1]-0.2978) < 1e-3);

   // 0:59, 1:83, 2:99, 3:67, 4:67, 5:49, 6:27, 7:22, 8:11, 9:6
   // 10:6, 11:3, 12:2, 13:3
   uint32_t x5[500] = {1,0,0,1,1,2,1,0,5,7,1,3,3,1,6,0,2,5,7,0,5,2,1,
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
   test_assert(mle_nm(x5, 1, 500, par));
   test_assert(fabs(par->alpha-3.0057) < 1e-3);
   test_assert(fabs(par->p[0]-0.4854) < 1e-3);
   test_assert(fabs(par->p[1]-0.5145) < 1e-3);

   free(par);
   par = NULL;

   // Run cases in two dimensions.
   par = new_zinm_par(2);
   if (par == NULL) {
      fprintf(stderr, "test error line %d\n", __LINE__);
      return;
   }

   // Data obtained with the R code shown below.
   //     set.seed(123)
   //     l = rgamma(200, shape=2.1)
   //     x = rpois(200, lambda=l)
   //     y = rpois(200, lambda=3.5*l)
   //     as.vector(t(cbind(x,y)))
   //
   // The data is pasted below.
   uint32_t x6[400] = {1,2,3,15,0,0,2,0,3,25,3,4,0,2,0,0,2,12,3,4,
      5,9,0,10,2,2,5,10,2,10,1,5,0,4,1,7,1,0,1,3,0,0,1,1,0,2,1,4,
      0,5,0,3,5,7,4,8,1,10,2,6,2,11,0,3,2,6,1,7,6,8,0,1,4,8,1,4,
      1,4,3,6,2,1,1,8,0,3,0,1,0,1,4,2,3,6,2,6,5,25,1,0,1,7,1,10,
      2,5,7,6,6,15,1,1,0,0,3,2,1,15,0,5,1,1,1,4,1,7,4,3,0,5,0,2,
      3,10,2,14,1,6,0,2,3,7,0,4,3,8,0,2,3,13,1,5,2,25,0,0,6,8,5,
      15,3,7,6,12,2,4,4,4,0,4,4,14,3,6,2,4,2,5,0,5,2,6,1,10,3,1,
      0,2,2,6,4,17,4,7,2,7,3,14,2,12,1,8,2,0,2,7,1,2,3,12,4,8,6,
      15,4,4,1,5,1,11,0,1,3,3,0,6,2,13,1,5,2,2,6,16,1,3,2,10,2,7,
      2,2,6,9,0,2,3,12,2,9,4,8,2,7,0,4,7,27,0,3,2,19,2,3,1,4,1,0,
      0,0,0,6,2,7,1,9,1,6,2,8,1,3,4,12,1,4,3,12,1,4,6,19,3,9,1,0,
      1,2,10,24,4,15,1,8,0,0,1,4,1,8,1,19,2,10,2,7,0,3,0,3,3,6,1,
      5,4,22,2,7,1,3,2,7,0,0,0,8,2,8,1,5,1,1,3,8,0,1,0,5,1,3,1,4,
      3,6,3,9,2,2,1,4,2,10,0,4,0,3,5,3,0,0,0,1,0,2,2,0,8,17,5,14,
      2,9,2,5,2,8,4,12,0,2,1,4,6,12,7,16,7,15,0,6};

   test_assert(mle_nm(x6, 2, 200, par));
   // Real value of alpha is 2.1, and p is (.18, .18, .64).

   test_assert(fabs(par->alpha-1.937450) < 1e-3);
   test_assert(fabs(par->p[0]-0.181793) < 1e-3);
   test_assert(fabs(par->p[1]-0.189070) < 1e-3);
   test_assert(fabs(par->p[2]-0.629137) < 1e-3);
   
   free(par);
   return;

}


void
test_mle_zinm
(void)
{

   // Cases checked by simulated annealing.
   
   zinm_par_t *par = new_zinm_par(1);
   if (par == NULL) {
      fprintf(stderr, "test error line %d\n", __LINE__);
      return;
   }
   
   // 0:53, 1:8, 2:9, 3:8, 4:4, 5:7, 6:3, 7:3, 8:1, 9:3, 10:1
   uint32_t x1[100] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,3,7,4,5,3,5,1,5,7,3,6,1,2,9,1,10,6,2,2,2,3,2,1,1,5,0,2,3,
      9,4,2,9,3,5,7,3,5,2,1,0,6,1,4,2,3,4,5,8,1 };

   test_assert_critical(mle_zinm(x1, 1, 100, par) != NULL);
   test_assert(fabs(par->alpha-3.6855) < 1e-3);
   test_assert(fabs(par->p[0]-0.5044) < 1e-3);
   test_assert(fabs(par->pi-0.5110) < 1e-3);

   // 0:73, 1:7, 2:7, 3:3, 4:4, 5:2, 7:1, 8:2, 9:1
   uint32_t x2[100] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,
      1,0,0,8,1,1,3,0,0,8,2,4,1,2,0,3,2,0,0,4,0,0,3,1,0,2,0,0,5,
      7,0,0,2,4,0,2,1,0,0,0,0,0,0,0,1,2,9,0,4 };

   test_assert_critical(mle_zinm(x2, 1, 100, par) != NULL);
   test_assert(fabs(par->alpha-1.8251) < 1e-3);
   test_assert(fabs(par->p[0]-0.4109) < 1e-3);
   test_assert(fabs(par->pi-0.3363) < 1e-3);

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
      {"tabulate", test_tabulate},
      {"compute_means", test_compute_means},
      {"eval_nb_f", test_eval_nb_f},
      {"eval_nb_dfda", test_eval_nb_dfda},
      {"eval_zinm_f", test_eval_zinm_f},
      {"eval_zinm_g", test_eval_zinm_g},
      {"eval_zinm_dfda", test_eval_zinm_dfda},
      {"eval_zinm_dfdp", test_eval_zinm_dfdp},
      {"eval_zinm_dgda", test_eval_zinm_dgda},
      {"ll_zinm", test_ll_zinm},
      {"nb_est_alpha", test_nb_est_alpha},
      {"mle_nm", test_mle_nm},
      {"mle_zinm", test_mle_zinm},
      {NULL, NULL}
   };

   return run_unittest(argc, argv, test_cases);

}
