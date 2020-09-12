#ifndef _curve25519_base_256_25_h_
#include "../../_common_libs/radix25_vec.h"
#define _curve25519_base_256_25_h_

/* general functions (to be called by benchmarking) */
void scalar_mult_opt_all_25(__m256i base_point_xlo, __m256i base_point_xhi, __m256i base_point_ylo, __m256i base_point_yhi, char *sk, __m256i *Pkxlo, __m256i *Pkxhi, __m256i *Pkylo, __m256i *Pkyhi);
void public_key_gen_opt_all_25(char *sk, __m256i *xlo, __m256i *xhi, __m256i *ylo, __m256i *yhi);
int validate_opt_all_25(char **sk_A, char **sk_B, int num_test_cases);
int checkPk_opt_all_25(char **sk, char **pk, int num_test_cases);

#endif // _curve25519_base_256_25_h_