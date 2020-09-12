#ifndef _curve25519_base_256_17_h_
#include "../../_common_libs/radix17_vec.h"
#define _curve25519_base_256_17_h_

/* general functions (to be called by benchmarking) */
void scalar_mult_base_256_17(__m256i base_point_xlo, __m256i base_point_xhi, __m256i base_point_ylo, __m256i base_point_yhi, char *sk, __m256i *Pkxlo, __m256i *Pkxhi, __m256i *Pkylo, __m256i *Pkyhi);
void public_key_gen_base_256_17(char *sk, __m256i *xlo, __m256i *xhi, __m256i *ylo, __m256i *yhi);
int validate_base_256_17(char **sk_A, char **sk_B, int num_test_cases);
int checkPk_base_256_17(char **sk, char **pk, int num_test_cases);

#endif // _curve25519_base_256_17_h_