#ifndef _curve25519_base_256_51_h_
#define _curve25519_base_256_51_h_

/* general functions (to be called by benchmarking) */
void public_key_gen_opt_all_51(char *sk, long Pkx[5], long Pky[5]);
void scalar_mult_opt_all_51(long base_point_x[5], long base_point_y[5], char *sk, long Pkx[5], long Pky[5]);
int validate_opt_all_51(char **sk_A, char **sk_B, int num_test_cases);
int checkPk_opt_all_51(char **sk, char **pk, int num_test_cases);

#endif // _curve25519_base_256_51_h_