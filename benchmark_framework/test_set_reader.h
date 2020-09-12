#ifndef __TEST_SET_READER_H__
#include <stdio.h>
#include <string.h>
#define __TEST_SET_READER_H__

char** secrets_A;
char** secrets_B;
char** public_keys_dec;
char** public_keys_r17;
char** public_keys_r51;
char** public_keys_r25;

int read_test_set( char *rel_path_to_benchmarking_dir, int test_set_no, int do_read_pks, int *num_test_cases, int *secret_key_length );

#endif /* __TEST_SET_READER_ */