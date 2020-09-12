#ifndef __BENCHMARKING_H__
#include <stdio.h>
#include <string.h>
#define __BENCHMARKING_H__

typedef int ( pubKeyGen_func )( char**, int );            // 1. char**: secret keys, 2. int: number of test cases
typedef int  ( validate_func  )( char**, char**, int );    // 1. char**: secret keys Alice, 2. char**: secret keys Bob, 3. int: number of test cases
typedef int  ( checkPk_func   )( char**, char**, int );    // 1. char**: secret keys Alice, 2. char**: public keys 3. int: number of test cases

int run_benchmarking( char *rel_path_to_benchmarking_dir, int test_set_no, char run_mode, int be_verbose, char *res_dir_name, pubKeyGen_func *f, validate_func *v, checkPk_func *c, char *pubKeyFormat, char *name, char *id );

#endif /* __BENCHMARKING_H__ */