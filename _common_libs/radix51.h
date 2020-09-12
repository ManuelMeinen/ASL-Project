#ifndef __RADIX51_H__
#include <stdio.h>
#include <string.h>
#define __RADIX51_H__

typedef struct {
    signed long long d[5];
} radix51_t;

void print_radix51(radix51_t *x, char *string);
void radix51_toString( char* dest, radix51_t x );

#endif /* __RADIX51_H__ */