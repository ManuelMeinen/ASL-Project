#ifndef __RADIX25_H__

#include <stdio.h>
#include <string.h>
#include <assert.h>
#define __RADIX25_H__

typedef struct {
    signed int d[10];
} radix25_t;

void print_radix25(radix25_t *x, char *string);
void radix25_toString( char* dest, radix25_t x );

#endif /* __RADIX25_H__ */