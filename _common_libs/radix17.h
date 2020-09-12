#ifndef __RADIX17_H__
#include <stdio.h>
#include <string.h>
#include <assert.h>
#define __RADIX17_H__

typedef struct {
    signed int d[16];
} radix17_t;
 

void print_radix17(radix17_t *x, char *string);
void radix17_toString( char* dest, radix17_t x );
#endif /* __RADIX17_H__ */