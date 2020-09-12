#include "bignum.h"
#include "../../_common_libs/radix25.h"
#include "test.c"
void print(radix25_t r){

    for(int i =0;i < 10; i++){
        printf("%llu, ",r.d[i]);
    }
    printf("\n");
}

int main(int argc, char const *argv[])
{
    test();

    return 0;
}

