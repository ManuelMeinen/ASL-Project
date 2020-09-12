#include "bignum.h"
#include "../../_common_libs/radix25_vec.h"
#include "../../_common_libs/vec_utils.h"
#include "test.c"


int main(int argc, char const *argv[])
{
    init_vec_constants();
    test();

    return 0;
}

