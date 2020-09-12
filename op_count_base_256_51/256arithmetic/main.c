#include "bignum.h"
#include <gmp.h>

int main(int argc, char const *argv[])
{
    init();
    char sk[257] = "0110100000111100010111001110100101110000011001011011010011111000100110000111011100010001001011110110110001100100010110010000101000101101001011000000110010101011000110100111101010001110000010100001010000011100100110001011000110010100110100101100111100010110"; //change this

    public_key_generation(sk);

    validate();

    mpz_t x, y, inv, p, res, res1;
    mpz_init(y);
    mpz_init(res);
    mpz_init(res1);
    mpz_init(p);
    mpz_init(x);
    mpz_init(inv);
    mpz_set_str(x, "1D6DE2C7E3CB8C623C1909854A16406DBBC21A6E02670C8B6D720B1B03656E98", 16);
    mpz_set_str(p, "7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffed", 16);

    mpz_set_str(y, "18BD503B49EBC6B87DA21D2DF8E3487161BC9CB7AB5662AFDAFE0E62AC516CC2", 16);

    mpz_invert(inv, y, p);
    mpz_mul(res, x, inv);
    mpz_mod(res1, res, p);
    mpz_clear(inv);

    gmp_printf("%s: %Zd\n", "gmp: ", res1);

    radix51_t a = {575160602422936, 1618811318005166, 97863898957929, 1742854989071013, 517724694265016};
    radix51_t b = {1704666820865218, 1958688903789407, 999601117524594, 945849795254012, 435221983174332};
    radix51_t r = {0, 0, 0, 0, 0};

    div_256(&r, &a, &b);

    printf("256 result: %llu, %llu, %llu, %llu, %llu\n", r.d[0], r.d[1], r.d[2], r.d[3], r.d[4]);

    return 0;
}