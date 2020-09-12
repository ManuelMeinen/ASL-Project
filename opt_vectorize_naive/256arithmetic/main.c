#include "bignum.h"
#include "../../_common_libs/radix17.h"
#include <gmp.h>
#include "test.c"
void print(radix17_t r){

    for(int i =0;i < 15; i++){
        printf("%d, ",r.d[i]);
    }
    printf("\n");
}

int main(int argc, char const *argv[])
{
    // init();
    // char sk[257] = "0110100000111100010111001110100101110000011001011011010011111000100110000111011100010001001011110110110001100100010110010000101000101101001011000000110010101011000110100111101010001110000010100001010000011100100110001011000110010100110100101100111100010110"; //change this

    // //public_key_generation(sk);
    
    // //validate();

     mpz_t x, y, inv, p, res, res1;
    // mpz_init(x);
    mpz_init(y);
     mpz_init(res);
     mpz_init(res1);
        mpz_init(p);
    mpz_init(x);
    mpz_init(inv);
    mpz_set_str(x, "1D6DE2C7E3CB8C623C1909854A16406DBBC21A6E02670C8B6D720B1B03656E98", 16);
    mpz_set_str(p, "7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffed", 16);

    mpz_set_str(y, "18BD503B49EBC6B87DA21D2DF8E3487161BC9CB7AB5662AFDAFE0E62AC516CC2", 16);
    

    mpz_invert(inv, x, p);
    // mpz_mul(res, x, inv);
    // mpz_mod(res1, res, p);
    mpz_clear(inv);

    //print_mpz("gmp result: ", inv);
    gmp_printf("%s: %Zd\n", "gmp: ", inv);

    radix17_t a = {93848, 98738, 33478, 93614, 28872, 94227, 67689, 56183, 5696, 49829, 67138, 101447, 15544, 71231, 30135};
    radix17_t a_2 = {93848, 98738, 33478, 93614, 28872, 94227, 67689, 56183, 5696, 49829, 67138, 101447, 15544, 71231, 30135};
    radix17_t b = {93378, 87592, 99224, 129887, 91690, 114010, 62066, 58051, 58184, 38652, 92295, 55055, 40636, 33242, 25333};
    radix17_t r = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    test();

    return 0;
}

