#include "../../_common_libs/radix25_curveParams.h"
#include "../../_common_libs/radix25.h"
#include "curve25519_base_256_25.h"
#include "../../_common_libs/common.h"
#include <string.h>
#include <stdlib.h>

int main(){
    init_curve25519_params_rad25(256);
    char** sk_A = (char**)malloc(sizeof(char*));
    char** sk_B = (char**)malloc(sizeof(char*));
    sk_A[0] = (char*)malloc(sizeof(char)*257);
    sk_B[0] = (char*)malloc(sizeof(char)*257);
    char keyA[257] = "0110100000111100010111001110100101110000011001011011010011111000100110000111011100010001001011110110110001100100010110010000101000101101001011000000110010101011000110100111101010001110000010100001010000011100100110001011000110010100110100101100111100010110";
    char keyB[257] = "0000000001100010110001011101101101111001001010000100011010111010110001100000011100001111100111011001000100010111000000110000100001010011010101000011011111100010000110011001011011111100011110010000010011010010001001101000111111110011110001100011110100111011";
    set_sk_length(strlen(keyA));
    radix25_t Px, Py;
    public_key_gen_base_256_25(keyA, &Px, &Py);

    strcpy(sk_A[0], keyA);
    strcpy(sk_B[0], keyB);
    if(validate_base_256_25(sk_A,sk_B,1)){
        printf("OK\n");
    }else{
        printf("Not OK\n");
    }
    //printf("%llX %llX %llX %llX %llX\n", Px.d[4], Px.d[3], Px.d[2], Px.d[1], Px.d[0]);
    //printf("%llX %llX %llX %llX %llX\n", Py.d[4], Py.d[3], Py.d[2], Py.d[1], Py.d[0]);
    return 0;
}