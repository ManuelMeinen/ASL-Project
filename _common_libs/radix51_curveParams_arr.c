#include "radix51_curveParams_arr.h"


void init_curve25519_params_rad51_arr() {
    /*
    init A_rad51_arr, B_rad51_arr, Gx_rad51_arr, Gy_rad51_arr
    */
    A_rad51_arr [0] = 486662;
    B_rad51_arr [0] = 1;
    Gx_rad51_arr[0] = 9;
    for (int i = 1; i < 5; i++){
    	A_rad51_arr [i] = 0;
    	B_rad51_arr [i] = 0;
    	Gx_rad51_arr[i] = 0;
    }

    Gy_rad51_arr[0] = (signed long long) 498776679568345;
    Gy_rad51_arr[1] = (signed long long) 2196478692312381;
    Gy_rad51_arr[2] = (signed long long) 1446154986911029;
    Gy_rad51_arr[3] = (signed long long) 943862271874619;
    Gy_rad51_arr[4] = (signed long long) 574914022836744;    
}