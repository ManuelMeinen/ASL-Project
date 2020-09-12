#include "../../_common_libs/vec_utils.h"
#include <immintrin.h>
#include <emmintrin.h>
#include "bignum.h"

void test(){

    //radix17_t a = {93848, 98738, 33478, 93614, 28872, 94227, 67689, 56183, 5696, 49829, 67138, 101447, 15544, 71231, 30135};

    __m256i a_lo = _mm256_setr_epi32(93848, 98738, 33478, 93614, 28872, 94227, 67689, 56183);
    __m256i a_hi = _mm256_setr_epi32(5696, 49829, 67138, 101447, 15544, 71231, 30135, 0);         // most upper zero 
    
    //radix17_t b = {93378, 87592, 99224, 129887, 91690, 114010, 62066, 58051, 58184, 38652, 92295, 55055, 40636, 33242, 25333};
    __m256i b_lo = _mm256_setr_epi32(93378, 87592, 99224, 129887, 91690, 114010, 62066, 58051);
    __m256i b_hi = _mm256_setr_epi32(58184, 38652, 92295, 55055, 40636, 33242, 25333, 0);

    //radix17_t res_a = {56154, 55259, 1631, 92430, 120563, 77165, 129756, 114234, 63880, 88481, 28361, 25431, 56181, 104473, 55468};
    __m256i res_a_lo = _mm256_setr_epi32(56154, 55259, 1631, 92430, 120563, 77165, 129756, 114234);
    __m256i res_a_hi = _mm256_setr_epi32(63880, 88481, 28361, 25431, 56181, 104473, 55468, 0);  

    //radix17_t res_s = {470, 11146, 65326, 94798, 68253, 111288, 5622, 129204, 78583, 11176, 105915, 46391, 105980, 37988, 4802};
    __m256i res_s_lo = _mm256_setr_epi32(470, 11146, 65326, 94798, 68253, 111288, 5622, 129204);
    __m256i res_s_hi = _mm256_setr_epi32(78583, 11176, 105915, 46391, 105980, 37988, 4802, 0); 

    //radix17_t res_m = {15444, 107877, 102751, 41041, 78630, 125767, 116504, 3903, 61671, 38547, 109296, 123445, 79710, 54756, 120986};
    __m256i res_m_lo = _mm256_setr_epi32(15444, 107877, 102751, 41041, 78630, 125767, 116504, 3903);
    __m256i res_m_hi = _mm256_setr_epi32(61671, 38547, 109296, 123445, 79710, 54756, 120986, 0); 

    //radix17_t res_d = {122199, 16319, 38675, 84396, 22535, 46458, 65972, 92663, 84285, 119002, 123765, 89652, 70501, 110994, 53028};
    __m256i res_d_lo = _mm256_setr_epi32(122199, 16319, 38675, 84396, 22535, 46458, 65972, 92663);
    __m256i res_d_hi = _mm256_setr_epi32(84285, 119002, 123765, 89652, 70501, 110994, 53028, 0); 

    //radix17_t res_i = {74578, 47944, 46177, 120778, 24311, 102199, 92456, 69701, 85441, 78250, 100548, 6090, 100362, 101960, 9513};
    __m256i res_i_lo = _mm256_setr_epi32(74578, 47944, 46177, 120778, 24311, 102199, 92456, 69701);
    __m256i res_i_hi = _mm256_setr_epi32(85441, 78250, 100548, 6090, 100362, 101960, 9513, 0); 

    //radix17_t res_sq = {74901, 3972, 81019, 112246, 67178, 78665, 127341, 67663, 63403, 41862, 87376, 130945, 47594, 395, 53540};
    __m256i res_sq_lo = _mm256_setr_epi32(74901, 3972, 81019, 112246, 67178, 78665, 127341, 67663);
    __m256i res_sq_hi = _mm256_setr_epi32(63403, 41862, 87376, 130945, 47594, 395, 53540, 0);
    
    __m256i ra_lo,ra_hi, rs_lo, rs_hi, rm_lo,rm_hi,rd_lo, rd_hi, ri_lo, ri_hi, rsq_lo, rsq_hi;
    set_ui_256(&ra_lo, &ra_hi, 0);
    set_ui_256(&rs_lo, &rs_hi, 0);
    set_ui_256(&rm_lo, &rm_hi, 0);
    set_ui_256(&rd_lo, &rd_hi, 0);
    set_ui_256(&ri_lo, &ri_hi, 0);
    set_ui_256(&rsq_lo, &rsq_hi, 0);

    add_256(&ra_lo, &ra_hi, a_lo, a_hi, b_lo, b_hi);
    sub_256(&rs_lo, &rs_hi, a_lo, a_hi, b_lo, b_hi);
    mul_256(&rm_lo, &rm_hi, a_lo, a_hi, b_lo, b_hi);
    div_256(&rd_lo, &rd_hi, a_lo, a_hi, b_lo, b_hi);
    inverse_256(&ri_lo, &ri_hi, a_lo, a_hi);
    square_256(&rsq_lo, &rsq_hi, a_lo, a_hi);

    if(cmp_coord_256(ra_lo, ra_hi, res_a_lo, res_a_hi)){
        printf("Add passed\n");
    }else{
        print_vec_int(ra_lo, "ra_lo");
        print_vec_int(ra_hi, "ra_hi");
        print_vec_int(res_a_lo, "res_a_lo");
        print_vec_int(res_a_hi, "res_a_hi");
        printf("Add failed\n");
    }
    if(cmp_coord_256(rs_lo, rs_hi, res_s_lo, res_s_hi)){
        printf("Sub passed\n");
    }
    if(cmp_coord_256(rm_lo, rm_hi, res_m_lo, res_m_hi)){
        printf("Mul passed\n");
    }else{
        print_vec_int(rm_lo, "rm_lo");
        print_vec_int(rm_hi, "rm_hi");
        print_vec_int(res_m_lo, "res_m_lo");
        print_vec_int(res_m_hi, "res_m_hi");
        printf("Mul failed\n");
    }

    if(cmp_coord_256(rsq_lo, rsq_hi, res_sq_lo, res_sq_hi)){
        printf("Square passed\n");
    }else{
        print_vec_int(rsq_lo, "rsq_lo");
        print_vec_int(rsq_hi, "rsq_hi");
        print_vec_int(res_sq_lo, "res_sq_lo");
        print_vec_int(res_sq_hi, "res_sq_hi");
        printf("Square failed\n");
    }

    if(cmp_coord_256(rd_lo, rd_hi, res_d_lo, res_d_hi)){
        printf("Div passed\n");
    }

     if(cmp_coord_256(ri_lo, ri_hi, res_i_lo, res_i_hi)){
        printf("Inv passed\n");
    }
}
