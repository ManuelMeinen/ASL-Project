#include "../../_common_libs/radix25_vec.h"
#include "../../_common_libs/radix25.h"
#include "../../_common_libs/vec_utils.h"
#include "bignum.h"

void test(){

/*  VALUES from BASE IMPLEMENTATION

    radix25_t a = {56979096, 8570560, 26308014, 24122168, 49219689, 1458285, 8700581, 25970563, 8273080, 7714699};//13311224122668955738146552185873637781933051739291637575097427091064598392472
    radix25_t b = {5336258, 25401515, 5634911, 29186739, 25621106, 14895217, 17733372, 14094260, 62168764, 6485312};//11189996199370894582589062329635885928504840703623090769010532534865250839746
    
    radix25_t res_a = {62315354, 417643, 31942926, 19754475, 7731932, 16353503, 26433953, 6510391, 3332981, 14200012}; //24501220322039850320735614515509523710437892442914728344107959625929849232218
    radix25_t res_s = {51642838, 16723477, 20673102, 28489861, 23598582, 20117500, 58076072, 11876302, 13213180, 1229386}; //2121227923298061155557489856237751853428211035668546806086894556199347552726
    radix25_t res_m = {46808148, 26304466, 38576209, 32196505, 41928472, 15787783, 31495827, 31602133, 63518558, 30972522};//53441124274223813656485740086574774531223107401464235683382378494160525737044
    radix25_t res_d = {58711383, 9900831, 1001900, 11893292, 65995188, 21577140, 49008858, 22951153, 52761445, 13575384};//23423465588619757183253169998465881899659595894318503666586588616682836581719
    radix25_t res_i = {43066194, 11821405, 32495562, 26162991, 9136424, 21873032, 25768362, 1559236, 9537546, 2435527};//4202347638909874817463795626263865895879347174602030051662579480283193615186
    radix25_t res_sq = {50930837, 20740871, 14005878, 20138371, 10482029, 16231300, 44082054, 33522090, 51821034, 13706240};//23649249294429756376717901120500673641464092031992322068868721792945861764245

*/
    //radix25_t a = {56979096, 8570560, 26308014, 24122168, 49219689, 1458285, 0, 0, 0, 0};
    // __m256i a_lo = _mm256_setr_epi32(56979096, 8570560, 26308014, 24122168, 49219689, 1458285, 0, 0);
    // __m256i a_hi = _mm256_setr_epi32(0, 0, 0, 0, 0, 0, 0, 0);

	__m256i a_lo = _mm256_setr_epi32(56979096, 8570560, 26308014, 24122168, 49219689, 1458285, 8700581, 25970563);
    __m256i a_hi = _mm256_setr_epi32(8273080, 7714699, 0, 0, 0, 0, 0, 0);
    //debut sqr
    //50106187, 21029885, 57791077, 12795838, 16475153, 29972673, 59712616, 11747338, 22499141, 11357625,
    // __m256i a_lo = _mm256_setr_epi32(50106187, 21029885, 57791077, 12795838, 16475153, 29972673, 59712616, 11747338);
    // __m256i a_hi = _mm256_setr_epi32(22499141, 11357625, 0, 0, 0, 0, 0, 0);


    // __m256i a_lo = _mm256_setr_epi32(1, 2, 3, 4, 5, 6, 0, 1);
    // __m256i a_hi = _mm256_setr_epi32(1, 1, 0, 0, 0, 0, 0, 0);

    //radix25_t b = {5336258, 25401515, 5634911, 29186739, 25621106, 14895217, 0, 0, 0, 0};
    // __m256i b_lo = _mm256_setr_epi32(5336258, 25401515, 5634911, 29186739, 25621106, 14895217, 0, 0);
    // __m256i b_hi = _mm256_setr_epi32(0, 0, 0, 0, 0, 0, 0, 0);

	__m256i b_lo = _mm256_setr_epi32(5336258, 25401515, 5634911, 29186739, 25621106, 14895217, 17733372, 14094260);
    __m256i b_hi = _mm256_setr_epi32(62168764, 6485312, 0, 0, 0, 0, 0, 0);

    // __m256i b_lo = _mm256_setr_epi32(1, 2, 3, 4, 5, 6, 0, 1);
    // __m256i b_hi = _mm256_setr_epi32(1, 1, 0, 0, 0, 0, 0, 0);

    //radix25_t res_a = {62315354, 417643, 31942926, 19754475, 7731932, 16353503, 0, 0, 0, 0};
    __m256i res_a_lo = _mm256_setr_epi32(62315354, 417643, 31942926, 19754475, 7731932, 16353503, 26433953, 6510391);
    __m256i res_a_hi = _mm256_setr_epi32(3332981, 14200012, 0, 0, 0, 0, 0, 0);

    //radix25_t res_s = {51642819, 16723477, 20673102, 28489861, 23598582, 20117500, 67108863, 33554431, 67108863, 33554431};
    //__m256i res_s_lo = _mm256_setr_epi32(51642819, 16723477, 20673102, 28489861, 23598582, 20117500, 67108863, 33554431);
    //__m256i res_s_hi = _mm256_setr_epi32(67108863, 33554431, 0, 0, 0, 0, 0, 0);
	
	__m256i res_s_lo = _mm256_setr_epi32(51642838, 16723477, 20673102, 28489861, 23598582, 20117500, 58076072, 11876302);
    __m256i res_s_hi = _mm256_setr_epi32(13213180, 1229386, 0, 0, 0, 0, 0, 0);

    //radix25_t res_m = {3328876, 15174893, 49551188, 488918, 32091215, 18874244, 9011022, 22400410, 19166965, 28210162};
    __m256i res_m_lo = _mm256_setr_epi32(46808148, 26304466, 38576209, 32196505, 41928472, 15787783, 31495827, 31602133);
    __m256i res_m_hi = _mm256_setr_epi32(63518558, 30972522, 0, 0, 0, 0, 0, 0);
    //  __m256i res_m_lo = _mm256_setr_epi32(1939, 460, 964, 438, 545, 94, 205, 116);
    // __m256i res_m_hi = _mm256_setr_epi32(169, 72, 0, 0, 0, 0, 0, 0);

    __m256i res_d_lo = _mm256_setr_epi32(58711383, 9900831, 1001900, 11893292, 65995188, 21577140, 49008858, 22951153);
    __m256i res_d_hi = _mm256_setr_epi32(52761445, 13575384, 0, 0, 0, 0, 0, 0);

    __m256i res_i_lo = _mm256_setr_epi32(43066194, 11821405, 32495562, 26162991, 9136424, 21873032, 25768362, 1559236);
    __m256i res_i_hi = _mm256_setr_epi32(9537546, 2435527, 0, 0, 0, 0, 0, 0);

    //__m256i res_sq_lo = _mm256_setr_epi32(50930837, 20740871, 14005878, 20138371, 10482029, 16231300, 44082054, 33522090);
    //__m256i res_sq_hi = _mm256_setr_epi32(51821034, 13706240, 0, 0, 0, 0, 0, 0);
    //debug sqr
    //39907819, 12753165, 44623054, 2652141, 28360218, 367167, 17971751, 18749786, 25524492, 26432548, 
    __m256i res_sq_lo = _mm256_setr_epi32( 27013257, 8056918, 29673231, 20923677, 40959805, 5476750, 62486760, 25945319);
    __m256i res_sq_hi = _mm256_setr_epi32( 38104616, 32223607, 0, 0, 0, 0, 0, 0);

    //radix25_t ra,rs,rm,rd,ri,rsq;
    __m256i ra_lo,rs_lo,rm_lo,rd_lo,ri_lo,rsq_lo,
            ra_hi,rs_hi,rm_hi,rd_hi,ri_hi,rsq_hi;
    set_ui_256(&ra_lo, &ra_hi, 0);
    set_ui_256(&rs_lo, &rs_hi, 0);
    set_ui_256(&rm_lo, &rm_hi, 0);
    set_ui_256(&rd_lo, &rd_hi, 0);
    set_ui_256(&ri_lo, &ri_hi, 0);
    set_ui_256(&rsq_lo, &rsq_hi, 0);

    //add_256(&ra_lo, &ra_hi, a_lo, a_hi, b_lo, b_hi); reduce_2_to_2( &ra_lo, &ra_hi );
    //sub_256(&rs_lo, &rs_hi, a_lo, a_hi, b_lo, b_hi);
    //mul_256(&rm_lo, &rm_hi, a_lo, a_hi, b_lo, b_hi);
    //div_256(&rd_lo, &rd_hi, a_lo, a_hi, b_lo, b_hi);
    inverse_256(&ri_lo, &ri_hi, a_lo, a_hi);
    //square_256(&rsq_lo, &rsq_hi, a_lo, a_hi);

    //if(cmp_coord_256(ra_lo, ra_hi, res_a_lo, res_a_hi)){
    //    printf("Add passed\n");
    //    //print_vec_int(ra_lo   , "  ra_lo   ");
    //    //print_vec_int(ra_hi   , "  ra_hi   ");
    //}else{
    //	printf( "Add failed: \n" );
    //    print_vec_int(ra_lo   , "  ra_lo   ");
    //    print_vec_int(res_a_lo, "  expected");
    //    print_vec_int(ra_hi   , "  ra_hi   ");
    //    print_vec_int(res_a_hi, "  expected");
    //}
    //
    //if(cmp_coord_256(rs_lo, rs_hi, res_s_lo, res_s_hi)){
    //    printf("Sub passed\n");
    //}else{
    //	printf( "Sub failed: \n" );
    //    print_vec_int(rs_lo   , "  rs_lo   ");
    //    print_vec_int(res_s_lo, "  expected");
    //    print_vec_int(rs_hi   , "  rs_hi   ");
    //    print_vec_int(res_s_hi, "  expected");
    //}
    //
    //if(cmp_coord_256(rm_lo, rm_hi, res_m_lo, res_m_hi)){
    //    printf("Mul passed\n");
    //}else{
    //	printf( "Mul failed: \n" );
    //    print_vec_int(rm_lo   , "  rm_lo   ");
    //    print_vec_int(res_m_lo, "  expected");
    //    print_vec_int(rm_hi   , "  rm_hi   ");
    //    print_vec_int(res_m_hi, "  expected");
    //}
//
    // if(cmp_coord_256(rsq_lo, rsq_hi, res_sq_lo, res_sq_hi)){
    //     printf("Square passed\n");
    // }else{
    //  	printf( "Square failed: \n" );
    //     print_vec_int(rsq_lo   , "  rsq_lo  ");
    //     print_vec_int(res_sq_lo, "  expected");
    //     print_vec_int(rsq_hi   , "  rsq_hi  ");
    //     print_vec_int(res_sq_hi, "  expected");
    // }
//
    //if(cmp_coord_256(rd_lo, rd_hi, res_d_lo, res_d_hi)){
    //    printf("Div passed\n");
    //}else{
    // 	printf( "Div failed: \n" );
    //    print_vec_int(rd_lo   , "  rd_lo   ");
    //    print_vec_int(res_d_lo, "  expected");
    //    print_vec_int(rd_hi   , "  rd_hi   ");
    //    print_vec_int(res_d_hi, "  expected");
    //}

    if(cmp_coord_256(ri_lo, ri_hi, res_i_lo, res_i_hi)){
      printf("Inv passed\n");
    }else{
       printf( "Inv failed: \n" );
       print_vec_int(ri_lo   , "  ri_lo   ");
       print_vec_int(res_i_lo, "  expected");
       print_vec_int(ri_hi   , "  ri_hi   ");
       print_vec_int(res_i_hi, "  expected");
    }

}