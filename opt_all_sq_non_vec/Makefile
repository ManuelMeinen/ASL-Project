all:
	gcc -mavx2 \
	-march=native \
	main_opt_all_sq_non_vec.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_17.c  \
	../_common_libs/vec_utils.c \
	../_common_libs/radix17_vec.c \
	../_common_libs/radix17_curveParams_vec.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

no-opt:
	gcc -mavx2 \
	-march=native \
	main_opt_all_sq_non_vec.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_17.c  \
	../_common_libs/vec_utils.c \
	../_common_libs/radix17_vec.c \
	../_common_libs/radix17_curveParams_vec.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

o3:
	gcc -O3 -mavx2 \
	-march=native \
	main_opt_all_sq_non_vec.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_17.c  \
	../_common_libs/vec_utils.c \
	../_common_libs/radix17_vec.c \
	../_common_libs/radix17_curveParams_vec.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

o2:
	gcc -O2 -mavx2 \
	-march=native \
	main_opt_all_sq_non_vec.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_17.c  \
	../_common_libs/vec_utils.c \
	../_common_libs/radix17_vec.c \
	../_common_libs/radix17_curveParams_vec.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

o3-no-vect:
	gcc -O3 -fno-tree-vectorize -mavx2 \
	-march=native \
	main_opt_all_sq_non_vec.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_17.c  \
	../_common_libs/vec_utils.c \
	../_common_libs/radix17_vec.c \
	../_common_libs/radix17_curveParams_vec.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

o3-avx2:
	gcc -O3 -mavx2 -mfma \
	-march=native \
	main_opt_all_sq_non_vec.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_17.c  \
	../_common_libs/vec_utils.c \
	../_common_libs/radix17_vec.c \
	../_common_libs/radix17_curveParams_vec.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

opt:
	gcc -O3 -mavx2 \
	-march=native \
	main_opt_all_sq_non_vec.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_17.c  \
	../_common_libs/vec_utils.c \
	../_common_libs/radix17_vec.c \
	../_common_libs/radix17_curveParams_vec.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out\
	 -fno-delete-null-pointer-checks -fno-schedule-insns2 -fno-optimize-sibling-calls -fno-loop-interchange -fno-rerun-cse-after-loop -fno-cse-skip-blocks -fno-unswitch-loops -fno-align-jumps -fno-expensive-optimizations -fno-tree-partial-pre -fno-ipa-cp -fno-tree-pre -fno-sched-interblock -fno-schedule-insns -fno-tree-loop-distribution -fno-ipa-sra -fno-strict-aliasing -fno-cse-follow-jumps -fno-caller-saves -fno-isolate-erroneous-paths-dereference -fno-optimize-strlen -fno-vect-cost-model -fno-partial-inlining -fno-devirtualize-speculatively -fno-align-labels -fno-tree-switch-conversion -fno-reorder-functions -fno-sched-spec -fno-loop-unroll-and-jam -fno-gcse-lm -fno-tree-loop-vectorize -fno-indirect-inlining -fno-thread-jumps -fno-crossjumping -fno-tree-tail-merge -fno-hoist-adjacent-loads -fno-tree-builtin-call-dce -fno-peel-loops -fno-devirtualize -fno-tree-slp-vectorize -fno-ipa-ra -fno-ipa-cp-clone -fno-ipa-icf -fno-reorder-blocks-and-partition -fno-predictive-commoning -fno-lra-remat -fno-tree-vrp -fno-align-loops -fno-gcse  

clean:
	rm out