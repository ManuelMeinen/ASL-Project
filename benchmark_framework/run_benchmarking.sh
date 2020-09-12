#!/bin/bash

# ------------------------------------------------
# Auxiliary functions
# ------------------------------------------------

function GetTestSetStr {
  typeset res=""
  for i in "${test_set_numbers[@]}"
  do
    res="$res $i"
  done
  echo $res
}

# ------------------------------------------------
# Run benchmarking framework
# ------------------------------------------------

# --- change to directory of script
script_dir_name=$(dirname "$BASH_SOURCE")
if [ $script_dir_name != "." ]; then cd $script_dir_name; fi
baseDir=$(pwd)

# --- specify benchmarking parameters
#     CHANGE THESE IF WANTING TO RUN FOR DIFFERENT COMPILER FLAGS OR TEST SETS
compiler_flag_ids=( "no-opt" "o2" "o3-no-vect" "o3")
test_set_numbers=( 1 )
func_ids=("base_256_17" "opt_scalar_replacement" "opt_scalar_repl_and_vectorize" "opt_curve" "opt_all" "base_256_25" "opt_all_25" "base_256_51" "opt_all_51")  # the values specified here must correspond to the directory names containing the implementation

# --- specify run mode of benchmarking
#     COMMENT OUT THE APPLICABLE VALUE!
#run_mode="v"  # --> VALIDATION-ONLY
run_mode="b"  # --> BENCHMARKING-ONLY
#run_mode="f"  # --> FULL BENCHMARKING (validation + benchmarking)

# --- specify verbosity
be_verbose="0" # false
#be_verbose="1" # true

# --- set up results directory
bench_out_dir_name="benchmark_results/$(date -u '+20%y%m%d.%H%M%S')"
mkdir -p "$bench_out_dir_name"

for cmp_flag in ${compiler_flag_ids[@]}; do
  mkdir -p "$bench_out_dir_name/$cmp_flag"
done

# --- run benchmarking
for fId in ${func_ids[@]}; do
	if [ "$be_verbose" == "1" ]; then
		echo "###"
		echo "### Running benchmarking for function '$fId' and test set(s) '$(GetTestSetStr)'"
		echo "###"
	fi
	# change to directory of function implementation
	if [ ! -d "../$fId" ]; then
		echo " -- ERROR: directory corresponding to implementation '$fId' does not exist"
		exit 1;
	fi
	cd "../$fId"
	# run benchmarking for each of the compiler flags specified
	for cmp_flag in ${compiler_flag_ids[@]}; do
		# compile with specified compiler flag
		# res=$(make $cmp_flag) #&> /dev/null
		#if [ "$?" != 0 ]; then echo "$res"; echo ""; echo " -- ERROR: failed to compile out file"; exit 1; fi
		make $cmp_flag &> /dev/null
		if [ "$?" != 0 ]; then echo ""; echo " -- ERROR: failed to compile out file"; exit 1; fi
		if [ "$be_verbose" == "1" ]; then
			echo "### Compiled with option '$cmp_flag'"; echo ""
		fi
		for num in ${test_set_numbers[@]}; do
			if [ "$be_verbose" == "1" ]; then
				./out $num $run_mode $be_verbose "../benchmark_framework/$bench_out_dir_name/$cmp_flag/"
			else
				perf=$(./out $num $run_mode $be_verbose "../benchmark_framework/$bench_out_dir_name/$cmp_flag/")
				echo "$perf,$num,$fId,$cmp_flag"
			fi
		done
		make clear &> /dev/null
	done
	if [ "$be_verbose" == "1" ]; then echo ""; fi
	cd "$baseDir"
done

exit 0;