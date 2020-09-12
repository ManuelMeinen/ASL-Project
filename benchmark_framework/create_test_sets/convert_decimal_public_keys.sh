#!/bin/bash

# ------------------------------------------------
# Convert the decimal public keys found in
# <test_set_num> test set to the desired format
# ------------------------------------------------

# --- get and check params
test_set_num="$1"
format_name="$2"

if [ -z "$test_set_num" ];   then echo "Missing parameter 'test_set_num'";   exit 1; fi
if [ -z "$format_name" ];    then echo "Missing parameter 'format_name'";    exit 1; fi

# --- change to directory of script
script_dir_name=$(dirname "$BASH_SOURCE")
if [ $script_dir_name != "." ]; then cd $script_dir_name; fi

# --- check if file exists containing the decimal public keys
bench_in_dir_name="../benchmark_in_files"
test_set_name="$bench_in_dir_name/test_set_${test_set_num}"
pubKey_fName="$test_set_name/public_keys.txt"

if [ ! -f "$pubKey_fName" ]; then
  echo " -- ERROR: No public key file exists for test set '$test_set_name'"
  exit 1;
fi

# --- determine conversion script name
case "$format_name" in
  "radix51" ) script_name="../../output_conversions/decimal_limbs51.py";;
  "radix17" ) script_name="../../output_conversions/decimal_limbs17.py";;
  "radix25" ) script_name="../../output_conversions/decimal_limbs25.py";;
  *) echo " -- ERROR: Unsupported format name '$format_name' found"; exit 1;;
esac

# --- handle existing output file
out_fName="$test_set_name/public_keys_${format_name}.txt"

if [ -f "$out_fName" ]; then rm "$out_fName"; fi

# --- convert outputs
while IFS= read -r line
do
  echo "$(python $script_name $line)" >> "$out_fName"
done < "$pubKey_fName"
