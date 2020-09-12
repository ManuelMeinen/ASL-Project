#!/bin/bash

# ------------------------------------------------
# Create a test set for the benchmarking framework
# ------------------------------------------------

# --- get and check params
test_set_num="$1"
num_test_cases="$2"
secret_length="$3"

if [ -z "$test_set_num" ];   then echo "Missing parameter 'test_set_num'";   exit 1; fi
if [ -z "$num_test_cases" ]; then echo "Missing parameter 'num_test_cases'"; exit 1; fi
if [ -z "$secret_length" ];  then echo "Missing parameter 'secret_length'";  exit 1; fi

# --- change to directory of script
script_dir_name=$(dirname "$BASH_SOURCE")
if [ $script_dir_name != "." ]; then cd $script_dir_name; fi

# --- set up test set directory
bench_in_dir_name="../benchmark_in_files"
test_set_name="$bench_in_dir_name/test_set_${test_set_num}"

if [ -d "$test_set_name" ]; then
  while true; do
    read -p "  Warning: test set $test_set_num already exists. Do you wish to continue anyway (previous test set will be deleted)? (yes/no) " answer
    case $answer in
        [Yy]* ) break;;
        [Nn]* ) exit;;
        * ) echo "Please answer 'yes' or 'no'.";;
    esac
  done
  rm -r "$test_set_name"
fi

mkdir -p "$test_set_name"

# --- add general parameters
echo "$num_test_cases" >  "$test_set_name/general_parameters.txt"
echo "$secret_length"  >> "$test_set_name/general_parameters.txt"

# --- create secret key files (for Alice and Bob)
for (( i=1; i<=$num_test_cases; i++ )); do
  echo "$(python rand_256bit_sk.py)" >> "$test_set_name/secret_keys.txt"
  echo "$(python rand_256bit_sk.py)" >> "$test_set_name/secret_keys_B.txt"
done

# --- generate public keys for Alice using the base GMP implementation
if [ -e "create_ground_truths" ]; then make clear &> /dev/null; fi
make &> /dev/null
if [ "$?" != 0 ]; then
  echo " -- ERROR: failed to compile c program to create ground truths"
  rm -r "$test_set_name"
  exit 1;
fi
./create_ground_truths "$test_set_num" > "$test_set_name/public_keys.txt"
if [ "$?" != 0 ]; then
  echo " -- ERROR: failed to generate the public keys. $(tail $test_set_name/public_keys.txt)"
  rm -r "$test_set_name"
  #rm "$test_set_name/public_keys.txt"
  exit 1;
fi

# --- convert decimal public keys to various other formats
./convert_decimal_public_keys.sh "$test_set_num" "radix51"
./convert_decimal_public_keys.sh "$test_set_num" "radix17"
./convert_decimal_public_keys.sh "$test_set_num" "radix25"

