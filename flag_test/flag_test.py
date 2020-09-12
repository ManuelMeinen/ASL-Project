import subprocess
import pandas as pd

path = '../opt_all/'
files_opt_all = ["main_opt_all.c", "../_common_libs/common.c", "256arithmetic/bignum.c", "curve25519/curve25519_base_256_17.c", "../_common_libs/radix17_vec.c", "../_common_libs/radix17_curveParams_vec.c", "../benchmark_framework/benchmarking.c", "../benchmark_framework/test_set_reader.c"]
mandatory_flags = ["-mavx2","-march=native"]
all_flags = []

def make(flags, files):
    '''
    Compile the files with the given flags
    IN: str[] flags str[] files
    '''
    cmd = "gcc -O3"
    for flag in mandatory_flags:
        cmd += " "+flag
    for file in files:
        cmd += " "+file
    for flag in flags:
        cmd += " "+flag
    cmd += " -o out"
    #print(cmd)
    try:
        out = subprocess.check_output(cmd, shell=True, cwd=path)
    except:
        print("something went wrong")
    

def is_valid_flag(flags, files):
    '''

    IN: str[] flags, str[] files
    OUT: True iff files compile with flags successfully (no errors and no warnings)
    '''
    cmd = "gcc"
    for flag in mandatory_flags:
        cmd += " "+flag
    for file in files:
        cmd += " "+file
    for flag in flags:
        cmd += " "+flag
    cmd += " -o out"
    try:
        out = subprocess.check_output(cmd, shell=True, cwd=path, stderr=subprocess.STDOUT).decode("utf-8")
        if len(out)==0:
            return True
        else:
            return False
    except:
        return False
    
def determine_valid_flags():
    '''
    Get all flags from the csv files, determine the valid ones and save them in valid_flags.csv
    '''
    o1, o2, o3 = get_flags()
    valid_o1_flags = []
    valid_o2_flags = []
    valid_o3_flags = []
    for flag in o3:
        if is_valid_flag([flag], files_opt_all):
            valid_o1_flags.append(flag)
    all_flags = valid_o1_flags
    for flag in o2:
        if is_valid_flag([flag], files_opt_all):
            valid_o2_flags.append(flag)
    all_flags += valid_o2_flags
    for flag in o3:
        if is_valid_flag([flag], files_opt_all):
            valid_o3_flags.append(flag)
    all_flags += valid_o3_flags
    print(all_flags)
    df = pd.DataFrame(all_flags)
    df.to_csv('valid_flags.csv', header=False, index=False)
    

def get_cycles(test_set):
    '''
    Run benchmarking on test set test_set and return the number of cycles needed
    IN: int test_set
    OUT: float cycles
    '''
    cycles = subprocess.check_output("./out "+str(test_set)+" b 0 res", shell=True, cwd=path)
    return float(cycles)

def get_flags():
    '''
    Read the flags for each optimization level from the csv files
    OUT: str[] o1, str[] o2, str[] o3
    '''
    o1 = pd.read_csv('flagsO1.csv')
    o1.columns = ['flag']
    o2 = pd.read_csv('flagsO2.csv')
    o2.columns = ['flag']
    o3 = pd.read_csv('flagsO3.csv')
    o3.columns = ['flag']
    return o1['flag'].values.tolist(), o2['flag'].values.tolist(), o3['flag'].values.tolist()

def get_valid_flags():
    '''
    Read the valid flags from a csv
    OUT: str[] flags
    '''
    flags = pd.read_csv('valid_flags.csv')
    flags.columns = ['flag']
    return flags['flag'].values.tolist()

def add_one_by_one_strategy():
    # run with no flags
    valid_flags = get_valid_flags()
    opt_flags = []
    make(opt_flags, files_opt_all)
    total_min_cycles = get_cycles(1)
    min_cycles = 0
    
    while min_cycles<total_min_cycles:
        min_cycles = total_min_cycles
        best_flag = ""
        for flag in valid_flags:
            if not flag in opt_flags:
                test_flags = opt_flags + [flag]
                make(test_flags, files_opt_all)
                cycles = get_cycles(1)
                if cycles<min_cycles:
                    best_flag =  flag
                    min_cycles = cycles
        opt_flags.append(best_flag)
    df = pd.DataFrame(opt_flags)
    df.to_csv('opt_flags.csv', index=False, header=False)
    print("The best result was produced by using the following flags:")
    print(df)
    print("The min no. of cycles is: "+str(total_min_cycles))
        
def run_opt():
    flags = pd.read_csv('valid_flags.csv')
    flags.columns = ['flag']
    flags = flags['flag'].values.tolist()
    make(flags, files_opt_all)
    cycles = get_cycles(1)
    print("The number of cycles for opt flags is: "+str(cycles))

def main():
    #determine_valid_flags()
    add_one_by_one_strategy()
    #run_opt()

if __name__ == '__main__':
    main()