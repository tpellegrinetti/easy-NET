#!/bin/bash

# Initialize variables with default values
folder="."
cutoff=60
relative=false

# Function to calculate relative abundance
calculate_relative() {
    awk -F'\t' -v cutoff="$cutoff" '{
        if (NR == 1) {
            print;
            next;
        }
        sum = 0;
        for (i = 2; i <= NF; i++) sum += $i;
        if (NR == 2) total = sum;
        else total += sum;
        line[NR - 1] = $0;
        sum_line[NR - 1] = sum;
    } END {
        for (i in line) {
            if (((sum_line[i] / total)*100) > cutoff) print line[i];
        }
    }' $1
}

# Process command-line arguments
while getopts "n:a:d:" opt; do
    case $opt in
        n) cutoff=$OPTARG ;;
        a) cutoff=$OPTARG
           relative=true ;;
        d) folder=$OPTARG ;;
        ?) echo "Usage: cmd [-n cutoff] [-a relative_cutoff] [-d directory]"
           exit 1 ;;
    esac
done

# Create the 'filtered' directory in the specified location
mkdir -p "${folder}/filtered"

# Iterate over the files in the specified directory
for file in "${folder}"/*.txt; do
    base_name=$(basename "$file")
    output_path="${folder}/filtered/${base_name}"
    if [ "$relative" = true ]; then
        calculate_relative "$file" > "$output_path"
    else
        awk -F'\t' -v cutoff="$cutoff" 'NR==1 {print; next} {sum=0; for(i=2; i<=NF; i++) sum += $i; if(sum > cutoff) print}' "$file" > "$output_path"
    fi
done

# Generate networks in the filtered folder
# Process networks in the filtered folder
cd ${folder}/filtered/
for file in *.txt; do
    base_name=$(basename "$file")
    mkdir ${base_name}_net
    mkdir ${base_name}_net/perm/
    mkdir ${base_name}_net/pvalues/
    python /SparCC/SparCC.py ${file} --cor_file=${base_name}_net/cor_sparcc.out
    python /SparCC/MakeBootstraps.py ${file} -n 100 -t permutation_#.txt -p ${base_name}_net/perm/ 
    for f in `seq 0 99`; do
        python /SparCC/SparCC.py ${base_name}_net/perm/permutation_${f}.txt -i 100 --cor_file=${base_name}_net/pvalues/perm_cor_${f}.txt
    done
    python /SparCC/PseudoPvals.py ${base_name}_net/cor_sparcc.out ${base_name}_net/pvalues/perm_cor_#.txt 100 -o ${base_name}_net/pvals_two_sided.txt -t two_sided 
    cd ${base_name}_net
    python /SparCC/get_significant_pairs.py 

    # Change header of selected_cor.txt
    sed -i '1s/.*/Source\tTarget\tWeight\tPval/' selected_cor.txt

    cd ..
done

