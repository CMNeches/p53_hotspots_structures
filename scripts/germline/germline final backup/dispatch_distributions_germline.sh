#!/usr/bin/bash
IFS=$'\n'

for tissue_type in $(cat ./short_topo_germline.txt)
do

for mutation_type in $(cat ./mutation_types.txt)
do

for bio_sex in $(cat ./bio_sex.txt)
do

echo "$tissue_type"$'\t'"$mutation_type"$'\t'"$bio_sex"
output_dir="../../results/figures_and_statistics/germline/distributions/"${tissue_type}"/"${mutation_type}"/"${bio_sex}""
mkdir -p "$output_dir"

python3 distributions_germline.py "$tissue_type" "$mutation_type" "$bio_sex" "$output_dir"
done
done
done

