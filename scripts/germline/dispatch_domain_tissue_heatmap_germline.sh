#!/usr/bin/bash
IFS=$'\n'

for mutation_type in $(cat ./mutation_types.txt)
do

for bio_sex in $(cat ./bio_sex.txt)
do

echo "$mutation_type"$'\t'"$bio_sex"

python3 domain_tissue_heatmap_germline.py "$mutation_type" "$bio_sex"

done
done
