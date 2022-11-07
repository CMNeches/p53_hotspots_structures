#!/usr/bin/bash
IFS=$'\n'

for tissue_type in $(cat ./short_topo.txt)
do

for mutation_type in $(cat ./mutation_types.txt)
do

for bio_sex in $(cat ./bio_sex.txt)
do

echo "$tissue_type"$'\t'"$mutation_type"$'\t'"$bio_sex"
output_dir="../../results/hotspots/somatic/"${tissue_type}"/"${mutation_type}"/"${bio_sex}""
mkdir -p "$output_dir"

python3 mutation_distribution_graphs_somatic.py "$tissue_type" "$mutation_type" "$bio_sex" "$output_dir"

done
done
done
