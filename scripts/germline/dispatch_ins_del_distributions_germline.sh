#!/usr/bin/bash
IFS=$'\n'

for tissue_type in $(cat ./short_topo.txt)
do

for bio_sex in $(cat ./bio_sex.txt)
do

echo "$tissue_type"$'\t'"$bio_sex"
output_dir="../../results/ins_del_fs_dists/germline/"
mkdir -p "$output_dir"

python3 ins_del_distributions_germline.py "$tissue_type" "$bio_sex" "$output_dir"
done
done

