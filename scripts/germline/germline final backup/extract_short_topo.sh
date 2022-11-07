#!/usr/bin/bash

data="../../data/Germline.tsv"
cut -d $'\t' -f 59 $data | sort -u | grep -v Short_topo  > ./short_topo_germline.txt
echo "all" >> ./short_topo_germline.txt
