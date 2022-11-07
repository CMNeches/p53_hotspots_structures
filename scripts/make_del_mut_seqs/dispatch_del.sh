#!/usr/bin/bash

#python3 make_SNV_mut_seqs.py
chars=( {a..z} )
n=15
for ((i=0; i<n; i++))
do
	name="${chars[i]}"
	isoform_name="isoform_""$name"
	echo "$isoform_name"
	python3 make_del_mut_seqs.py "$isoform_name"
done

#write a bash script that loops through the letters a through o, and runs the python script
