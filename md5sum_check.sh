#!/usr/bin/env bash


for i in $(ls -d /staging/wl/caim/RNA_seq_KO/*/*/); do
	cd $i
	file=$(ls *.txt)
	md5sum -c $file
done
