#!/bin/bash
#
# copy_from_proof_tmp.sh - script to fetch skim files from the
#                          proof work disk on /tmp/proof/...
# author: richard.t.jones at uconn.edu
# version: nov. 14, 2017

file=/tmp/proof/jonesrt/ntuple_real_skim.root
outf=/scratch/ntuple_real_skim-`hostname -s | sed 's/stat1//'`.root
if [[ -r $file ]]; then
    if [[ ! -r $outf ]]; then
		cp $file $outf
	fi
	size1=`ls -l $file | awk '{print $5}'`
    size2=`ls -l $outf | awk '{print $5}'`
    while [[ $size1 != $size2 ]]; do
		echo size mismatch between $file $size1 and $outf $size2
		cp $file $outf
		size1=`ls -l $file | awk '{print $5}'`
    	size2=`ls -l $outf | awk '{print $5}'`
	done
	echo "size matching!"
fi

file=/tmp/proof/jonesrt/mctuple_gen_skim.root
outf=/scratch/mctuple_gen_skim-`hostname -s | sed 's/nod//'`.root
if [[ -r $file ]]; then
    if [[ ! -r $outf ]]; then
		cp $file $outf
	fi
	size1=`ls -l $file | awk '{print $5}'`
    size2=`ls -l $outf | awk '{print $5}'`
    while [[ $size1 != $size2 ]]; do
		echo size mismatch between $file $size1 and $outf $size2
		cp $file $outf
		size1=`ls -l $file | awk '{print $5}'`
    	size2=`ls -l $outf | awk '{print $5}'`
	done
	echo "size matching!"
fi

file=/tmp/proof/jonesrt/mctuple_sim_skim.root
outf=/scratch/mctuple_sim_skim-`hostname -s | sed 's/nod//'`.root
if [[ -r $file ]]; then
    if [[ ! -r $outf ]]; then
		cp $file $outf
	fi
	size1=`ls -l $file | awk '{print $5}'`
	size2=`ls -l $outf | awk '{print $5}'`
    while [[ $size1 != $size2 ]]; do
		echo size mismatch between $file $size1 and $outf $size2
		cp $file $outf
		size1=`ls -l $file | awk '{print $5}'`
		size2=`ls -l $outf | awk '{print $5}'`
	done
	echo "size matching!"
fi
