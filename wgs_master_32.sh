#!/bin/bash
# FASTQ to data - manuscript1 - WGS ESBL
# Started 210412 Jacky Dwiyanto
# Completed V1 210522 Jacky Dwiyanto
# Modified for 32 thread server
# Modified for SRA files

# Run this script as ./wgs_master_32.sh 
# raw sequence data (gzipped) must be placed in a folder named raw
# Script needs the pathway to the raw sequence data and assumes the paired reads were differentiated with _1 and _2.

# # Compile filenames
# for f in `ls raw_st131/*_1.fastq.gz | sed 's/raw_st131\///' `; do echo $f | sed 's/_1.fastq.gz//' ; done > sample.txt

# # QC and filtering of raw reads
# mkdir -p fastp_out
# while read line; do fastp --thread 32 -i raw/${line}_1.fastq.gz -I raw/${line}_2.fastq.gz -o fastp_out/${line}_1.fastq.gz -O fastp_out/${line}_2.fastq.gz; done < sample.txt

# SRST2 analysis
# source ~/miniconda3/bin/activate srst2 # verify srst2 environment
# 
# mkdir -p srst2_out
# while read line; do srst2 --threads 32 --input_pe fastp_out/${line}_1.fastq.gz fastp_out/${line}_2.fastq.gz --output srst2_out/$line --log --mlst_db ~/analysis_data/reference_wgs/Escherichia_coli#1.fasta --mlst_definitions ~/analysis_data/reference_wgs/profiles_csv --mlst_delimiter _ ; done < sample.txt
# # while read line; do srst2 --threads 32 --input_pe fastp_out/${line}_1.fastq.gz fastp_out/${line}_2.fastq.gz --output srst2_out/$line --log --mlst_db ~/analysis_data/reference_wgs/Escherichia_coli#1.fasta --mlst_definitions ~/analysis_data/reference_wgs/profiles_csv --mlst_delimiter _ --gene_db ~/analysis_data/reference_wgs/ARGannot_r3.fasta ~/analysis_data/reference_wgs/PlasmidFinder.fasta ~/analysis_data/reference_wgs/Escherichia_VF_clustered.fasta ~/analysis_data/reference_wgs/Plasmid18Replicons.fasta; done < sample.txt
# compile SRST2 output
# 
# srst2 --prev_output srst2_out/*__mlst__* --output srst2_compiled_out/srst2_mlst_compiled
# srst2 --prev_output srst2_out/*_genes__ARGannot_r3__results.txt --output srst2_compiled_out/srst2_argannot_compiled
# srst2 --prev_output srst2_out/*_genes__Escherichia_VF_clustered__results.txt --output srst2_compiled_out/srst2_VF_compiled
# srst2 --prev_output srst2_out/*_genes__Plasmid18Replicons__results.txt --output srst2_compiled_out/srst2_plasmidreplicon_compiled
# srst2 --prev_output srst2_out/*_genes__PlasmidFinder__results.txt --output srst2_compiled_out/srst2_plasmidfinder_compiled

# Assembly with unicycler
source ~/miniconda3/bin/activate unicycler
mkdir -p unicycler_out
while read line; do unicycler -t 32 -1 ./fastp_out/${line}_1.fastq.gz -2 ./fastp_out/${line}_2.fastq.gz -o unicycler_out/$line/; done < sample.txt

# # Prokka and roary analysis
source ~/miniconda3/bin/activate prokka
while read line; do prokka --cpus 0 --outdir ./prokka_out/${line}/ ./unicycler_out/${line}/assembly.fasta; done < sample.txt

mkdir -p prokka_out/prokka_compiled # make directory to compiled gff files
while read line; do cp prokka_out/${line}/*gff prokka_out/prokka_compiled/prokka_${line}.gff; done < sample.txt # move gff to a single folder

# # run roary
roary -f roary_out -p 32 -e -n -v -r prokka_out/prokka_compiled/*gff

# # generate tree
FastTree -nt -gtr roary_out/core_gene_alignment.aln > roary_out/my_tree.newick

# # generate plot
python roary_plots.py roary_out/my_tree.newick roary_out/gene_presence_absence.csv

