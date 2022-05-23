#!/bin/bash

# Time the script
# Note: SECONDS is a bash builtin variable that gives the number of seconds since the script was started
SECONDS=0

# Following instructions on https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_mr#rat_6.0.0

#Download fasta
wget http://ftp.ensembl.org/pub/release-105/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz
gunzip Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz
 
#Download GTF
wget http://ftp.ensembl.org/pub/release-105/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.105.gtf.gz
gunzip Rattus_norvegicus.mRatBN7.2.105.gtf.gz

#Filter GTF
cellranger mkgtf \
  Rattus_norvegicus.mRatBN7.2.105.gtf Rattus_norvegicus.mRatBN7.2.105.filtered.gtf \
  --attribute=gene_biotype:protein_coding \
  --attribute=gene_biotype:lincRNA \
  --attribute=gene_biotype:antisense \
  --attribute=gene_biotype:IG_LV_gene \
  --attribute=gene_biotype:IG_V_gene \
  --attribute=gene_biotype:IG_V_pseudogene \
  --attribute=gene_biotype:IG_D_gene \
  --attribute=gene_biotype:IG_J_gene \
  --attribute=gene_biotype:IG_J_pseudogene \
  --attribute=gene_biotype:IG_C_gene \
  --attribute=gene_biotype:IG_C_pseudogene \
  --attribute=gene_biotype:TR_V_gene \
  --attribute=gene_biotype:TR_V_pseudogene \
  --attribute=gene_biotype:TR_D_gene \
  --attribute=gene_biotype:TR_J_gene \
  --attribute=gene_biotype:TR_J_pseudogene \
  --attribute=gene_biotype:TR_C_gene

#Run mkref
cellranger mkref \
  --genome=mRatBN7 \
  --fasta=Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa \
  --genes=Rattus_norvegicus.mRatBN7.2.105.filtered.gtf \
  --ref-version=1.0.0

echo "Rat reference built in $(( $SECONDS / 60 )) minutes."