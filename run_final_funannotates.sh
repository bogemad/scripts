#!/bin/bash

#common vars
aco_prots=/scratch/theileria/full_assemblies/aconoidasida_no_tori.proteins.fa
est=/scratch/theileria/full_assemblies/shintoku/tori_est.fasta
trin=/scratch/theileria/full_assemblies/shintoku/funannotate/training/funannotate_train.trinity-GG.fasta
out=funannotate_final_only_protein_evidence_no_tori_self_training

#ikeda
annodir=/scratch/theileria/full_assemblies/shintoku
fasta=tori_shintoku.remasked.fa.softmasked.no_organelle.fa
strain=Ikeda
isolate=Shintoku
name=MACM_
cd $annodir

funannotate predict -i $fasta \
-o $out \
-s "Theileria orientalis" \
--cpus 32 \
--isolate "$isolate" \
--strain $strain \
--name $name \
--busco_db alveolata_stramenophiles \
--organism other \
--protein_evidence $aco_prots /mnt/raid/db/uniprot_sprot.fasta \
--max_intronlen 10000 \
-w proteins:10 \
--min_training_models 100 \
--optimize_augustus \
--repeat_filter none

#chitose
annodir=/scratch/theileria/full_assemblies/chitose/annotation
fasta=chitose.trycycler.remasked.fa.softmasked.chr_fixed.fa
isolate="Fish Creek"
name=MACJ_
strain=Chitose
cd $annodir

funannotate predict -i $fasta \
-o $out \
-s "Theileria orientalis" \
--cpus 32 \
--isolate "$isolate" \
--strain $strain \
--name $name \
--busco_db alveolata_stramenophiles \
--organism other \
--protein_evidence $aco_prots /mnt/raid/db/uniprot_sprot.fasta \
--max_intronlen 10000 \
-w proteins:10 \
--min_training_models 100 \
--optimize_augustus \
--repeat_filter none

#buffeli
annodir=/scratch/theileria/full_assemblies/buffeli/annotation
fasta=buffeli.trycycler.remasked.fa.softmasked.chr_fixed.fa
isolate="Goon Nure"
name=MACK_
strain=Buffeli
cd $annodir

funannotate predict -i $fasta \
-o $out \
-s "Theileria orientalis" \
--cpus 32 \
--strain $strain \
--isolate "$isolate" \
--name $name \
--busco_db alveolata_stramenophiles \
--organism other \
--protein_evidence $aco_prots /mnt/raid/db/uniprot_sprot.fasta \
--max_intronlen 10000 \
-w proteins:10 \
--min_training_models 100 \
--optimize_augustus \
--repeat_filter none
