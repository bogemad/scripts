#!/bin/bash

raw_key=$1
original_annotation_gb=$2
manual_annotation_dir=$3
species="Theileria orientalis"
isolate="$4"
strain="$5"
prefix="${species// /_}_${strain}"
sbt=$6
outdir=$7
final_key=$outdir/$(basename $raw_key)

rm -r $outdir
mkdir -p $outdir

blanks_script=/mnt/raid/scratch_tmp/theileria/scripts/fill_in_locus_tag_blanks.py
fix_ltag_script=/mnt/raid/scratch_tmp/theileria/scripts/fix_locus_tags.py
run_tbl2asn_script=/mnt/raid/scratch_tmp/theileria/scripts/generate_sqn.sh

python $blanks_script $raw_key $manual_annotation_dir/$prefix.gbk $original_annotation_gb $final_key
python $fix_ltag_script $final_key $original_annotation_gb $manual_annotation_dir/$prefix.tbl $outdir/$prefix.tbl
cp -av $manual_annotation_dir/$prefix.scaffolds.fa $outdir/$prefix.fsa
cp -av $sbt $outdir/
sbt=$outdir/$(basename $sbt)
bash $run_tbl2asn_script $outdir $sbt $species $isolate $strain
mv $outdir/$prefix.tbl $outdir/$prefix.tbl.bak
funannotate gbk2parts -g $outdir/$prefix.gbk -o $outdir/$prefix
mv $outdir/$prefix.tbl.bak $outdir/$prefix.tbl
