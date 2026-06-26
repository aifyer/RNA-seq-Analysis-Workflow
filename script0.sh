
cd $reads_dir
while IFS= read -r line; do
       
    r1=$(ls $line*_1_*.fq | tr '\n' ',')
    r2=$(ls $line*_2_*.fq | tr '\n' ',')

    if [[ -z $r1 ]]; then
       :
    else
       cn1=$(printf $r1 | wc -c)
       reads1=$(printf $r1 | cut -b -$(( $cn1 - 1 )))  

       cn2=$(printf $r2 | wc -c)
       reads2=$(printf $r2 | cut -b -$(( $cn2 - 1 )))       
    
       STAR --runThreadN 8 \
       --genomeDir $genome_dir \
       --outFilterMultimapNmax 5 \
       --outSAMtype BAM SortedByCoordinate \
       --outBAMsortingThreadN 4 \
       --readFilesCommand zcat \
       --outFileNamePrefix $alignment_dir/$line- \
       --readFilesIn $reads1 $reads2
    
    fi

done < sample_ids.txt

cd $alignment_dir

ls *.bam > $qc_dir/bam_list.txt

while IFS= read -r line; do
    
    samtools index $line $line.bai

done < $qc_dir/bam_list.txt

while IFS= read -r line; do
    
    name=$(echo $line | cut -d- -f 1-2)
    samtools view -s 0.0001 -b $alignment_dir/$line > sampling.bam
    samtools index sampling.bam sampling.bam.bai

    geneBody_coverage.py -r $rseqc_reference -o $name -i sampling.bam
    inner_distance.py -r $rseqc_reference -k 1000 -o $name -i sampling.bam
    junction_annotation.py -r $rseqc_reference -o $name -i sampling.bam
    junction_saturation.py -r $rseqc_reference -s 50 -o $name -i sampling.bam

done < bam_list.txt

rm sampling*
