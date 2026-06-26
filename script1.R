
BiocManager::install("Rsubread")
library("Rsubread")

d=scan(file="bam_files.txt",what="")
bam_files=dput(d)

anno_gtf="annotation.gtf"

fCount=featureCounts(files=bam_files,annot.ext=anno_gtf,
                     isGTFAnnotationFile=TRUE,GTF.attrType="gene_id", nthreads=4, 
                     isPairedEnd=TRUE, countReadPairs=TRUE,countMultiMappingReads=TRUE, fraction=TRUE)

save(fCount,file="feature_counts.rda")

write.table(fCount$stat,file="featureCounts.summary",col.names=T,row.names=F,sep="\t") 
