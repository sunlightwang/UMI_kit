# UMI_kit

parsing RNA-seq data with Unique Molecular Identifiers (UMI)


## Usage
1. sequencing adapter trimming

2. UMI barcode, following G, and 3' poly-(A) trimming
```UMI_trim.pl yoursample.fastq.gz yoursample```

3. read mapping to genome

4. convert bam to bed
```bamToBed -bed12 -i yoursample.bam | awk -vOFS='\t' '{split($4,a,/=/); $4=a[2]; print $0}' | gzip > yoursample.bed.gz```

5. UMI collapse
```UMI_collapse.pl yoursample.bed.gz yoursample.UMI_collapsed.bed.gz```

6. UMI deduplicates
```UMI_dedup.pl yoursample.UMI_collapsed.bed.gz yoursample.UMI_dedup.bed.gz```

## Contact

xi.wang (at) dkfz-heidelberg.de

