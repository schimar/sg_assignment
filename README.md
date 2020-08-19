# Code to accompany the SG assignment 
Details on the actual task shall be omitted here. 


more infos here... Unix Bash 

<br/>

1) Read coverage 

First, we need to create a bed file with all regions on chromosome 19 that are not in the target region (note that we're only going to need the first 3 columns here):
```
bedtools intersect -v -a chr19.bam -b target_regions.bed -bed | cut -f1,2,3 > non_targets.bed

# and remove lines with '-1'
sed '/\-1/d' non_targets.bed > nt.bed
```

Then, we calculate the coverage for both the target & non-target regions and the coverage for chromosome 19:
```
bedtools coverage -a target_regions.bed -b chr19.bam -mean > MeanCoverageBEDtarget.bedgraph
bedtools coverage -a nt.bed -b chr19.bam -mean > MeanCoverageBEDnt.bedgraph
bedtools genomecov -ibam chr19.bam -bga -split | grep '^19' > CoverageTotal.bedgraph

```

Next, we calculate the coverage histograms for target, non-target and total regions (for this section of chr 19) and pull out the cumulative portion of the output files:
```
bedtools coverage -hist -a target_regions.bed -b chr19.bam | grep  '^all' > target_hist_all.cov
bedtools coverage -hist -a nt.bed -b chr19.bam | grep '^all' > nt_hist_all.cov
bedtools coverage -hist -abam chr19.bam -b target_regions.bed nt.bed | grep '^all' > all_hist_all.cov
```

Further, we calculate the nucleotide composition, to get the GC content in our region:
```
bedtools nuc -fi human_g1k_v37.fa -bed target_regions.bed > targets.nuc
bedtools nuc -fi human_g1k_v37.fa -bed nt.bed > nt.nuc 


# Let's quickly remove the '#' in the first line of both files:
for i in *.nuc; do
	id=${i%.nuc}
    sed 's/#//g' $i > ${id}_tmp.nuc
	mv ${id}_tmp.nuc}  
done

```

For the calculation of summary statistics and plotting, please see the accompanied *R* script *cov_plots.r* in the *bin* folder.

<br/>

2) Variant calling 

fb, bbtools & filtering
 variant calling settings and filtering (word so that it fits with the "make sure that you clarify" - first point
 details about matching of variants for ex. 3

 matching of variants with bed file !!! 

<br/>

3) Analytical performance 

For plotting of the derivations of the confusion matrices, please see the accompanied *R* script *confMatStats.r* in the *bin* folder.


