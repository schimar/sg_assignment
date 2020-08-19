# Code to accompany the SG assignment 
Details on the actual task shall be omitted here. 


more infos here... Unix Bash 


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
For plotting, the code can be found in the accompanied R script *cov_plots.r*

2) Variant calling 




3) Analytical performance 




