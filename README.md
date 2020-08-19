# Code to accompany the SG assignment 
Details on the actual task shall be omitted here. 


more infos here...


1) Read coverage 

First, we need to create a bed file with all regions on chromosome 19 that are not in the target region:
```
bedtools intersect -v -a chr19.bam -b target_regions.bed -bed > non_targets.bed

# and remove lines with '-1'
sed '/\-1/d' non_targets.bed > nt.bed
```

2) Variant calling 




3) Analytical performance 




