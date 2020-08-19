# Code to accompany the SG assignment 
Details on the actual task shall be omitted here. 


more infos here... Unix Bash

=====================================================================

<br/>

## 1) Read coverage 

First, we need to create a bed file with all regions on chromosome 19 that are not in the target region (note that we're only going to need the first 3 columns here):
```
bedtools intersect -v -a chr19.bam -b target_regions.bed -bed | cut -f1,2,3 > non_targets.bed

# and remove lines with '-1'
sed '/\-1/d' non_targets.bed > nt.bed
```

Now, we calculate the coverage for both the target & non-target regions and the coverage for chromosome 19:
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

=====================================================================

## 2) Variant calling 

We'll call variants with relaxed parameters, where 
- for *freebayes*, we'll relax the standard parameters by lowering the 'c' parameter, so that only one alternative allele has to be found to be considered, and 
- for *bbtools*, we'll decrease the rarity parameter from 1 to 0.01 (i.e. we'll consider variants with an alternative frequency of 1%) and we'll decrease the phred score for variants to be considered to 2 (where 20 would be the default).
```
freebayes/bin/freebayes -f human_g1k_v37.fa -r 19 chr19.bam -C 1 > fbvars19rlxd.vcf  

bbmap/callvariants.sh in=chr19.bam out=bbvars19rlxd2.vcf ref=human_g1k_v37.fa ploidy=2 -Xmx18g threads=8 -minscore=2.0 rarity=0.01

# further filter both vcf files:

bin\vcfFBfilt.py fbvars19rlxd.vcf > fbvars19rlxd_fltrd.vcf
bin/vcfBBfilt.py bbvars19rlxd2.vcf > bbvars19rlxd2fltrd.vcf
```
Note that the strict variant callsets were obtained by running *freebayes* without *C = 1* and *bbtools* with *minscore=20.0 and rarity=1.0*. 


<br/>

=====================================================================

## 3) Analytical performance 

For plotting of the derivations of the confusion matrices, please see the accompanied *R* script *confMatStats.r* in the *bin* folder.


