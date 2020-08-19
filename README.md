# Code to accompany the SG assignment 
Details on the actual task shall be omitted here. 

The following software is required to make this work:  
	- UNIX Bash  
	- Python 3.8  
	- R (tested on v.3.6.3)  
	- bedtools (tested on v.2.29.2)  
	- bbtools v.38.86  
	- freebayes v.1.3.2  
	- samtools (with v.>=1.10)  

The code was run on a Ubuntu 20.04 OS, and all of the below code was run in Bash. 
Please note that the path of respective scripts is not taken into account below. You'd have to append the *bin/* folder to make it work. Also, the files given in the assignment, as well as the human_g1k_v37.fa assembly (downloaded from the 1000genome project) need to be in your working directory. 

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

The vcf files for all four callsets can be found in the subfolder [*vars*](https://github.com/schimar/sg_assignment/tree/master/vars), with  [*freebayes*](http://github.com/schimar/sg_assignment/vars/fbvars19rlxd_fltrd.tar.gz) and [*bbtools*](http://github.com/schimar/sg_assignment/vars/bbvars19rlxd_fltrd.tar.gz) vcf files with the relaxed parameters and the remaining two vcf files are called with stricter parameters. 
We'll call variants with relaxed parameters, where 
- for *freebayes*, we'll relax the standard parameters by lowering the 'c' parameter, so that only one alternative allele has to be found to be considered, and 
- for *bbtools*, we'll decrease the rarity parameter from 1 to 0.01 (i.e. we'll consider variants with an alternative frequency of 1%) and we'll decrease the phred score for variants to be considered to 2 (where 20 would be the default).
```
freebayes/bin/freebayes -f human_g1k_v37.fa -r 19 chr19.bam -C 1 > fbvars19rlxd.vcf  

bbmap/callvariants.sh in=chr19.bam out=bbvars19rlxd2.vcf ref=human_g1k_v37.fa ploidy=2 -Xmx18g threads=8 -minscore=2.0 rarity=0.01

# further filter both vcf files (see Supporting Information in the report, as well as the respective *python* filtering scripts in the *bin* folder for details on filtering criteria):

bin/vcfFBfilt.py fbvars19rlxd.vcf > fbvars19rlxd_fltrd.vcf
bin/vcfBBfilt.py bbvars19rlxd2.vcf > bbvars19rlxd2fltrd.vcf
```
Note that the strict variant callsets were obtained by running *freebayes* without *C = 1* and *bbtools* with *minscore=20.0 and rarity=1.0*. 

<br/>

### Tables 1 and S1
To get the numbers of SNPs, DELs and INS for ground_truth.vcf, we can run:
```
grep -v '#' ground_truth.vcf | awk -f bin/nTypGT.awk
```
Similarly, for the *freebayes* and *bbtools* sets, we can do the following: 

```
# fbvars
egrep -o 'TYPE=[a-z]{3}' fbvars19rlxd_fltrd.vcf | grep -v '#' | wc -l
egrep -o 'TYPE=[a-z]{3}' fbvars19rlxd_fltrd.vcf | cut -f2 -d'=' | sort | uniq -c

# bbvars:
egrep -o 'TYP=[A-Z]{3}' bbvars19fltrd.vcf | cut -f2 -d'=' | sort | uniq -c
egrep -o 'TYP=[A-Z]{3}' bbvars19rlxd2fltrd.vcf | cut -f2 -d'=' | sort | uniq -c
```

And lastly, we count how many variants are present in our target-regions:
```
# fbvars
vcfInBed.py fbvars19rlxd_fltrd.vcf ../target_regions.bed | wc -l
vcfInBed.py fbvars19fltrd.vcf ../target_regions.bed | wc -l

# bbvars
vcfInBed.py bbvars19fltrd.vcf ../target_regions.bed | wc -l
vcfInBed.py bbvars19rlxd2fltrd.vcf ../target_regions.bed | wc -l

# ground_truth.vcf 
vcfInBed.py ground_truth.vcf ../target_regions.bed | wc -l
```

<br/>

=====================================================================

## 3) Analytical performance 

We first calculate the confusion matrix for the respective variant callset. Here, the ground-truth set always acts as the "actual" data, with the *freebayes* and *bbtools* variants as predicted variant set. 

```

### confusion matrices: 
for i in *bvars*fltrd.vcf; do  
	echo $i >> confMat.txt
	python3 ../shd/getConfusionMatrix.py ground_truth.vcf $i  >> confMat.txt
	printf "\n---------------------------------------------------------\n" >> confMat.txt
done

### this would give us the following (e.g. for bbvars19rlxd2fltrd.vcf)

#		 bbvars19rlxd2fltrd.vcf
#		 +------------+-----------------+----------------+
#		 |            |   predicted yes |   predicted no |
#		 +============+=================+================+
#		 | actual yes |             966 |          17647 |
#		 +------------+-----------------+----------------+
#		 | actual no  |           19118 |       14701826 |
#		 +------------+-----------------+----------------+ 
#		 
#		 Sensitivity = 5
#		 Precision = 4
#		 Specificity = 99


### derivations from confusion matrix (we pretty much do the same as above, but here, we only grep the derivations ffrom the output, to write them into a file): 

for i in *bvars*fltrd.vcf; do  
	echo $i >> confMatDerivat.txt
	python3 ../shd/getConfusionMatrix.py ground_truth.vcf $i | grep 'Sensitivity\|Precision\|Specificity' >> confMatDerivat.txt
	printf "\n---------------------------------------------------------\n" >> confMatDerivat.txt
done

```

To plot the derivations of the confusion matrices (i.e. sensitivity, precision and specificity), please see the accompanied *R* script *confMatStats.r* in the *bin* folder.


