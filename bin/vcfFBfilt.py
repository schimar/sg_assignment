#! /usr/bin/python3
#
# This script filters a vcf file based on overall sequence coverage, number of
# non-reference reads, number of alleles, and mapping quality.
# See below for default values, and to change them, if necessary. Additionally,

# Usage: ./vcfBBfilt.py <variant file>.vcf > outfile.vcf

from sys import argv
import re

### stringency variables, edit as desired
minCoverage = 0 # minimum number of seqs; DP
minAltRds = 0 # minimum number of sequences with the alternative allele; AD
#fixed = [0.0, 1.0] # removes loci fixed for alt; RAF
mapQual = 0 # minimum mapping quality
typls = ["del", "ins", "snp"]   # other types in freebayes: [com: complex events (composite insertion and substitution events), mnp: multi-nucleotide polymorphisms]

if __name__ == "__main__":


    n_seqs_retained = int()
    with open(argv[1], 'rb') as file:
        for line in file:
            line = line.decode('utf-8')
            if line[0:1] == '#':
                #continue
                print(line.strip('\n'))
            else:
                lspl = line.strip('\n').split('\t')
                dp = int(re.findall('DP=[0-9]+', line)[0].split('=')[1])
                ac = int(re.findall('AC=[0-9]+', line)[0].split('=')[1])
                af = float(re.findall('AF=[0.0-9.0]+', line)[0].split('=')[1])
                mqm = int(re.findall('MQM=[0-9]+', line)[0].split('=')[1])
                typ = str(re.findall('TYPE=[a-z]+', line)[0].split('=')[1])
                #
                #if (af not in fixed and dp >= minCoverage and mqm >= mapQual and ac >= minAltRds):
                if (dp >= minCoverage and mqm >= mapQual and ac >= minAltRds and typ in typls):
                    print(line.strip('\n'))


        file.close()

#print '#Retained %i variable loci' % n_seqs_retained
