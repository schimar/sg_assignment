#! /usr/bin/python3
#
# This script filters a vcf file based on overall sequence coverage, number of
# non-reference reads, number of alleles, and mapping quality.
# See below for default values, and to change them, if necessary. Additionally,

# Usage: ./vcfBBfilt.py <variant file>.vcf > outfile.vcf

from sys import argv
import re

### stringency variables, edit as desired
minCoverage = 32 # minimum number of seqs; DP
minAltRds = 2 # minimum number of sequences with the alternative allele; AD
fixed = [0.0, 1.0] # removes loci fixed for alt; RAF
mapQual = 20 # minimum mapping quality
typls = ["DEL", "INS", "SUB"]




if __name__ == "__main__":


    n_seqs_retained = int()
    with open(argv[1], 'rb') as file:
        for line in file:
            line = line.decode('utf-8')
            if line[0] == '#':
                #continue
                print(line.strip('\n'))
            else:
                lspl = line.strip('\n').split('\t')
                test = lspl[6]
                altAll = lspl[4]
                if test == 'PASS':
                    dp = int(re.findall('DP=[0-9]+', line)[0].split('=')[1])
                    ac = int(re.findall('AD=[0-9]+', line)[0].split('=')[1])
                    raf = float(re.findall('RAF=[0.0-9.0]+', line)[0].split('=')[1])
                    mqm = int(re.findall('MQM=[0-9]+', line)[0].split('=')[1])
                    typ = str(re.findall('TYP=[A-Z]+', line)[0].split('=')[1])
                    #
                    if (raf not in fixed and dp >= minCoverage and typ in typls and mqm >= mapQual and ac >= minAltRds):
                        print(line.strip('\n'))
        file.close()

