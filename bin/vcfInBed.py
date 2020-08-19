#! /bin/python3

# this script reads through a bed file with target regions and a vcf file with variants, and outputs only the variants that lie in the respective target regions.

# Usage: ./vcfInBed.py <file>.vcf <regions>.bed

import numpy as np
from sys import argv



with open(argv[2], 'rb') as bed:
    bedpos = np.array(list())#np.empty(shape=1, dtype=np.int8)
    for line in bed:
        line = line.decode('utf-8')
        chrom, start, stop = map(int, line.strip('\n').split('\t'))
        bedpos = np.append(bedpos, np.arange(start+1, stop+2))  # bed file is 0-indexed, while vcf is not
bed.close()


with open(argv[1], 'rb') as vcf:
    for line in vcf:
        line = line.decode('utf-8')
        if line[0:1] == '#':
            #continue
            print(line.strip('\n'))
        else:
            vcfpos = int(line.split('\t')[1])
            if vcfpos in bedpos:
                print(line.strip('\n'))
vcf.close()






