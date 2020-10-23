#! /usr/bin/env python

import sys

vcf = file(sys.argv[1])

#python program for filtering VCF files

print_flag = False

for line in vcf:
    if "#CHROM" in line:
        print_flag = True
    if print_flag:
        print(line)

