# coding=utf-8

import os
import sys
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('--infile', help='input file', action='store', dest='infile', default='', required=True)
parser.add_argument('--outfile', help='', action='store', dest='outfile', default='', required=True)
args = parser.parse_args()


def main():
    out = open(args.outfile,'w')
    header = ['chr', 'locs', 'Iseq', 'svlen', 'all_reads', 'alt_reads', 'predict']
    out.write('\t'.join(header) + '\n')
    vcf = open(args.infile, 'r')
    for line in vcf:
        if line.startswith('Chr'):
            continue
        array = line.strip().split('\t')
        chr_t = array[0]
        locs = array[1]
        Iseq = array[4]
        func_ref = array[5]
        #if func_ref == 'intronic':
        #    continue
        otherinfo = array[-1]
        otherinfos = otherinfo.split(':')[1]
        otherinfos1 = otherinfos.split(',')
        allreads = float(otherinfos1[0]) + float(otherinfos1[1])
        mutreads = otherinfos1[1]
        info = array[-3]  # 'SVLEN=24;AF=0.01655;SVTYPE=FLT3-ITD'
        if re.search(r"SVLEN=(?P<length>\d+);AF=(?P<vaf>\d+\.\d+)", info):
            match = re.search(r"SVLEN=(?P<length>\d+);AF=(?P<vaf>\d+\.\d+)", info)
            vaf = match.group('vaf')
            if float(vaf) < 0.01:
                vaf = 0.01
            length = match.group('length')
            out.write('\t'.join([chr_t, str(locs), Iseq, str(length), str(allreads), str(mutreads), str(vaf)]) + '\n')


if __name__ == '__main__':
    main()
