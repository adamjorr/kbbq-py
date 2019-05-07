#!/usr/bin/env python3

import kbbq.compare_reads as cr
import kbbq.benchmark as bm
import argparse
import sys

#helper commands that pass arguments to the proper functions

def recalibrate(args):
    pass

def benchmark(args):
    bm.benchmark(bamfile = args.bam, fafile = args.reference,
        vcffile = args.vcf, fastq = args.fastq, label = args.label)

def plot(args):
    pass

def main():
    parser = argparse.ArgumentParser(description = 'K-mer Based Base Quality score recalibration')
    subparsers = parser.add_subparsers(title='command', help="valid commands"), #choices=['recalibrate','benchmark','plot'])
    
    #recalibrate command TODO
    recalibrate_parser = subparsers.add_parser('recalibrate', description = 'Recalibrate a BAM or FASTQ file')
    
    #benchmark command
    benchmark_parser = subparsers.add_parser('benchmark', description = 'Benchmark a BAM or FASTQ file using a truth set')
    benchmark_parser.add_argument('-b', '--bam', required = True,
        help = 'Truth set BAM file. Differences from the reference at nonvariable sites will be interpreted as errors.')
    benchmark_parser.add_argument('-r', '--reference', required = True, help = 'FASTA file containing the reference genome')
    benchmark_parser.add_argument('-v', '--vcf', required = True, help = 'VCF file containing variable sites')
    benchmark_parser.add_argument('-f', '--fastq', default = None, required = False, help = 'fastq file to benchmark')
    benchmark_parser.add_argument('-l', '--label', default = None, required = False, help = 'label to use for label column')
    benchmark_parser.set_defaults(command=benchmark)

    #plot command
    plot_parser = subparsers.add_parser('plot', description = 'Plot data output from the benchmark command')

    #parse args
    args = parser.parse_args()
    args.command(args)

if __name__ == '__main__':
    main()
