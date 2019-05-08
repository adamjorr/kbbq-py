#!/usr/bin/env python3

import kbbq
import compare_reads as cr
import benchmark as bm
import argparse
import sys

#helper commands that pass arguments to the proper functions

def recalibrate(args):
    pass

def benchmark(args):
    bm.benchmark(bamfile = args.bam, fafile = args.reference,
        vcffile = args.vcf, fastqfile = args.fastq,
        label = args.label, use_oq = args.use_oq)

def plot(args):
    pass

def main():
    parser = argparse.ArgumentParser(description = 'K-mer Based Base Quality score recalibration')
    parser.add_argument('-v', '--version', action = 'version', version = kbbq.__version__)
    parser.set_defaults(command = lambda args: parser.print_help())
    subparsers = parser.add_subparsers(title='command', description="valid commands")
    
    #recalibrate command TODO
    recalibrate_parser = subparsers.add_parser('recalibrate', description = 'Recalibrate a BAM or FASTQ file')
    
    #benchmark command
    benchmark_parser = subparsers.add_parser('benchmark', description = 'Benchmark a BAM or FASTQ file using a truth set')
    benchmark_reqd = benchmark_parser.add_argument_group(title = 'required arguments')
    benchmark_reqd.add_argument('-b', '--bam', required = True,
        help = 'Truth set BAM file. Differences from the reference at nonvariable sites will be interpreted as errors.')
    benchmark_reqd.add_argument('-r', '--reference', required = True, help = 'FASTA file containing the reference genome')
    benchmark_reqd.add_argument('-v', '--vcf', required = True, help = 'VCF file containing variable sites')
    benchmark_parser.add_argument('-f', '--fastq', default = None, required = False, help = 'fastq file to benchmark')
    benchmark_parser.add_argument('-l', '--label', default = None, required = False, help = 'label to use for label column')
    benchmark_parser.add_argument('-u', '--use-oq', action = 'store_true', help = 'Use the OQ tag when benchmarking a BAM file. Does nothing if a fastq file is provided.')
    benchmark_parser.set_defaults(command=benchmark)

    #plot command
    plot_parser = subparsers.add_parser('plot', description = 'Plot data output from the benchmark command')

    #parse args
    args = parser.parse_args()
    args.command(args)

if __name__ == '__main__':
    main()
