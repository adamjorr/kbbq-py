#!/usr/bin/env python3

import kbbq
import kbbq.compare_reads as cr
import kbbq.benchmark as bm
import kbbq.recalibrate as re
import kbbq.plot
import argparse
import sys

#helper commands that pass arguments to the proper functions

def recalibrate(args):
    re.recalibrate(bam = args.bam, fastq = args.fastq, infer_rg = args.infer_rg,
    use_oq = args.use_oq, set_oq = args.set_oq, gatkreport = args.gatkreport)

def benchmark(args):
    bm.benchmark(bamfile = args.bam, fafile = args.reference,
        vcffile = args.vcf, fastqfile = args.fastq,
        label = args.label, use_oq = args.use_oq,
        bedfh = args.bedfile)

def plot(args):
    kbbq.plot.plot_benchmark(fhin = args.file, outfile = args.outfile, plottype = args.type)

def main():
    parser = argparse.ArgumentParser(description = 'K-mer Based Base Quality score recalibration')
    parser.add_argument('-v', '--version', action = 'version', version = kbbq.__version__)
    parser.set_defaults(command = lambda args: parser.print_help())
    subparsers = parser.add_subparsers(title='command', description="valid commands")
    
    #reused strings
    oq_help = 'Use the OQ tag to get quality scores when working with a BAM file. Does nothing if a fastq file is provided.'

    #recalibrate command
    recalibrate_parser = subparsers.add_parser('recalibrate', description = 'Recalibrate a BAM or FASTQ file')

    #currently we're making inputs a required input and only allowing 1 input.
    #In the future, we should allow as many inputs as desired, and if
    #there are no inputs provided read from STDIN and try to infer the filetype.
    #They should be positional arguments too, probably.
    recalibrate_input = recalibrate_parser.add_mutually_exclusive_group(required = True)
    recalibrate_input.add_argument('-b','--bam', help = 'BAM to recalibrate')
    recalibrate_input.add_argument('-f', '--fastq', nargs=2, help = 'FASTQ file to recalibrate and a corrected version from your favorite error corrector.')

    recalibrate_parser.add_argument('-u','--use-oq', action = 'store_true', help = oq_help)
    recalibrate_parser.add_argument('-s','--set-oq', action = 'store_true', help = 'Set the \'OQ\' flag prior to recalibration. Only works when producing BAM output.')
    recalibrate_parser.add_argument('-g','--gatkreport',
        help = 'If the given path points to an existing GATK report, \
        load the model from the report instead of calculating it. \
        If the file doesn\'t exist, save the calculated model to the given path.')
    recalibrate_parser.add_argument('--infer-rg', action = 'store_true',
        help = 'Attempt to infer the read group from a FASTQ read. Only works with FASTQ input. \
        The default behavior is to treat each input FASTQ file as its own read group.')
    #TODO: method, model, lighter options, prefix, output, gatkreport
    recalibrate_parser.set_defaults(command=recalibrate)

    #benchmark command
    benchmark_parser = subparsers.add_parser('benchmark', description = 'Benchmark a BAM or FASTQ file using a truth set')
    benchmark_reqd = benchmark_parser.add_argument_group(title = 'required arguments')
    benchmark_reqd.add_argument('-b', '--bam', required = True,
        help = 'Truth set BAM file. Differences from the reference at nonvariable sites will be interpreted as errors.')
    benchmark_reqd.add_argument('-r', '--reference', required = True, help = 'FASTA file containing the reference genome')
    benchmark_reqd.add_argument('-v', '--vcf', required = True, help = 'VCF file containing variable sites')
    benchmark_parser.add_argument('-f', '--fastq', default = None, required = False, help = 'fastq file to benchmark')
    benchmark_parser.add_argument('-l', '--label', default = None, required = False, help = 'label to use for label column')
    benchmark_parser.add_argument('-u', '--use-oq', action = 'store_true', help = oq_help)
    benchmark_parser.add_argument('-d', '--bedfile', type = argparse.FileType('r'),
        help = 'BED file of confident regions. Sites outside the given regions will be skipped.')
    benchmark_parser.set_defaults(command=benchmark)

    #plot command
    plot_parser = subparsers.add_parser('plot', description = 'Plot data output from the benchmark command')
    plot_reqd = plot_parser.add_argument_group(title = 'required arguments')
    plot_parser.add_argument('-t', '--type', default = 'calibration', choices = ['calibration', 'sample-size'], help = 'Type of plot to produce')
    plot_parser.add_argument('file', nargs = '?', type = argparse.FileType('r'), default = sys.stdin, help = 'Input file')
    plot_reqd.add_argument('-o', '--outfile', required = True, help = 'file name to save plot as')
    plot_parser.set_defaults(command=plot)

    #parse args
    args = parser.parse_args()
    args.command(args)

if __name__ == '__main__':
    main()
