# kbbq

[![codecov](https://codecov.io/gh/adamjorr/kbbq-py/branch/master/graph/badge.svg)](https://codecov.io/gh/adamjorr/kbbq-py)
[![CircleCI](https://circleci.com/gh/adamjorr/kbbq-py.svg?style=svg)](https://circleci.com/gh/adamjorr/kbbq-py)
[![Documentation Status](https://readthedocs.org/projects/kbbq/badge/?version=latest)](https://kbbq.readthedocs.io/en/latest/?badge=latest)

k-mer based base quality recalibration

This is a python prototype of kbbq; the C++ version can be found at https://www.github.com/adamjorr/kbbq

The full documentation can be found here: https://kbbq.readthedocs.io/

## Installation

kbbq can be installed via pip or your favorite pip-compatible program,
like [Pipenv](https://docs.pipenv.org/).

To install the latest release with pip::

	pip install git+https://github.com/adamjorr/kbbq.git

## Recalibrating Reads

kbbq uses a FASTQ file and a corrected FASTQ file to produce a recalibrated
FASTQ file. The output file has the original read sequences but modified quality
scores.

If the reads are in a file called `reads.fq` and `reads.corrected.fq`, you can
recalibrate `reads.fq` to produce a `reads.recalibrated.fq` file
with the `recalibrate` command like so:

```
kbbq recalibrate -f reads.fq reads.cor.fq > reads.recalibrated.fq
```

The FASTQ files must be interleaved. If you have read group information in the
read names you can use the ``--infer-rg`` flag to infer which read group each
read belongs to. You can read more about the naming conventions required to
infer pairing and read group data here: [FASTQ Input Help](https://kbbq.readthedocs.io/en/latest/cli/fastq_input.html#fastq-input-help) .
Without ``--infer-rg``, it is assumed all reads in the input file come from a single
read group.

Benchmark Quality Scores
------------------------
kbbq can benchmark quality scores given a BAM alignment, a fasta reference,
a VCF, and optionally a BED file. Bases in the alignment that:

 - are in a variable site or
 - are outside the regions specified in the BED file or
 - are soft clipped

are skipped. If a base is not skipped but does not match the reference, it is
assumed to be an error.

The output file will be a tab separated file with 4 columns:

1. assigned quality scores that appear in the data set
2. the actual quality score of bases that were assigned the score
3. an optionally-provided label (modified with ``--label``)
4. the number of bases that were assigned that quality score in the dataset

An example invokation is:

```
kbbq benchmark -b alignment.bam -r reference.fa -v variants.vcf -d good-sites.bed > calibration.tsv
```

and an example output is:

```
2       10      alignment.bam    6460879
3       3       alignment.bam    170777
4       3       alignment.bam    173821
5       5       alignment.bam    246709
6       5       alignment.bam    304227
7       6       alignment.bam    463211
8       7       alignment.bam    416902
9       8       alignment.bam    491355
```

A FASTQ file of reads can be provided with the `-f` flag if the reads in the file are properly
named such that they can be found in the corresponding BAM file. The BAM file will be used to determine
which bases are errors. For more detail on the required naming scheme, see [FASTQ Input Help](https://kbbq.readthedocs.io/en/latest/cli/fastq_input.html#fastq-input-help) .

Plotting Benchmarked Scores
---------------------------

The output of the `benchmark` command can be plotted with the `plot` command.
This makes a basic plot; the output type will be inferred by matplotlib based on the
file extension. To plot the calibration of our file created above:

```
kbbq plot -o calibration.pdf calibration.tsv
```

You can use `-t sample-size` to plot the number of bases instead of the quality
score calibration::

```
kbbq plot -t sample-size -o sample-size.pdf calibration.tsv
```

For further reading on program features, check out https://kbbq.readthedocs.io/en/latest/cli/cli.html.

