Example Recalibration Workflow
==============================
.. correcting reads and recalibrating

This tutorial will guide you through a full recalibration workflow
with ``kbbq`` using an example data set. The data set comes from 
[Li_2018]_. It is a variant-calling benchmarking data set of sequenced
human hydatidiform moles.

.. [Li_2018] Li H, Bloom JM, Farjoun Y, Fleharty M, Gauthier L, Neale B, MacArthur D (2018) A synthetic-diploid benchmark for accurate variant-calling evaluation. Nat Methods, 15:595-597. [PMID:30013044]

Software Requirements
---------------------

.. highlight:: bash

This tutorial will require slightly more software than ``kbbq`` alone.
These are somewhat common tools that appear in many bioinformatic workflows
and will be used to download and manipulate the data. We'll need:

#. samtools
#. bcftools
#. seqtk
#. lighter

Samtools and bcftools follow the download, ``./configure``,
``make``, ``make install`` paradigm. seqtk does not require ``./configure``,
and lighter will require manual installation to a directory on your PATH.

.. note::

	These installation instructions will install the programs to the
	``~/.local`` directory, which is where ``pip`` installs programs
	with the ``--user`` flag. If you want to install somewhere else,
	modify the appropriate ``--prefix``, ``BINDIR`` or ``install``
	argument.

To download and install samtools::

	wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
	bunzip samtools-1.9.tar.bz2
	./configure --prefix=~/.local
	make
	make install

To install bcftools::

	wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
	bunzip bcftools-1.9.tar.bz2
	./configure --prefix=~/.local
	make
	make install

To download and install seqtk::

	git clone https://github.com/lh3/seqtk.git
	cd seqtk
	make
	make install BINDIR=~/.local/bin

To download and install lighter::

	git clone https://github.com/mourisl/Lighter.git
	cd Lighter
	make
	install lighter ~/.local/bin #install is like cp

Finally, you'll need to have ``kbbq`` installed.
The fast way to do this is::

	pip install --user git+https://github.com/adamjorr/kbbq.git

Follow the instructions in :doc:`installation` for more information.

.. note::

	These instructions install to the ``.local/bin`` directory in your home directory.
	This follows the convention ``pip install --user`` uses, and is probably already on your PATH.
	If it isn't, edit ``~/.profile`` and add the following line::

		export PATH=${PATH}:~/.local/bin

	Then run::

		. ~/.profile

	to reload the variable.

Downloading the Data
--------------------

Start by making a new directory we can use to store the files
we'll create in this tutorial::

	mkdir kbbq_tutorial
	cd kbbq_tutorial

As an example dataset, we'll use the [Li_2018]_ CHM1/13 data set.
This is a variant benchmarking dataset that includes:

#. bam alignment
#. BED file designating confident regions
#. VCF file containing variable sites

.. note::

	The only data required to recalibrate are the reads (in this case,
	from the BAM file.) We use the BED to subset the BAM to a reasonable
	size. If you aren't going to benchmark, download and subset your bam
	and skip directly to :ref:`correcting_reads`

We'll start by downloading the CHM1/13 evaluation kit and subsetting it
to a more reasonable size for a tutorial. The following lines download
the tar archive, extract the BED file, uncompress it, and take the first
25 lines::

	wget https://github.com/lh3/CHM-eval/releases/download/v0.5/CHM-evalkit-20180222.tar
	tar -xnf CHM-evalkit-20180222.tar --to-stdout CHM-eval.kit/full.37m.bed.gz | \
	zcat | \
	head -n 25 > confident.bed

This will extract the first 10 confident regions to a file called confident.bed.
We will limit our analyses to these 10 regions.

The evaluation kit also includes a VCF of variable sites. We will use ``bcftools view``
to subset the whole file to just the confident regions we picked earlier::

	tar -xnf CHM-evalkit-20180222.tar --to-stdout CHM-eval.kit/full.37m.vcf.gz | \
	bcftools view -T confident.bed -Oz -o confident.vcf.gz

.. note::

	A reference is **only** required for the :ref:`benchmark` command.
	If you're not interested in benchmarking, feel free to skip the reference
	download.

We also need the reference for the benchmark command.
We only need the chromosomes specified in our region
file. To avoid downloading the whole file, we can download the reference index
first and then just specify the chromosomes that are present in the regions file.
To get the chromosomes in the regions::

	cut -f1 confident.bed | sort | uniq > chromosomes.txt

We can cat the contents of that file to xargs to have xargs put each chromosome
on the end of the command line. So to download and subset the reference::

	wget http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta.fai
	cat chromosomes.txt | \
	xargs samtools faidx http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta > ref.fa

Now we need to download the BAM file.
In my experience, downloading from the internet is much faster
if we specify the regions on the command line instead of using the
``-L`` flag, so we'll use ``awk`` to convert the regions in the BED
to the proper format for samtools and ``xargs`` to add the regions
to the end of the command line::

	cat confident.bed | \
	awk 'BEGIN{OFS=""}{print $1,":",($2+1),"-",$3}' > confident.regions

Since this is a huge amount of data,
it's better to download the index first and use samtools to just grab the
reads that intersect the confident regions we want. 
Even though we do this, the download will still take a minute or so::

	wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR134/ERR1341796/CHM1_CHM13_2.bam.bai
	cat confident.regions | \
	xargs samtools view -h -b -M -q 1 -F 3844 -o confident.bam ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR134/ERR1341796/CHM1_CHM13_2.bam

Once this download finishes you should have 4 important files:

#. confident.bed
#. confident.vcf.gz
#. confident.bam
#. ref.fa

.. note::

	You only need ``confident.bam`` if you plan to skip the benchmarking example.

Recalibrating a BAM file
-------------------------

.. warning::

	Recalibrating a BAM file is not yet supported.
	We hope to support this feature soon in an upcoming release.
	Once this feature is enabled, instructions to do so will be here.

.. _correcting_reads:

Correcting Reads
-----------------

The current implementation detects read errors by comparing the original read
with a corrected one output from an error correction method. In this tutorial
we use ``lighter`` [Song_2014]_, but any sufficiently-accurate method will work.

.. [Song_2014] Song, L., Florea, L. and Langmead, B., Lighter: Fast and Memory-efficient Sequencing Error Correction without Counting. Genome Biol. 2014 Nov 15;15(11):509.

.. note::

	We intend to implement a similar algorithm to make it easier to recalibrate
	data without having to correct and do other processing steps. Please look
	forward to it.

First, we need to extract the reads in the BAM. For more information
about this, check out :ref:`bam_to_fastq`. To do this,
sort the reads by name and use the ``samtools fastq`` command.
We will discard non-paired and singleton reads by setting the output
files for those reads to ``/dev/null``::

	samtools sort -@ 4 -n -o confident.nsorted.bam -O bam confident.bam
	samtools fastq -t -N -F 3844 -O -0 /dev/null -s /dev/null -1 reads.1.fq -2 reads.2.fq confident.nsorted.bam

We will then merge the fastq files to interleaved format.
We will also convert the whitespace to '_' characters using
the ``tr`` command::

	seqtk mergepe reads.1.fq reads.2.fq | tr ' ' _ > reads.fq

Now we correct the reads with ``lighter``. ``lighter`` requires
an output directory, so we'll need to create that before running
the command::

	mkdir lighter

Additionally, we need to calculate the best parameters
to use for the correction. Lighter requires 3 parameters: a k-mer
size, a genome size, and a sampling rate. I tend to use k=32 because
that seems to work well. If you don't know your approximate genome
size, take a guess. To calculate the sampling rate, the ``lighter``
authors recommend using ``7 / avg. read depth``. You can estimate
the depth by multiplying the number of reads by the read length
and dividing by the genome size.

You can also use ``samtools depth`` and a short ``awk`` script
to calculate these parameters::

	samtools depth -b confident.bed -m0 confident.bam | \
	awk '{x+=$3}END{print "bases:", NR, "\ndepth:", x/NR, "\nalpha:", 7/(x/NR)}'

Which should output something like::

	bases: 73465 
	depth: 49.1876 
	alpha: 0.142312

Now we can use lighter::

	lighter -r reads.fq -k 32 70000 .14 -od lighter

The output file we're interested in is ``lighter/reads.cor.fq``,
which we'll copy to our tutorial directory while removing the space
characters and replacing them with underscores as before::

	cat lighter/reads.cor.fq | tr ' ' _ > reads.cor.fq


Calibrating FASTQ Reads
-----------------------

The reference-free and alignment-free recalibration workflow works in three phases:

#. Find errors in reads
#. Build a model
#. Correct reads with the model

This is similar to GATK's method, but errors are detected by comparing corrected
and uncorrected reads rather than by comparing the aligned read to the reference.

Now that all the data has been properly prepared, to recalibrate::

	kbbq recalibrate --infer-rg -f reads.fq reads.cor.fq > reads.recalibrated.fq

The ``--infer-rg`` flag will ensure each read is assigned to its proper read group.
For more information on how it does this, check out :ref:`infer_rg`.

Benchmarking Calibration
-------------------------

In this case, the dataset has been validated and can serve as a truth set.
In a research environment, especially with nonmodel organisms, you won't
likely have a truth set you can use to check your calibration. :ref:`benchmark`
requires a set of variable sites in a VCF file and a set of confident regions
in a BED file. It assumes any site in the confident regions that are not variable
but don't match the reference are errors.

To check the calibration, we'll use the :ref:`benchmark` and :ref:`plot`
commands. Benchmark requires that the names in the FASTQ files match up with
the reads in the BAM file that comprises the truth set. Read
:ref:`read_matching` for more information.
Benchmark can also benchmark reads in a BAM file if a FASTQ file isn't provided.
To show both usages and to make a more interesting plot, we'll benchmark
the original qualities, the updated qualities
in the bam, and the qualities produced by KBBQ. We use the ``-l`` option to
specify the label that will go in the benchmark file and on the plot.

First, to benchmark the BAM reads::

	kbbq benchmark -l BAM -b confident.bam -r ref.fa -v confident.vcf.gz > bam.benchmark.txt

To benchmark the original qualities from the BAM::

	kbbq benchmark -l Original --use-oq -b confident.bam -r ref.fa -v confident.vcf.gz > oq.benchmark.txt

And to benchmark the qualities KBBQ assigned from the new fastq files::

	kbbq benchmark -l KBBQ -b confident.bam -r ref.fa -v confident.vcf.gz -f reads.recalibrated.fq > kbbq.benchmark.txt

Plotting the Benchmark
----------------------
:ref:`plot` is designed to make it as simple as possible to plot the output of :ref:`benchmark`.
Since we have multiple benchmark files, we can plot them individually by calling plot multiple
times or concatenate them into one file to plot all 3 lines on one plot. To do this::

	cat bam.benchmark.txt oq.benchmark.txt kbbq.benchmark.txt | \
	kbbq plot -o plot.pdf

The output name is passed directly to :func:`matplotlib.pyplot.savefig`, so the output
type is determined by the name.

The default type of plot made is ``calibration``, however, plot can also plot sample sizes,
which may be informative if you want to know how many bases of each quality score are in your
dataset. To plot these, use the ``-t`` option to change the plot type to ``sample-size``::

	cat bam.benchmark.txt oq.benchmark.txt kbbq.benchmark.txt | \
	kbbq plot -t sample-size -o sample-size.pdf

