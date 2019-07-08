Command Line Help
*****************

::


    usage: kbbq [-h] [-v] {recalibrate,benchmark,plot} ...
    
    K-mer Based Base Quality score recalibration
    
    optional arguments:
      -h, --help            show this help message and exit
      -v, --version         show program's version number and exit
    
    command:
      valid commands
    
      {recalibrate,benchmark,plot}

.. _recalibrate:

recalibrate
-----------

.. highlight:: console

.. warning::

   Recalibrating reads from a BAM file is not yet supported.
   Use the :code:`samtools fastq` command to convert your reads
   to fastq format. This should go something like::

      samtools sort -n -O bam input.bam > namesorted.bam
      samtools fastq -t -N -F 3844 -O -0 /dev/null -s /dev/null -1 reads.1.fq -2 reads.2.fq namesorted.bam

   Reads must then be interleaved for input.
   This can be done with a tool like seqtk::

      seqtk mergepe reads.1.fq reads.2.fq > reads.merged.fq

   See :doc:`fastq_input` for more information.

::


    usage: kbbq recalibrate [-h] (-b BAM | -f FASTQ FASTQ) [-u] [-s]
                            [-g GATKREPORT] [--infer-rg]
    
    Recalibrate a BAM or FASTQ file
    
    optional arguments:
      -h, --help            show this help message and exit
      -b BAM, --bam BAM     BAM to recalibrate
      -f FASTQ FASTQ, --fastq FASTQ FASTQ
                            FASTQ file to recalibrate and a corrected version from
                            your favorite error corrector.
      -u, --use-oq          Use the OQ tag to get quality scores when working with
                            a BAM file. Does nothing if a fastq file is provided.
      -s, --set-oq          Set the 'OQ' flag prior to recalibration. Only works
                            when producing BAM output.
      -g GATKREPORT, --gatkreport GATKREPORT
                            If the given path points to an existing GATK report,
                            load the model from the report instead of calculating
                            it. If the file doesn't exist, save the calculated
                            model to the given path.
      --infer-rg            Attempt to infer the read group from a FASTQ read.
                            Only works with FASTQ input. The default behavior is
                            to treat each input FASTQ file as its own read group.

.. _benchmark:

benchmark
---------

::


    usage: kbbq benchmark [-h] -b BAM -r REFERENCE -v VCF [-f FASTQ] [-l LABEL]
                          [-u] [-d BEDFILE]
    
    Benchmark a BAM or FASTQ file using a truth set
    
    optional arguments:
      -h, --help            show this help message and exit
      -f FASTQ, --fastq FASTQ
                            fastq file to benchmark
      -l LABEL, --label LABEL
                            label to use for label column
      -u, --use-oq          Use the OQ tag to get quality scores when working with
                            a BAM file. Does nothing if a fastq file is provided.
      -d BEDFILE, --bedfile BEDFILE
                            BED file of confident regions. Sites outside the given
                            regions will be skipped.
    
    required arguments:
      -b BAM, --bam BAM     Truth set BAM file. Differences from the reference at
                            nonvariable sites will be interpreted as errors.
      -r REFERENCE, --reference REFERENCE
                            FASTA file containing the reference genome
      -v VCF, --vcf VCF     VCF file containing variable sites

.. _plot:

plot
----

::


    usage: kbbq plot [-h] [-t {calibration,sample-size}] -o OUTFILE [file]
    
    Plot data output from the benchmark command
    
    positional arguments:
      file                  Input file
    
    optional arguments:
      -h, --help            show this help message and exit
      -t {calibration,sample-size}, --type {calibration,sample-size}
                            Type of plot to produce
    
    required arguments:
      -o OUTFILE, --outfile OUTFILE
                            file name to save plot as



