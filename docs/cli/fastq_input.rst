FASTQ Input Help
================

Currently we only support FASTQ input reads for the ``recalibrate`` command.
We also support only one FASTQ file (and its corrected variant) at a time.
Because of this, there are some onerous requirements for reads names and length.

We plan to support input reads from any number of BAM or FASTQ files before release.
We also plan to remove the requirement of externally correcting the reads before input.

FASTQ Input Requirements
------------------------
TODO

Converting from BAM to FASTQ
----------------------------

.. highlight:: console

Recalibrating reads from a BAM file is not yet supported.
Use the :code:`samtools fastq` command to convert your reads
to fastq format. This should go something like::

  samtools sort -n -O bam input.bam > namesorted.bam
  samtools fastq -t -N -F 3844 -O -0 /dev/null -s /dev/null -1 reads.1.fq -2 reads.2.fq namesorted.bam

This will sort your input by read name (required for :code:`samtools fastq`),
:code:`-t` will append RG tags to the output read name (currently required for recalibration),
:code:`-N` add ``/1`` or ``/2`` to the output read name (also currentlyrequired for recalibration),
:code:`-F` will exclude unmapped, not primary alignment, QC failed, optical duplicate, and supplementary alignment reads,
and :code:`-O` will use OQ tags to obtain the read quality if available.
Reads without a READ1 or READ2 flag and singletons will not be output, while READ1 reads will be output to ``reads.1.fq``
and READ2 reads will be output to ``reads.2.fq``.

Reads must then be interleaved for input.
This can be done with a tool like seqtk::

  seqtk mergepe reads.1.fq reads.2.fq > reads.merged.fq
