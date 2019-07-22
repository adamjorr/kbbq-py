FASTQ Input Help
================

Currently we only support FASTQ input reads for the :ref:`recalibrate` command.
We also support only one FASTQ file (and its corrected variant) at a time.
Because of this, there are some onerous requirements for reads names and length.

We plan to support input reads from any number of BAM or FASTQ files before release.
We also plan to remove the requirement of externally correcting the reads before input.

FASTQ Input Requirements
------------------------

+----------------+-------------------------------------------------------------+
| Feature        | Requirement                                                 |
+================+=============================================================+
| ``--infer-rg`` | Add ``RG:`` to ID and append to name delimited with ``_``   |
+----------------+-------------------------------------------------------------+
| Recalibrate    | Corrected read name begins with uncorrected read name       |
+----------------+-------------------------------------------------------------+
| Paired reads   | First ``_`` delimited field of name ends with ``/2`` for    |
|                | second in pair reads                                        |
+----------------+-------------------------------------------------------------+
| Benchmark      | First ``_`` delimited field of name (without ``/1`` or      |
|                | ``/2`` suffix) must match a read name in the BAM            |
+----------------+-------------------------------------------------------------+

Only one FASTQ file of reads is currently supported. You must supply the
file and a corrected version of the file on the command line. For example,
if your reads are in ``reads.fq`` and a corrected version is in ``reads.cor.fq``,
you would use the :ref:`recalibrate` command like::

	kbbq recalibrate -f reads.fq reads.cor.fq

To ensure the reads in the uncorrected file and corrected file are properly
aligned, :ref:`recalibrate` will check that the corrected version of each read
has a name that begins with the name of the uncorrected read. This requirement
is satisfied when ``lighter`` is used as the error corrector, but should be
satisfied by other error correctors as well.

.. _infer_rg:

Inferring Read Groups
*********************

Currently the only other command line option to :ref:`recalibrate` that applies
to FASTQ files is ``--infer-rg``, which will attempt to parse the read name for
the read group the read belongs to. If you don't want to mess with renaming your
reads or already have them in sets of 1 read group per file, you can safely use
:ref:`recalibrate` on each file separately; the default behavior without the
``--infer-rg`` flag is to treat all reads in the file as belonging to the same
read group.

To properly infer the read group, the read must have the read group appended
with a ``_`` character to the end of the read name. The read group should also
include a ``RG:`` prefix, and may include other information (such as a type
indicator) delimited by a ``:`` before the read group ID.
Thus a read name like ``@HJCMTCCXX160113:5:1101:7760:55965/1`` in read group
``HJCMT.5`` is made suitable for parsing by adding ``_RG:HJCMT.5`` to the end,
forming the full name ``@HJCMTCCXX160113:5:1101:7760:55965/1_RG:HJCMT.5``. It
is also equally valid to include additional information before the actual ID.
For example, it's OK to add the ``Z`` type indicator from samtools indicating
the RG is a string type, and in that case the full read name would be::

	@HJCMTCCXX160113:5:1101:7760:55965/1_RG:Z:HJCMT.5

First or Second in Pair
***********************

The program will interpret any read where the last 2 characters of the first field
of a ``_`` delimited name are ``/2`` as a second-in-pair read. Any other name and the
read will be interpreted as first in pair. If your reads are paired and match this naming
scheme or if your reads are unpaired, there is nothing you must do to mark your reads.

All the reads below will be interpreted as 2nd in pair:

	- ``@HJCMTCCXX160113:5:1101:7760:55965/2_RG:Z:HJCMT.5``
	- ``@HJCMTCCXX160113:5:1101:7760:55965/2``
	- ``@HJCMTCCXX160113:5:1101:7760:55965/2_foo``

All the reads below will be interpreted as **not** 2nd in pair:

	- ``@HJCMTCCXX160113:5:1101:7760:55965/1_RG:Z:HJCMT.5``
	- ``@HJCMTCCXX160113:5:1101:7760:55965/1``
	- ``@HJCMTCCXX160113:5:1101:7760:55965``
	- ``@HJCMTCCXX160113:5:1101:7760:55965/3``
	- ``@HJCMTCCXX160113:5:1101:7760:55965_2``
	- ``@HJCMTCCXX160113:5:1101:7760:55965_foo_/2``

If your reads **are** paired, but the second-in-pair reads are not properly marked,
recalibration *may* be less effective, though I haven't seen data to indicate that.

.. _read_matching:

Benchmark Read Matching
***********************

For :ref:`benchmark` to properly match the read in a FASTQ with the read
in the BAM file, the first ``_`` delimited field must match the read name,
minus any ``/1`` or ``/2`` parts of the read name. For example, to match
the 2nd in pair read ``HK2WYCCXX160124:1:1219:24545:4315``, the FASTQ represented
by FASTQ lines,

	- ``@HJCMTCCXX160113:5:1101:7760:55965/2_RG:Z:HJCMT.5``
	- ``@HJCMTCCXX160113:5:1101:7760:55965/2``
	- ``@HJCMTCCXX160113:5:1101:7760:55965/2_foo``

will all work. However, to match a READ2-flagged read, :ref:`benchmark`
must be able to determine that the read is 2nd in pair. Read the discussion
above for how the program determines this. Reads flagged as READ1 or without
either READ1 or READ2 flags set in the bam can both be matched with FASTQ
reads that aren't interpreted as 2nd in pair.

.. _bam_to_fastq:

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

If you used the ``-t`` option to add RG tags, you'll want to remove the spaces inserted by ``samtools``
as many error correctors won't support them. Currently ``kbbq`` enforces a ``_`` character delimiter,
but this requirement will be eased in future releases.
The ``tr`` command can efficiently replace the spaces with ``_`` like this::

	cat reads.merged.fq | tr ' ' _ > reads.fixed.fq
