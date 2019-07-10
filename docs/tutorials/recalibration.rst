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

To download and install samtools::

	wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
	bunzip samtools-1.9.tar.bz2
	./configure --prefix=~/.local
	make
	make install

To install bcftools, follow the same instructions but replace all occurrences
of ``samtools`` with ``bcftools``

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

Finally, you'll need to have ``kbbq`` installed. Follow the instructions in :doc:`installation` for more information.

.. note::

	These instructions install to the ``.local/bin`` directory in your home directory.
	This follows the convention ``pip install --user`` uses, and is probably already on your PATH.
	If it isn't, edit ``~/.bash_profile`` and add the following line::

		export PATH=${PATH}:~/.local/bin

	Then run::

		. ~/.bash_profile

	to reload the variable.

Downloading and Preparing Data
-------------------------------

Recalibrating a BAM file
-------------------------

.. warning::

	Recalibrating a BAM file is not yet supported. We hope to support this feature soon in an upcoming release.
	Once this feature is enabled, instructions to do so will be here.

Correcting Reads
-----------------

Calibrating FASTQ Reads
-----------------------

Benchmarking Calibration
-------------------------

In this case, the dataset has been validated and can serve as a truth set.
In a research environment, especially with nonmodel organisms, you won't
likely have a truth set you can use to check your calibration. The :ref:`benchmark`
requires a set of variable sites in a VCF file and a set of confident regions
in a BED file. It assumes any site in the confident regions that are not variable
but don't match the reference are errors.

To check the calibration, we'll use the :ref:`benchmark` and :ref:`plot`
commands. Benchmark requires that the names in the FASTQ files match up with
the reads in the BAM file that comprises the truth set. It will also work
with just the reads in the BAM file. 
