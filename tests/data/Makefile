#Makefile to download and subset test data

#Variables
NUMREGIONS = 3 #Number of confident regions to use before sampling
NUMRGS = 2 #Only take sequences from the first NUMRGS rgs in the file
READSUBSAMPLING = 666.1 #The argument for samtools view -s. See the samtools manual.

#File names
FULLBEDNAME = CHM-eval.kit/full.37m.bed.gz
FULLVCFNAME = CHM-eval.kit/full.37m.vcf.gz
REFNAME = ref.fa.gz
REFIDX = ref.fa.gz.fai

CONFBEDNAME = conf_regions.bed
CONFVCFNAME = conf_regions.vcf.gz
CONFVCFIDX = $(addsuffix .csi, $(CONFVCFNAME))
CONFBAMNAME = conf_regions.bam
RGSFILE = rgs.txt

FIRSTREADS = reads.1.fq
SECONDREADS = reads.2.fq
FULLREADS = reads.fq
LIGHTERREADS = lighter.fq

RECALTABLE = conf_regions.recal.txt
RECALBAM = conf_regions.recal.bam
VARSITES = variable_sites.txt

#note this TAR cmd will fail unless the BED is made first.
define TARTOCMD =
'if [ $${TAR_FILENAME} == "$(FULLBEDNAME)" ]; then \
	zcat - | head -n $(NUMREGIONS) >$(CONFBEDNAME) ; \
elif [ $${TAR_FILENAME} == "$(FULLVCFNAME)" ]; then \
	bcftools view -T $(CONFBEDNAME) -Oz -o $(CONFVCFNAME) - ; \
fi'
endef

.PHONY: data

.SECONDARY:

data: $(RECALTABLE) $(RECALBAM) $(VARSITES)

%.vcf.gz.csi: %.vcf.gz
	bcftools index -c $<

%.fa.gz.fai: %.fa.gz
	samtools faidx $<

$(CONFBEDNAME) $(CONFVCFNAME):
	wget -qO- https://github.com/lh3/CHM-eval/releases/download/v0.5/CHM-evalkit-20180222.tar | \
	tar -xnf- --occurrence=1 --to-command=$(TARTOCMD) $(FULLBEDNAME) $(FULLVCFNAME)

$(REFNAME): $(CONFBEDNAME)
	samtools faidx http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta $(shell cut -f1 $< | sort | uniq) | \
	bgzip > $@

$(CONFBAMNAME): $(CONFBEDNAME)
	samtools view -h -L $< -u ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR134/ERR1341796/CHM1_CHM13_2.bam | \
	samtools view -h -u -F 3844
	samtools sort -@ 4 -O bam -o $@

$(VARSITES): $(CONFVCFNAME) $(CONFVCFIDX) $(CONFBEDNAME)
	bcftools query -R $(word 3, $^) -f '%CHROM\t%POS\n' $< >$@

$(FIRSTREADS) $(SECONDREADS): $(CONFBAMNAME)
	samtools sort -@ 4 -n -O bam $< | \
	samtools fastq -t -N -F 3844 -O -0 /dev/null -s /dev/null -1 $(FIRSTREADS) -2 $(SECONDREADS)

$(FULLREADS): $(FIRSTREADS) $(SECONDREADS)
	seqtk mergepe $< $(word 2, $^) | tr ' ' _ > $@

$(LIGHTERREADS): $(FULLREADS)
	LIGHTERTMP=`mktemp -d lighter.XXXXXX.d` &&\
	lighter -r $< -k 32 3000000 0.2 -t 16 -od $${LIGHTERTMP} &&\
	cat $${LIGHTERTMP}/$(basename $<).cor$(suffix $<) | tr ' ' _ > $@ &&\
	rm -rf $${LIGHTERTMP}

$(RECALTABLE): $(CONFBAMNAME) $(REFNAME) $(CONFVCFNAME)
	gatk BaseRecalibrator -I $< -R $(word 2, $^) --known-sites $(word 3, $^) --use-original-qualities -O $@

$(RECALBAM): $(CONFBAMNAME) $(REFNAME) $(RECALTABLE)
	gatk ApplyBQSR -I $< -R $(word 2, $^) --bqsr-recal-file $(word 3, $^) --use-original-qualities -O $@
