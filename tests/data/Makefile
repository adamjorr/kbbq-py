#Makefile to download and subset test data

#Variables
NUMREGIONS = 10 #Number of confident regions to use before sampling
FULLBEDNAME = CHM-eval.kit/full.37m.bed.gz
FULLVCFNAME = CHM-eval.kit/full.37m.vcf.gz

CONFBEDNAME = conf_regions.bed
CONFVCFNAME = conf_regions.vcf.gz
CONFBAMNAME = conf_regions.bam

FIRSTREADS = reads.1.fq
SECONDREADS = reads.2.fq
FULLREADS = reads.fq
LIGHTERREADS = lighter.fq

#note this TAR cmd will fail unless the BED is made first.
define TARTOCMD =
if [ ${TAR_FILENAME} == "$(FULLBEDNAME)" ]; then
	zcat - | head -n 10 >$(CONFBEDNAME)
elif [ ${TAR_FILENAME} == "$(FULLVCFNAME)" ]; then
	bcftools view -T $(CONFBEDNAME) -Oz -o $(CONFVCFNAME) -
fi
endef

.PHONY: data

data: $(CONFBEDNAME) $(CONFVCFNAME)

$(CONFBEDNAME) $(CONFVCFNAME):
	wget -qO- https://github.com/lh3/CHM-eval/releases/download/v0.5/CHM-evalkit-20180222.tar | \
	tar -xnf- --occurrence=1 --to-command='$(TARTOCMD)' $(FULLBEDNAME) $(FULLVCFNAME)

$(CONFBAMNAME): $(CONFBEDNAME)
	samtools view -h -L $< -F 3844 ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR134/ERR1341796/CHM1_CHM13_2.bam | \
	samtools sort -O bam -o $@

$(FIRSTREADS) $(SECONDREADS): $(CONFBAMNAME)
	samtools sort -n -O bam $< | \
	samtools fastq -t -N -F 3844 -O -0 /dev/null -s /dev/null -1 $(FIRSTREADS) -2 $(SECONDREADS)

$(FULLREADS): $(FIRSTREADS) $(SECONDREADS)
	seqtk mergepe $< $(word 2, $^) | tr ' ' _ > $@

LIGHTERTMP := $(shell mktemp -d lighter_XXXXXX)
$(LIGHTERREADS): $(FULLREADS)
	lighter -r $< -k 32 3000000 0.2 -t 16 -od $(LIGHTERTMP) &&\
	cat $(LIGHTERTMP)/$(basename $<).cor$(suffix $<) | tr ' ' _ > $@ &&\
	rm -rf $(LIGHTERTMP)

#TODO:
#	"only_confident.sorted.recal.bam"
#	"../chr1.renamed.fa"
#	"variable_sites.txt"
