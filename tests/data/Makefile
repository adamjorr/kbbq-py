#Makefile to download and subset test data

#Variables
NUMREGIONS = 3 #Number of confident regions to use before sampling
NUMRGS = 2 #Only take sequences from the first NUMRGS rgs in the file
READSUBSAMPLING = -s 666.1 #The argument for samtools view -s. See the samtools manual.
BAMURL = ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR134/ERR1341796/CHM1_CHM13_2.bam
#BAMACCESSION = ERR1341796
TARURL = https://github.com/lh3/CHM-eval/releases/download/v0.5/CHM-evalkit-20180222.tar
REFURL = http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta

#File names
FULLBEDNAME = CHM-eval.kit/full.37m.bed.gz
FULLVCFNAME = CHM-eval.kit/full.37m.vcf.gz
REFNAME = ref.fa
REFIDX = ref.fa.fai

CONFBEDNAME = conf_regions.bed
CONFVCFNAME = conf_regions.vcf.gz
CONFVCFIDX = $(addsuffix .csi, $(CONFVCFNAME))
FULLBAMNAME = CHM1_CHM13_2.bam
SUBSETBAMNAME = unsorted.conf_regions.bam
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

BEDTOREGIONS='BEGIN{OFS=""; ORS=" "}{print $$1,":",($$2+1),"-",$$3}'


.PHONY: data

.SECONDARY:

data: $(RECALTABLE) $(RECALBAM) $(VARSITES)

%.vcf.gz.csi: %.vcf.gz
	bcftools index -c $<

%.fa.fai: %.fa
	samtools faidx $<

$(CONFBEDNAME) $(CONFVCFNAME):
	wget -qO- $(TARURL) | \
	tar -xnf- --occurrence=1 --to-command=$(TARTOCMD) $(FULLBEDNAME) $(FULLVCFNAME)

$(REFNAME): $(CONFBEDNAME)
	samtools faidx $(REFURL) $(shell cut -f1 $< | sort | uniq) > $@

$(RGSFILE): $(FULLBAMNAME)
	samtools view -H $< | \
	grep ^@RG | cut -f 2 | cut -d: -f2 | \
	head -n $(NUMRGS) > $@

#$(FULLBAMNAME): $(CONFBEDNAME)
#	sam-dump -r -1 $(addprefix --aligned-region=,$(shell cut -f1 $< | sort | uniq)) $(BAMACCESSION) | \
#	samtools view -b -o $@

#$(FULLBAMNAME):
#	lftp -c "pget -n 10 $(BAMURL) -o $@"

$(FULLBAMNAME).bai:
	lftp -c "pget -n 10 $(BAMURL).bai -o $@"

$(SUBSETBAMNAME): $(FULLBAMNAME).bai $(CONFBEDNAME) $(RGSFILE)
	samtools view -h -b -M -o $@ $(BAMURL) $(shell cat $(word 2, $^) | awk $(BEDTOREGIONS))

#options cut to just see download speed:
#-F 3844 $(READSUBSAMPLING) -R $(word 3, $^)

#$(SUBSETBAMNAME): $(FULLBAMNAME) $(CONFBEDNAME) $(RGSFILE) $(FULLBAMNAME).bai
#	samtools view -h -L $(word 2, $^) -F 3844 -R $(word 3, $^) $(READSUBSAMPLING) -b -o $@ $<

$(CONFBAMNAME): $(SUBSETBAMNAME)
	samtools sort -@ 4 -O bam -o $@ $<

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

$(RECALTABLE): $(CONFBAMNAME) $(REFNAME) $(CONFVCFNAME) $(REFIDX)
	gatk BaseRecalibrator -I $< -R $(word 2, $^) --known-sites $(word 3, $^) --use-original-qualities -O $@

$(RECALBAM): $(CONFBAMNAME) $(REFNAME) $(RECALTABLE) $(REFIDX)
	gatk ApplyBQSR -I $< -R $(word 2, $^) --bqsr-recal-file $(word 3, $^) --use-original-qualities -O $@
