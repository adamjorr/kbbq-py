#Makefile to download and subset test data

#Variables
NUMREGIONS = 3 #Number of confident regions to use before sampling
NUMRGS = 2 #Only take sequences from the first NUMRGS rgs in the file
READSUBSAMPLING = -s 666.1 #The argument for samtools view -s. See the samtools manual.
BAMURL = ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR134/ERR1341796/CHM1_CHM13_2.bam
#BAMACCESSION = ERR1341796
TARURL = https://github.com/lh3/CHM-eval/releases/download/v0.5/CHM-evalkit-20180222.tar
REFURL = http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta
FULLREFIDXURL = $(REFURL).fai
export SHELL = /bin/bash

#File names
FULLBEDNAME = CHM-eval.kit/full.37m.bed.gz
FULLVCFNAME = CHM-eval.kit/full.37m.vcf.gz
REFNAME = ref.fa
REFIDX = $(REFNAME).fai
REFDICT = $(basename $(REFNAME)).dict

CONFBEDNAME = conf_regions.bed
CONFVCFNAME = conf_regions.vcf.gz
CONFVCFIDX = $(addsuffix .tbi, $(CONFVCFNAME))
FULLBAMNAME = CHM1_CHM13_2.bam
CONFBAMNAME = conf_regions.bam
RGSFILE = rgs.txt

FIRSTREADS = reads.1.fq
SECONDREADS = reads.2.fq
FULLREADS = reads.fq
LIGHTERREADS = lighter.fq

RECALTABLE = conf_regions.recal.txt
RECALBAM = conf_regions.recal.bam
VARSITES = variable_sites.txt

BEDTOREGIONS='BEGIN{OFS=""}{print $$1,":",($$2+1),"-",$$3}'


.PHONY: data clean


data: $(RECALTABLE) $(RECALBAM) $(VARSITES) $(FULLREADS)

clean:
	rm -f $(RECALTABLE) $(RECALBAM) $(VARSITES) $(FULLBAMNAME).bai \
	$(CONFBEDNAME) $(REFNAME) $(REFDICT) $(REFIDX) \
	$(CONFVCFNAME) $(CONFVCFIDX) $(CONFBAMNAME) $(RGSFILE) \
	$(FIRSTREADS) $(SECONDREADS) $(FULLREADS) $(LIGHTERREADS) \
	$(basename $(RECALBAM)).bai $(notdir $(FULLREFIDXURL))

# Generic Indexing Rules
########################

%.vcf.gz.csi: %.vcf.gz
	bcftools index -cf $<

%.vcf.gz.tbi: %.vcf.gz
	bcftools index -tf $<

%.fa.fai: %.fa
	samtools faidx $<

%.dict: %.fa
	rm -rf $@ && picard CreateSequenceDictionary R=$<

#

$(CONFBEDNAME):
	wget -qO- $(TARURL) | \
	tar -xnf- --occurrence=1 --to-stdout $(FULLBEDNAME) | \
	zcat - | head -n $(NUMREGIONS) > $(CONFBEDNAME)

$(CONFVCFNAME): $(CONFBEDNAME)
	wget -qO- $(TARURL) | \
	tar -xnf- --occurrence=1 --to-stdout $(FULLVCFNAME) | \
	bcftools view -T $(CONFBEDNAME) -Oz -o $(CONFVCFNAME) -

$(notdir $(FULLREFIDXURL)):
	wget -q $(FULLREFIDXURL) -O $@

$(REFNAME): $(CONFBEDNAME) $(notdir $(FULLREFIDXURL))
	samtools faidx $(REFURL) $(shell cut -f1 $< | sort | uniq) > $@

$(RGSFILE): | $(FULLBAMNAME).bai
	samtools view -H $(BAMURL) | \
	grep ^@RG | cut -f 2 | cut -d: -f2 | \
	head -n $(NUMRGS) > $@

$(FULLBAMNAME).bai:
	lftp -c "pget -n 10 $(BAMURL).bai -o $@"

$(CONFBAMNAME): $(CONFBEDNAME) $(RGSFILE) | $(FULLBAMNAME).bai
	samtools view -h -b -M -q 1 -F 3844 -R $(word 2, $^) $(READSUBSAMPLING) -o $@ $(BAMURL) $(shell cat $< | awk $(BEDTOREGIONS))

$(VARSITES): $(CONFVCFNAME) $(CONFVCFIDX) $(CONFBEDNAME)
	bcftools query -R $(word 3, $^) -f '%CHROM\t%POS\n' $< >$@

$(FIRSTREADS) $(SECONDREADS): $(CONFBAMNAME)
	samtools sort -@ 4 -n -O bam $< | \
	samtools fastq -t -N -F 3844 -O -0 /dev/null -s /dev/null -1 $(FIRSTREADS) -2 $(SECONDREADS) -

$(FULLREADS): $(FIRSTREADS) $(SECONDREADS)
	seqtk mergepe $< $(word 2, $^) | tr ' ' _ > $@

$(LIGHTERREADS): $(FULLREADS)
	LIGHTERTMP=`mktemp -d lighter.XXXXXX.d` &&\
	lighter -r $< -k 32 3000000 0.2 -t 16 -od $${LIGHTERTMP} &&\
	cat $${LIGHTERTMP}/$(basename $<).cor$(suffix $<) | tr ' ' _ > $@ &&\
	rm -rf $${LIGHTERTMP}

$(RECALTABLE): $(CONFBAMNAME) $(REFNAME) $(CONFVCFNAME) $(CONFVCFIDX) $(REFIDX) $(REFDICT)
	gatk BaseRecalibrator -I $< -R $(word 2, $^) --known-sites $(word 3, $^) --use-original-qualities -O $@

$(RECALBAM): $(CONFBAMNAME) $(REFNAME) $(RECALTABLE) $(REFIDX) $(REFDICT)
	gatk ApplyBQSR -I $< -R $(word 2, $^) --bqsr-recal-file $(word 3, $^) --use-original-qualities -O $@
