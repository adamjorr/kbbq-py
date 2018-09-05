#!/usr/bin/env python3
import pysam
import numpy as np
from sklearn.linear_model import LogisticRegression as LR
from sklearn.isotonic import IsotonicRegression as IR
import importlib.util
spec = importlib.util.spec_from_file_location("ek", "/home/ajorr1/bin/jelly/error_kmers.py")
ek = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ek)
import os.path
import sys
import re

def new_count_qual_scores(samfile, ref, vcf, calibrated_quals):
    print(ek.tstamp(), "Counting Base Quality Scores . . .", file=sys.stderr)
    numerrors = np.zeros(43, dtype = np.uint64)
    numtotal = np.zeros(43, dtype = np.uint64)
    caliberrs = np.zeros(43, dtype = np.uint64)
    calibtotal = np.zeros(43, dtype = np.uint64)
    for read in samfile:
        if (read.is_unmapped or read.mate_is_unmapped):
            continue
        if (read.is_secondary or read.is_supplementary):
            continue

        refchr = read.reference_name
        readname = read.query_name
        assert read.is_read1 or read.is_read2
        suffix = ("1" if read.is_read1 else "2")
        fullname = readname + "/" + suffix
        if not fullname in calibrated_quals.keys():
            continue

        refpositions = read.get_reference_positions(full_length=True) #these should be 1-based positions but are actually 0-based
        errorpositions = np.array([i for i,pos in enumerate(refpositions) if pos is None or (read.query_sequence[i] != ref[refchr][pos] and pos+1 not in vcf[refchr])], dtype=np.intp)
        quals = np.array(read.query_qualities, dtype=np.intp)
        cquals = calibrated_quals[fullname]
        np.add.at(numtotal,quals,1)
        np.add.at(numerrors,quals[errorpositions],1)
        np.add.at(calibtotal, cquals,1)
        np.add.at(caliberrs, cquals[errorpositions],1)
    return numerrors, numtotal, caliberrs, calibtotal

def count_quals_from_plp(plpfilename, var_pos, calibrated_quals, uncalibrated_quals, suffix):
    numerrors = np.zeros(43, dtype = np.uint64)
    numtotal = np.zeros(43, dtype = np.uint64)
    caliberrs = np.zeros(43, dtype = np.uint64)
    calibtotal = np.zeros(43, dtype = np.uint64)
    numfinder = re.compile('[\+-](\d+)')
    with open(plpfilename, 'r') as infh:
        varidx = {k : 0 for k in var_pos.keys()}
        for line in infh:
            chrom, pos, refbase, depth, bases, quals, qpos, qname = line.split('\t')
            pos = int(pos)
            if var_pos[chrom][varidx[chrom]] == pos:
                #increment unless it's the last position in the chromosome
                varidx[chrom] = (varidx[chrom] + 1 if varidx[chrom] < len(var_pos[chrom]) - 1 else varidx[chrom])
                continue
            if int(depth) == 0:
                continue
            assert refbase != 'N'
            assert var_pos[chrom][varidx[chrom]] > pos
            bases = re.sub('\^.', '', bases)
            bases = bases.replace('$','')
            bases = bases.replace('<','')
            bases = bases.replace('>','')
            
            match = numfinder.search(bases)
            while match:
                try:
                    bases = bases[:match.start()] + bases[(match.end() + int(match[1])):]
                except TypeError:
                    print("Bases:",bases)
                    print("Match:",match)
                    raise
                match = numfinder.search(bases)
            
            #bases = re.sub('\+[0-9]+[ACGTNacgtn]+', '', bases)
            #bases = re.sub('-[0-9]+[ACGTNacgtn]+', '', bases)
            try:
                assert len(bases) == int(depth)
            except AssertionError:
                print("Bases:", bases)
                print("Depth:", depth)
                print("File:", plpfilename)
                print("Line:", line)
                raise
                
            assert len(bases) == len(quals)
            bases = np.array(list(bases), dtype = np.unicode_)
            erroneous = np.array(np.logical_and(bases != '.', bases != ','))

            qpos = np.array(qpos.split(','), dtype = np.int) - 1
            qname = qname.rstrip()
            qname = np.array(qname.split(','), dtype = np.unicode_)
            qname = np.core.defchararray.add(qname, suffix)
            try:
                quals = np.array([uncalibrated_quals[qname[i]][qpos[i]] for i in range(len(qname))], dtype = np.int)
                calibquals = np.array([calibrated_quals[qname[i]][qpos[i]] for i in range(len(qname))], dtype = np.int)
            except KeyError:
                continue
            np.add.at(numerrors, quals[erroneous], 1)
            np.add.at(numtotal, quals, 1)
            np.add.at(caliberrs, calibquals[erroneous], 1)
            np.add.at(calibtotal, calibquals, 1)
    return numerrors, numtotal, caliberrs, calibtotal

def load_positions(posfile):
    d = dict()
    with open(posfile, 'r') as infh:
        for line in infh:
            chrom, pos = line.rstrip().split()
            d.setdefault(chrom, list()).append(int(pos))
    return d

def find_rcorrected_sites(uncorrfile, corrfile):
    print(ek.tstamp(), "Finding Rcorrected sites", file=sys.stderr)
    uncorr_set = list(pysam.FastxFile(uncorrfile))
    corr_set = list(pysam.FastxFile(corrfile))
    #verify the sequences are the same and can be accessed by index
    for i in range(len(corr_set)):
        try:
            assert corr_set[i].name.startswith(uncorr_set[i].name)
        except AssertionError:
            print("Corr_set[i]:",corr_set[i])
            print("Uncorr_Set[i]:", uncorr_set[i])
            raise

    

def main():
    np.seterr(all = 'raise')
    print(ek.tstamp(), "Starting . . .", file=sys.stderr)
    uncorrfile = "reads.fq"
    corrfile = "nospace.reads.cor.fq"

    uncorr_set = list(pysam.FastxFile(uncorrfile))
    corr_set = list(pysam.FastxFile(corrfile))
    #verify the sequences are the same and can be accessed by index
    for i in range(len(corr_set)):
        try:
            assert corr_set[i].name.startswith(uncorr_set[i].name)
        except AssertionError:
            print("Corr_set[i]:",corr_set[i])
            print("Uncorr_Set[i]:", uncorr_set[i])
            raise
    print(ek.tstamp(),"1")
    uncorr_reads, corr_reads = uncorr_set, corr_set
    
    names = np.array([r.name for r in uncorr_reads], dtype = np.unicode_)
    uncorr_seqs = np.array([str(r.sequence) for r in uncorr_reads], dtype = np.unicode_)
    uncorr_quals = np.array([r.get_quality_array() for r in uncorr_reads], dtype = np.int)
    corr_seqs = np.array([str(r.sequence) for r in corr_reads], dtype = np.unicode_)

    print(ek.tstamp(),"2")

    nonerroneous_scores = np.zeros(43, dtype = np.int)
    erroneous_scores = np.zeros(43, dtype = np.int)
    for i in range(len(corr_seqs)):
        uncorr_s = uncorr_seqs[i]
        corr_s = corr_seqs[i]
        sites_corrected = np.array([uncorr_s[j] != corr_s[j] for j in range(len(uncorr_s))])
        np.add.at(nonerroneous_scores, uncorr_quals[i][np.logical_not(sites_corrected)].flatten(), 1)
        np.add.at(erroneous_scores, uncorr_quals[i][sites_corrected].flatten(), 1)

    print(ek.tstamp(), "3")

    #build "training" data for logit regression
    x_nonerror = np.repeat(np.arange(43), nonerroneous_scores)
    x_error = np.repeat(np.arange(43), erroneous_scores)
    x = np.concatenate([x_nonerror, x_error])
    y = np.repeat([0,1], [np.sum(nonerroneous_scores),np.sum(erroneous_scores)])
    assert len(x) == len(y)
    lr = LR(tol = 1e-8)
    lr.fit( x.reshape(-1,1), y)

    print(ek.tstamp(), "4")

    #ir = IR( out_of_bounds = 'clip' )
    #ir.fit( x, y)
    # how to get recalibrated probability: 
    #newp = lr.predict_proba( p_test.reshape(-1,1))[:,1]

    recalibratedfile = "reads.recalibrated.fq"
    qtostr = np.arange(43, dtype = np.uint32) + 33
    qtostr = qtostr.view('U1')
    uncalibrated_quals = dict()
    calibrated_quals = dict()
    with open(recalibratedfile, 'w') as fout:
        for entry in uncorr_set:
            oldp = np.array(entry.get_quality_array(), dtype = np.int)
            newp = np.array(lr.predict_proba(oldp.reshape(-1,1))[:,1], dtype = np.longdouble)
            #newp = np.array(ir.transform(oldp), dtype = np.longdouble)
            q = -10.0*np.log10(newp)
            quals = np.array(np.rint(q), dtype=np.int)
            quals = np.clip(quals, 0, 43)
            entry.quality = ''.join(qtostr[quals])
            fout.write(str(entry))
            uncalibrated_quals[entry.name] = np.array(oldp)
            calibrated_quals[entry.name] = np.array(quals, dtype = np.int)

    print(ek.tstamp(), "5")

    bad_positions = load_positions("variable_sites.txt")
    numerrs1, numtotal1, caliberrs1, calibtotal1 = count_quals_from_plp('only_confident.1.plp', bad_positions, calibrated_quals, uncalibrated_quals, "/1")
    numerrs2, numtotal2, caliberrs2, calibtotal2 = count_quals_from_plp('only_confident.2.plp', bad_positions, calibrated_quals, uncalibrated_quals, "/2")
    print(ek.tstamp(), "6")
    numerrs = numerrs1 + numerrs2
    numtotal = numtotal1 + numtotal2
    caliberrs = caliberrs1 + caliberrs2
    calibtotal = calibtotal1 + calibtotal2
    print(ek.tstamp(), "\nNumerrs:", numerrs, "\nNumtotal", numtotal, "\nCaliberrs", caliberrs, "\nCalibtotal", calibtotal, file=sys.stderr)
    ek.plot_qual_scores(numerrs, numtotal, "qualscores.png", "Raw Reads")
    ek.plot_qual_scores(caliberrs, calibtotal, "calibrated.png", "After Calibration")
    assert np.all(numerrs <= numtotal)
    assert np.all(caliberrs <= calibtotal)

if __name__ == '__main__':
    main()

