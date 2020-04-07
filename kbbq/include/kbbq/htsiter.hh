#include <iterator>
#include <string>
#include <random>
#include "include/htslib/htslib/hts.h"
#include "minionrng/minion.hpp"
#include "readutils.hh"
#include <zlib.h>
KSEQ_INIT(gzFile, gzread);


//unsigned char* s = bam_get_seq(bamrecord);
namespace htsiter{

class HTSFile{
public:
	HTSFile();
	~HTSFile();
	virtual int next()=0;
	virtual std::string next_str()=0;
	virtual readutils::CReadData get()=0;
	virtual void recalibrate(std::vector<int> qual)=0;
	virtual int open_out(std::string filename)=0; //open an output file so it can be written to later.
	virtual int write()=0; //write the current read to the opened file.
}

class BamFile: public HTSFile{
public:
	samFile *sf;
	hts_idx_t *idx;
	sam_hdr_t *h;
	hts_itr_t *itr;
	bam1_t *r;
	samFile *of;
	BamFile(std::string filename){
		sf = sam_open(filename.c_str(), "r");
	   	idx = sam_index_load(sf, filename.c_str());
	    h = sam_hdr_read(sf);
	    itr = sam_itr_queryi(idx, HTS_IDX_START, 0, 0); //iterate over whole file;
	}
	~BamFile(){
		bam_destroy1(record);
		sam_itr_destroy(itr);
		sam_hdr_destroy(h);
		hts_idx_destroy(idx);
		sam_close(sf);
		if(outfile != NULL){sam_close(outfile);}
	}

	// to use: while (ret = BamFile.next() >= 0){//do something with this->r}
	int next(){return sam_itr_next(sf, itr, r);}
	// return next read as a string. if there are no more, return the empty string.
	std::string next_str(){return this->next() >= 0 ? std::string(bam_get_seq(r)) : ""}
	//
	readutils::CReadData get(){return readutils::CReadData(this->r);}
	//
	void recalibrate(std::vector<int> qual){char* q = (char *)bam_get_qual(this->r); for(int i = 0; i < this->r->core.l_qseq; ++i){q[i] = (char)qual[i]};}
	// TODO:: add a PG tag to the header
	int open_out(std::string filename){this->of = sam_open(filename.c_str(), "wb"); return sam_hdr_write(this->of, this->h);}
	//
	int write(){return sam_write1(this->of, this->h, this->r);}
	//
}; //end of BamFile class

class FastqFile: public HTSFile
{
public:
	gzFile fh;
	kseq_t* r;
	gzFile ofh;
	FastqFile(std::string filename){
		fh = gzopen(filename.c_str(),"r")
		r = kseq_init(fh);
	};
	~FastqFile(){
		kseq_destroy(r);
		gzclose(fh);
		if(ofh != Z_NULL){gzclose(ofh);}
	}
	int next(){return kseq_read(r);}
	std::string next_str(){return this->next() >= 0? std::string(this->r->seq.s): "";}
	readutils::CReadData get(){return readutils::CReadData(this->r);}
	void recalibrate(std::vector<int> qual){for(int i = 0; i < this->r->qual.l; ++i){this->r->qual.s[i] = (char)qual[i];}}
	int open_out(std::string filename){ofh = gzopen(filename.c_str(),"w"); return ofh == Z_NULL ? -1 : 0;}
	int write(){return gzputs(ofh, std::string(this->r->name.s) + "\n" + std::string(this->r->seq.s) + "\n" + std::string(this->r->comment.s) + "\n" + std::string(this->r->qual.s) + "\n");}
};

class KmerSubsampler{
public:
	HTSFile* file;
	minion::Random rng;
	std::bernoulli_distribution d;
	std::string readseq = "";
	std::vector<uint64_t> kmers;
	size_t cur_kmer = 0;
	int k;
	KmerSubsampler(HTSFile* file): KmerSubsampler(file, 32){}
	KmerSubsampler(HTSFile* file, int k): KmerSubsampler(file, k, .15){}
	KmerSubsampler(HTSFile* file, int k, double alpha): KmerSubsampler(file, k, alpha, minion::create_uint64_seed(minion::create_seed_seq())){}
	KmerSubsampler(HTSFile* file, int k, double alpha, uint64_t seed): file(file), k(k), d(alpha) {rng.seed(seed);}
	
	//return the next kmer
	//once the file is finished iterating and there are no remaining kmers,
	//return 0. That means you should check readseq.empty() if you get a 0!
	uint64_t next_kmer(){
		if(cur_kmer < kmers.length()){
			return kmers[cur_kmer++]; //return the current kmer and advance
		} else {
			readseq = file.next_str();
			if(readseq.empty()){
				return 0; //no more sequences
			} else {
				kmers = bloom::hash_seq(readseq, k); //get new vector of hashes
				cur_kmer = 0; //reset current kmer
				return this->next_kmer(); //try again
			}
		}
	}

	//return the next kmer that survives sampling
	//once there are no more kmers return 0.
	// check to see if readseq.empty() if you get a 0 result!
	uint64_t next(){
		uint64_t kmer = this->next_kmer();
		if(!readseq.empty()){
			if(d(rng)){ //sampled
				return kmer;
			}
			else{ //try again
				return this->next();
			}
		}
		return 0; // empty
	}
}

}