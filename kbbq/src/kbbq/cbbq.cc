#include "htsiter.hh"
#include "bloom.hh"
#include "covariateutils.hh"
#include "recalibrateutils.hh"

using htsiter, recalibrateutils, bloom, covariateutils

int main(int argc, char* argv[]){
	BamFile file(std::string(argv[0]));
	KmerSubsampler subsampler(&file);
	//load subsampled bf x
	kmer_cache_t subsampled = subsample_kmers(subsampler);

	//calculate thresholds
	//TODO

	//get trusted kmers bf using subsampled bf
	file = BamFile(std::string(argv[0]));
	kmer_cache_t trusted = find_trusted_kmers(file, subsampled, thresholds);

	//use trusted kmers to find errors
	file = BamFile(std::string(argv[0]));
	CCovariateData data = get_covariatedata(file, trusted);

	//recalibrate reads and write to file
	dq_t dqs = data.get_dqs();
	file = BamFile(std::string(argv[0]));
	recalibrate_and_write(file, dqs, "-");
}


