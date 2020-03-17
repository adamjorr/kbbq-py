#ifndef READUTILS_H
#define READUTILS_H

#include <vector>
#include <unordered_map>
#include <string>
//
#include "include/htslib/htslib/sam.h"
// #include "include/htslib/htslib/kseq.h"
#include "include/htslib/htslib/bgzf.h"
#include "pysam_stream.h"
// KSEQ_INIT(BGZF *, bgzf_read)

namespace readutils{
	class CReadData{
		private:
			static std::unordered_map<std::string, std::string> rg_to_pu;
			static std::unordered_map<std::string, int> rg_to_int;
		public:
			CReadData(){}
			CReadData(bam1_t* bamrecord);
			// hello?
			CReadData(kseq_t* fastqrecord, std::string rg = "", int second = 2, std::string namedelimiter = "_");
			std::string seq;
			std::vector<int> qual;
			std::vector<bool> skips;
			std::string name;
			std::string rg;
			bool second;
			std::vector<bool> errors;
			std::string str_qual();
			std::string canonical_name();
			int get_rg_int();
			std::string get_pu();
			std::vector<bool> not_skipped_errors();
			static void load_rgs_from_bamfile(bam_hdr_t* header);

	};
}

#endif
