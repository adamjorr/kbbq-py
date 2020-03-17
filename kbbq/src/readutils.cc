#include "readutils.hh"
#include "include/htslib/htslib/sam.h"
#include "include/htslib/htslib/kseq.h"

#include <algorithm>

namespace readutils{

	static int8_t complement[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

	CReadData::CReadData(bam1_t* bamrecord){
		unsigned char* s = bam_get_seq(bamrecord);
		for(int i = 0; i < bamrecord->core.l_qseq; ++i){
			char c;
			if(bam_is_rev(bamrecord)){
				c = seq_nt16_str[complement[bam_seqi(s, i)]];
			}
			else{
				c = seq_nt16_str[bam_seqi(s, i)];
			}
			this->seq.push_back(c);
		}
		this->qual.assign(bam_get_qual(bamrecord), bam_get_qual(bamrecord) + bamrecord->core.l_qseq);
		if(bam_is_rev(bamrecord)){
			std::reverse(this->seq.begin(), this->seq.end());
			std::reverse(this->qual.begin(), this->qual.end());
		}
		this->skips.resize(bamrecord->core.l_qseq);
		this->name = bam_get_qname(bamrecord);
		this->rg = bam_aux2Z(bam_aux_get(bamrecord, "RG"));
		if(rg_to_pu.count(this->rg) == 0){
			rg_to_int[this->rg] = rg_to_int.size();
			rg_to_pu[this->rg] = rg; //when loaded from the header this is actually a PU.
		}
		this->second = bamrecord->core.flag & 0x80;
		this->errors.resize(bamrecord->core.l_qseq);
	}
	// if second is >1, that means infer.

	CReadData::CReadData(kseq_t* fastqrecord, std::string rg, int second, std::string namedelimiter){
		this->seq = std::string(fastqrecord->seq.s);
		// this->qual.assign(fastqrecord->qual.s, fastqrecord->qual.l);
		std::string quals(fastqrecord->qual.s);
		this->qual.assign(quals.begin(), quals.end());
		this->skips.resize(this->seq.length());

		std::string fullname(fastqrecord->name.s);
		size_t delim = fullname.find(namedelimiter);
		std::string first_name = fullname.substr(0, delim);

		while(rg == "" && delim != std::string::npos){
			// if we need to find rg
			fullname = fullname.substr(delim); //get the right part
			delim = fullname.find(namedelimiter); //reset the delimiter
			if(fullname.substr(0,3) == "RG:"){
				size_t last_colon = fullname.find_last_of(":", delim); // delim is last char to search
				rg = fullname.substr(last_colon, delim);
			}
		}
		this->rg = rg;

		if(second > 1){
			std::string tail = first_name.substr(first_name.length() - 2);
			second = (tail == "/2");
			if(second || tail == "/1"){
				first_name = first_name.substr(0, first_name.length() - 2);
			}
		}
		this->name = first_name;
		this->second = second;
		this->errors.resize(this->seq.length());

		if(rg_to_pu.count(this->rg) == 0){
			rg_to_int[this->rg] = rg_to_int.size();
			rg_to_pu[this->rg] = rg; //when loaded from the header this is actually a PU.
		}
	}

	void CReadData::load_rgs_from_bamfile(bam_hdr_t* header){
		std::string hdrtxt(header->text);
		size_t linedelim = hdrtxt.find('\n');
		// kstring_t val;
		// while (sam_hdr_find_tag_pos(header, "RG", i, "ID", val) == 0){
		// 	std::string rgid(val.s);
		// 	sam_hdr_find_tag_pos(header, "RG", i, "PU", val);
		// 	std::string pu(val.s);
		// 	if(rg_to_pu.count(rgid) == 0){
		// 		rg_to_int[rgid] = rg_to_int.size();
		// 		rg_to_pu[rgid] = pu;
		// 	}
		// }
		while(linedelim != std::string::npos){
			if(hdrtxt.substr(0,3) == "@RG"){
				std::string id("");
				std::string pu("");
				std::string line = hdrtxt.substr(0, linedelim);
				size_t tokendelim = line.find_first_of("\t ");
				while(tokendelim != std::string::npos){
					if(line.substr(0,3) == "ID:"){
						//id
						id = line.substr(4, tokendelim);
					} else if (line.substr(0, 3) == "PU:"){
						//pu
						pu = line.substr(4, tokendelim);
					}
					line = line.substr(tokendelim);
					tokendelim = line.find_first_of("\t ");
				}
				if(id != ""){
					rg_to_int[id] = rg_to_int.size();
					rg_to_pu[id] = pu;
				}
			}
			hdrtxt = hdrtxt.substr(linedelim);
			linedelim = hdrtxt.find('\n');
		}
	}

	std::string CReadData::str_qual(){
		std::string str_qual;
		for(size_t i = 0; i < this->qual.size(); i++){
			str_qual.push_back(this->qual[i] + 33);
		}
		return str_qual;
	}

	std::string CReadData::canonical_name(){
		std::string suffix;
		if (this->second){
			suffix = "/2";
		}
		else{
			suffix = "/1";
		}
		return this->name + suffix;
	}

	inline int CReadData::get_rg_int(){
		return this->rg_to_int[this->rg];
	}

	inline std::string CReadData::get_pu(){
		return this->rg_to_pu[this->rg];
	}

	std::vector<bool> CReadData::not_skipped_errors(){
		std::vector<bool> unskipped = this->skips;
		unskipped.flip();
		for(size_t i = 0; i < unskipped.size(); i++){
			unskipped[i] = (unskipped[i] && this->errors[i]);
		}
		return unskipped;
	}
}

