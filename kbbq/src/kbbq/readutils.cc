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

	CReadData::CReadData(kseq_t* fastqrecord, std::string rg = "", int second = 2, std::string namedelimiter = "_"){
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

	void CReadData::infer_read_errors(bloom::bloom_ary_t b, std::vector<int> thresholds, int k){
		std::array<std::vector<int>,2> overlapping = bloom::overlapping_kmers_in_bf(this->seq, b, k);
		std::vector<int> in = overlapping[0];
		std::vector<int> possible = overlapping[1];
		for(int i = 0; i < errors.length(); ++i){
			this->errors[i] = (in[i] <= thresholds[possible[i]] || this->qual[i] <= 6);
		}
	}

	int CReadData::correct_one(bloom::bloom_ary_t& t, int k){
		int best_fix_len = 0;
		char best_fix_base;
		size_t best_fix_pos = 0;
		for(size_t i = 0; i < this->seq.length(); ++i){
			std::string original_seq(this->seq); //copy original
			for(char c : {"A","C","G","T"}){
				original_seq[i] = c;
				start_pos = (size_t) std::min(std::max(i-k/2+1,0), this->seq.length()-k);
				int n_in = bloom::nkmers_in_bf(this->seq.substr(start_pos),t,k);
				if(n_in > best_fix_len){
					best_fix_base = c;
					best_fix_pos = i;
					best_fix_len = n_in;
				}
			}
		}
		this->seq[best_fix_pos] = best_fix_base;
		return best_fix_len;
	}

	//this is a chonky boi
	void CReadData::get_errors(bloom::bloom_ary_t& trusted, int k, int minqual = 6){
		std::string original_seq(this->seq);
		std::array<size_t,2> anchor = bloom::find_longest_trusted_seq(this->seq, trusted, k);
		if(anchor[0] == std::string::npos){ //no trusted kmers in this read.
			if(this->correct_one(trusted, k) == 0){
				return;
			} else {
				anchor = bloom::find_longest_trusted_seq(this->seq, trusted, k);
			}
		}
		if(anchor[0] == 0 && anchor[1] == std::string::npos){ //all kmers are trusted
			return;
		}
		//we're guaranteed to have a valid anchor now.
		bool corrected; //whether there were any corrections
		bool multiple; //whether there were any ties
		//right side
		for(int i = anchor[1] + 1; i < this->seq.length();){
			int start = i - k + 1; //seq containing all kmers that are affected
			// size_t end = i + k < this->seq.length() ? i + k : std::string::npos; //find_longest_fix will take care of this
			std::pair<std::vector<char>, int> lf = bloom::find_longest_fix(this->seq.substr(start), trusted, k);
			std::vector<char> fix = std::get<0>(lf);
			int fixlen = std::get<1>(lf);
			if(fixlen != 0){
				this->seq[i] = fix[0];
				this->errors[i] = true;
				i += fixlen;
				corrected = true;
				if(fix.size() > 1){
					multiple = true;
				}
			} else {
				//couldn't find a fix; skip ahead and try again if long enough
				i += k-1; //move ahead and make i = k-1 as the first base in the new kmer
				if(this->seq.length() - i + k <= (seq.length()/2) || this->seq.length() - i + k <= 2*k ){
					//sequence not long enough. end this side.
					break;
				}
			}
		}
		//left side
		//the bad base is at anchor[0]-1, then include the full kmer for that base.
		std::string sub = this->seq.substr(0, anchor[0] - 1 + k);
		std::string revcomped;
		for(auto it = sub.rbegin(); it != sub.rend(); it++){
				int c = seq_nt4_table[*it];
				char d = c < 4 ? seq_nt16_str[seq_nt16_table[3 - c]] : c;
				revcomped.push_back(d); //complement then turn to 4-bit encoded then turn to char
			}
		for(int i = anchor[0] - 1; i > 0;){ //iterating in original seq space
			int j = anchor[0] - 1 + k - 1 - i; //iterating in reversed space
			int end = i + k; //seq containing all kmers that are affected is [0, end) in original space
			//but [j -k + 1, npos) in reverse space. 
			std::string sub = revcomped.substr(j - k + 1); //get the right subsequence
			std::pair<std::vector<char>, int> lf = bloom::find_longest_fix(sub, trusted, k);
			std::vector<char> fix = std::get<0>(lf);
			int fixlen = std::get<1>(lf);
			if(fixlen != 0){
				int c = seq_nt4_table[fix[0]];
				char d = c < 4 ? seq_nt16_str[seq_nt16_table[3 - c]] : c;
				this->seq[j] = d;
				this->errors[j] = true;
				i -= fixlen;
				corrected = true;
				if(fix.size() > 1){
					multiple = true;
				}
			} else {
				//couldn't find a fix; skip ahead and try again if long enough
				i -= k+1;
				if(i + k <= (seq.length()/2) || i + k <= 2*k ){
					//sequence not long enough. end this side.
					break;
				}
			}
		}
		// check for overcorrection and fix it
		if(corrected){
			bool adjust = !multiple; //if there were any ties, don't adjust
			int ocwindow = 20;
			int threshold = 4;
			double occount = 0;
			//check for overcorrection
			std::vector<int> overcorrected_idx; //push_back overcorrected indices in order
			//then from overcorrected.begin() - k to overcorrected.end() + k should all be reset.
			for(int i = 0; i < this->seq.length(); ++i){
				if(this->errors[i] && seq_nt4_table[this->seq[i]] < 4){ //increment correction count
					if(this->qual[i] <= minqual){
						occount += 0.5;
					} else {
						++occount;
					}
				}
				if(i > ocwindow && this->errors[i-ocwindow] && seq_nt4_table[this->seq[i-ocwindow]] < 4){ //decrement count for not in window
					if(this->qual[i-ocwindow] <= minqual){
						occount -= 0.5;
					} else {
						--occount;
					}
				}
				//set threshold
				if(adjust && i >= ocwindow){
					threshold++;
				}
				if(adjust && i + ocwindow - 1 < this->seq.length()){
					threshold--;
				}
				//determine if overcorrected
				if(occount > threshold && this->errors[i]){
					overcorrected_idx.push_back(i);
				}

			}

			if(overcorrected_idx.size() > 0){
				int start = overcorrected_idx[0]-k+1; //the beginningmost position to check
				start = start >= 0? start : 0;
				int end = *overcorrected_idx.end()+k-1; //the endmost position to check
				end < this->seq.length() ? end : this->seq.length();
				//we start iteration but we need to unfix anything within k of an overcorrected position
				//OR within k of one of those fixed positions.
				for(int i = start; i < end; ++i){
					if(this->errors[i]){
						this->errors[i] = false;
						if(i-k+1 < start){ //go back a bit if we need to
							i = i-k+1 >= 0 ? i-k+1 : 0;
							start = i;
						}
						if(i+k-1 > end){ //change the end if we need to
							end = i+k-1 < this->seq.length() ? i+k-1 : this->seq.length();
						}
					}
				}
			}
		}
		this->seq = original_seq;
	}

	std::vector<int> CReadData::recalibrate(dq_t dqs, int minqual = 6){
		std::vector<int> recalibrated(this->seq.length());
		int rg = this->get_rg_int();
		for(int i = 0; i < this->seq.length(); ++i){
			int q = this->qual[i];
			if(q >= minqual){
				recalibrated[i] = dqs.meanq[rg] + dqs.rgdq[rg] + dqs.qscoredq[rg][q] +
				dqs.cycledq[rg][q][this->second][i];
				if(i > 0){
					int first = seq_nt4_table[this->seq[i-1]];
					int second = seq_nt4_table[this->seq[i]];
					if(first < 4 && second < 4){
						int8_t dinuc = dinuc_to_int(first, second);
						recalibrated[i] += dqs.dinucdq[rg][q][dinuc];
					}
				}
			}
		}
		return recalibrated;
	}
}

