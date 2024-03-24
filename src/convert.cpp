#include "convert.h"
#include "util.h"
#include <typeinfo>
#include <htslib/faidx.h>

//Path query Dot
string PathqueryDot(string filename){
	string newfilename;
	size_t lastDotPos = filename.find_last_of('.');
	if(lastDotPos != string::npos and filename.substr(lastDotPos + 1).compare("bed") == 0)	return filename;
	else
		newfilename = filename.substr(0, lastDotPos) + ".bed";

	return newfilename;
}
//Path query backslash
string Pathquerybackslash(string filename){
	string newfilename;
	size_t lastSlashPos = filename.find_last_of('/');
	if (lastSlashPos != std::string::npos and lastSlashPos < filename.length() - 1) {
		newfilename =  filename.substr(lastSlashPos + 1);
		return newfilename;
	}

	return filename;
}

// convert data
void convertBed(const string &infilename, const string &outfilename, const string &reffilename, string &mate_filename, string &snv_filename){
	ifstream infile;
	string line, chrname, chrname2, sv_type_str, ref_seq, alt_seq, reg_str;
	vector<string> str_vec, sv_type_vec;
	size_t start_pos, endpos, start_pos2, endpos2;
	int32_t sv_len = 0, refseq_len_tmp;
	vector<int32_t> sv_len_vec;
	vector<SV_item*> sv_item_vec;
	SV_item *sv_item;
	char *seq;
	faidx_t *fai;

	fai = fai_load(reffilename.c_str());

	infile.open(infilename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << infilename << endl;
		exit(1);
	}

	// convert
	while(getline(infile, line)){
		if(line.size()){
			if(line.at(0)!='#'){
				str_vec = split(line, "\t");

				chrname = str_vec.at(0);
				start_pos = stoi(str_vec.at(1));
				endpos = stoi(str_vec.at(2));

				// get sv type and sv_len
				sv_type_vec = getSVType(str_vec);
				//sv_len_vec = getSVLen(str_vec, sv_type_str);
				sv_len_vec = getSVLen(str_vec, sv_type_vec);

				sv_type_str = sv_type_vec.at(0);
				sv_len = sv_len_vec.at(0);

				ref_seq = str_vec.at(3);
				alt_seq = str_vec.at(4);

				if(sv_type_str.compare("DEL")==0 and sv_len_vec.at(0)>0) sv_len = -sv_len;
				if(sv_type_str.compare("TRA")==0 or sv_type_str.compare("BND")==0){
					chrname2 = str_vec.at(3);
					start_pos2 = stoi(str_vec.at(4));
					endpos2 = stoi(str_vec.at(5));
					ref_seq = alt_seq = "-";
				}else{
					chrname2 = "";
					start_pos2 = endpos2 = 0;
				}

				if(ref_seq.compare("N")==0 or ref_seq.compare(".")==0 or ref_seq.compare("0")==0){
					reg_str = chrname + ":" + to_string(start_pos) + "-" + to_string(endpos);
					seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
					if(strlen(seq)>=Max_SeqLen) ref_seq="-";
					else ref_seq = seq;
					free(seq);
				}
				if(alt_seq.compare("N")==0 or alt_seq.compare(".")==0 or alt_seq.compare("0")==0 or alt_seq.size()>=Max_SeqLen)
					alt_seq = "-";

				sv_item = allocateSVItem(chrname, start_pos, endpos, chrname2, start_pos2, endpos2, sv_type_str, sv_len, ref_seq, alt_seq);
				sv_item_vec.push_back(sv_item);
			}
		}
	}
	infile.close();
	fai_destroy(fai);

	// remove mate sv items
	if(mate_filename.size()) removemateSVItems(mate_filename, sv_item_vec);

	// remove snv items
	if(snv_filename.size()) removeSNVItems(snv_filename, sv_item_vec);

	// save to file
	output2File(outfilename, sv_item_vec, outConvertScreenFile);

	// release memory
	releaseSVItemVec(sv_item_vec);
}

//Determine if it is present in chromosome collections
bool isExistChromosomeSet(string &chrname){
	auto it = find(chromosomeSet.begin(), chromosomeSet.end(), chrname);
	if (it != chromosomeSet.end()) 	return true;
	else	return false;
}

//To determine if it's decoy
bool isDecoyChr(string &chrname){
	bool flag = false;
	size_t pos;

	pos = chrname.find(DECOY_PREFIX);
	if(pos==0) flag = true;
	else{
		pos = chrname.find(DECOY_PREFIX2);
		if(pos==0) flag = true;
	}

	return flag;
}
// convert vcf result
void convertVcf(const string &infilename, const string &outfilename, const string &reffilename, string &mate_filename, string &snv_filename){
	ifstream infile;
	string line, chrname, start_pos_str, endpos_str, endpos_str_1, endpos_str_2, chrname2, start_pos_str2, endpos_str2, sv_type_str, sv_type_str1, sv_type_str2, sv_len_str, sv_len_str1, sv_len_str2, bnd_str, bnd_pos_str, str_tmp;
	string ref_seq, alt_seq, reg_str;
	size_t start_pos, endpos, endpos_1, endpos_2, start_pos2, endpos2;
	size_t i, end_pos_ref, end_pos_ref1, end_pos_ref2;
	int32_t sv_len, sv_len1, sv_len2, ref_seq_len, refseq_len_tmp;
	bool is_svlen_flag;
	vector<int32_t> sv_len_vec;
	vector<string> str_vec, sv_type_vec, info_vec, sub_info_vec, bnd_pos_vec;
	vector<SV_item*> sv_item_vec;
	SV_item *sv_item, *sv_item1, *sv_item2;
	char *seq;
	faidx_t *fai;

	fai = fai_load(reffilename.c_str());
	is_svlen_flag = false;

	infile.open(infilename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << infilename << endl;
		exit(1);
	}

	// convert
	while(getline(infile, line)){
		if(line.size()){
			if(line.at(0)!='#'){

				str_vec = split(line, "\t");

				// get locations
				chrname = str_vec.at(0);
				if(chromosomeSet.size()>0){
					if(isExistChromosomeSet(chrname) == false) continue;
				}else{
					if(isDecoyChr(chrname) == true) continue;
				}
				start_pos_str = str_vec.at(1);
//				if(start_pos_str!="39059")
//					continue;
//				cout<<"start_pos_str:"<<start_pos_str<<endl;
				ref_seq_len = str_vec.at(3).size();
				ref_seq = str_vec.at(3);
				alt_seq = str_vec.at(4);

				if(alt_seq.at(0) == '<') alt_seq = "-";

				chrname2 = endpos_str = sv_type_str = sv_len_str = "";
				info_vec = split(str_vec.at(7), ";");

				for(i=0; i<info_vec.size(); i++){
					sub_info_vec = split(info_vec.at(i), "=");
					if(sub_info_vec.at(0).compare("CHR2")==0)
						chrname2 = sub_info_vec.at(1);
					if(sub_info_vec.at(0).compare("END")==0)
						endpos_str = sub_info_vec.at(1);
					if(sub_info_vec.at(0).compare("SVTYPE")==0)
						sv_type_str = sub_info_vec.at(1);
					if(sub_info_vec.at(0).compare("SVLEN")==0){
						sv_len_str = sub_info_vec.at(1);

						if(sv_len_str.find_first_not_of("-0123456789")!=string::npos) sv_len_str = "";
						else{ // ignore too large invalid variant length, e.g. SVLEN=999999999 in Sniffles output
							sv_len = stoi(sv_len_str);
							is_svlen_flag = true;
							if(sv_len>MAX_VALID_SVLEN) { sv_len_str = ""; sv_len = 0;}
						}
					}
				}

				//Skip mutations with incomplete information
//				if(sv_type_str.size()==0 and sv_len_str.size()==0 and endpos_str.size()==0)	continue;

				// get sv type and length
				if(sv_type_str.size()==0){
					sv_type_vec = getSVType(str_vec);
					if(sv_type_vec.size()==1){
						sv_type_str = sv_type_vec.at(0);
					}else{
							sv_type_str1 = sv_type_vec.at(0);
							sv_type_str2 = sv_type_vec.at(1);
					}
				}
				if(sv_len_str.size()==0){
					sv_len_vec = getSVLen(str_vec, sv_type_str);
					if(sv_len_vec.size()==1){
						sv_len_str = to_string(sv_len_vec.at(0));
					}else{
						if(sv_len_vec.size()!=0){
							sv_len_str1 = to_string(sv_len_vec.at(0));
							sv_len_str2 = to_string(sv_len_vec.at(1));
						}else continue;
					}
				}
				if(endpos_str.size()==0){
					if(sv_len_vec.size()==1){
							end_pos_ref = stoi(start_pos_str) + ref_seq_len - 1;
							endpos_str = to_string(end_pos_ref);
					}else{
						if(sv_len_str.size()>0)	endpos_str = to_string(stoi(start_pos_str) + stoi(sv_len_str) - 1);
						else endpos_str = to_string(stoi(start_pos_str) + abs(sv_len_vec.at(0)) - 1);
						end_pos_ref1 = stoi(start_pos_str) + ref_seq_len - 1;
						end_pos_ref2 = stoi(start_pos_str) + ref_seq_len - 1;
						endpos_str_1 = to_string(end_pos_ref1);
						endpos_str_2 = to_string(end_pos_ref2);
					}
				}

				if(sv_type_str.size()>0){
					if(sv_type_str.compare("TRA")==0 or sv_type_str.compare("BND")==0){ // TRA or BND
						if(sv_type_str.compare("BND")==0){ // BND
							bnd_str = str_vec.at(4);
							if(bnd_str.at(0)==']' or bnd_str.at(0)=='[')
								bnd_pos_str = bnd_str.substr(1, bnd_str.size()-3);
							else if(isBase(bnd_str.at(0)))
								bnd_pos_str = bnd_str.substr(2, bnd_str.size()-3);
							bnd_pos_vec = split(bnd_pos_str, ":");
							chrname2 = bnd_pos_vec.at(0);
							endpos_str = bnd_pos_vec.at(1);
						}

						start_pos_str2 = endpos_str;
						endpos_str2 = start_pos_str2;
						endpos_str = start_pos_str;
					}else { // INS, DEL, INV, DUP, UNC
						start_pos_str2 = "0";
						endpos_str2 = "0";
					}
					if(sv_len_vec.size()==1 or is_svlen_flag == true){
						start_pos = stoi(start_pos_str);
						endpos = stoi(endpos_str);
						start_pos2 = stoi(start_pos_str2);
						endpos2 = stoi(endpos_str2);
						if(sv_len_str.size()>0) sv_len = stoi(sv_len_str);
						else sv_len = 0;

						if(sv_type_str.compare("DEL")==0) sv_len = abs(sv_len);

						if(ref_seq.compare("N")==0 or ref_seq.compare(".")==0 or ref_seq.compare("0")==0){
							if(endpos < (size_t)stoi(start_pos_str))	endpos = (size_t)stoi(start_pos_str);
							reg_str = chrname + ":" + start_pos_str + "-" + to_string(endpos);
							seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
							if(seq == nullptr or strlen(seq) >= Max_SeqLen ) ref_seq = "-";
							else ref_seq = seq;
							free(seq);
						}
						if(alt_seq.size()>=Max_SeqLen or sv_type_str.compare("TRA")==0 or sv_type_str.compare("BND")==0 or alt_seq.compare(".")==0) alt_seq="-";
						sv_item = allocateSVItem(chrname, start_pos, endpos, chrname2, start_pos2, endpos2, sv_type_str, sv_len, ref_seq, alt_seq);
						sv_item_vec.push_back(sv_item);
					}else{
						start_pos = stoi(start_pos_str);
						endpos_1 = stoi(endpos_str_1);
						endpos_2 = stoi(endpos_str_2);
						start_pos2 = stoi(start_pos_str2);
						endpos2 = stoi(endpos_str2);
						sv_len1 = stoi(sv_len_str1);
						sv_len2 = stoi(sv_len_str2);

						if(sv_type_str1.compare("DEL")==0) sv_len1 = abs(sv_len1);
						if(sv_type_str2.compare("DEL")==0) sv_len2 = abs(sv_len2);

						ref_seq = alt_seq = "-";

						sv_item1 = allocateSVItem(chrname, start_pos, endpos_1, chrname2, start_pos2, endpos2, sv_type_str1, sv_len1, ref_seq, alt_seq);
						sv_item2 = allocateSVItem(chrname, start_pos, endpos_2, chrname2, start_pos2, endpos2, sv_type_str2, sv_len2, ref_seq, alt_seq);
						sv_item_vec.push_back(sv_item1);
						sv_item_vec.push_back(sv_item2);
					}
				}else{
					cout << line << endl;
					cout << "missing SVTYPE information" << endl;
				}
			}
		}
	}
	infile.close();
	fai_destroy(fai);

	// remove mate sv items
	if(mate_filename.size()) removemateSVItems(mate_filename, sv_item_vec);

	// remove snv items
	if(snv_filename.size()) removeSNVItems(snv_filename, sv_item_vec);

	// save to file
	output2File(outfilename, sv_item_vec, outConvertScreenFile);

	// release memory
	releaseSVItemVec(sv_item_vec);
}

// convert data
void convertCsv(const string &infilename, const string &outfilename, const string &reffilename, string &mate_filename, string &snv_filename){
	ifstream infile;
	string line, chrname, chrname2, sv_type_str, ref_seq, alt_seq, reg_str;
	vector<string> str_vec, sv_type_vec;
	size_t start_pos, endpos, start_pos2, endpos2;
	int32_t sv_len = 0, refseq_len_tmp;
	vector<int32_t> sv_len_vec;
	vector<SV_item*> sv_item_vec;
	SV_item *sv_item;
	char *seq;
	faidx_t *fai;

	fai = fai_load(reffilename.c_str());

	infile.open(infilename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << infilename << endl;
		exit(1);
	}

	// convert
	while(getline(infile, line)){
		if(line.size()>0){

			if(line.at(0)!='#'){
				str_vec = split(line, ",");
				if(str_vec.at(0).size()>0){

					chrname = str_vec.at(0);
					if(isDigitString(str_vec.at(1))==false or isDigitString(str_vec.at(2))==false)
						continue;
					start_pos = stoi(str_vec.at(1));
					endpos = stoi(str_vec.at(2));

					// get sv type and length
					sv_type_vec = getSVType(str_vec);
					sv_len_vec = getSVLen(str_vec, sv_type_vec);

					ref_seq = str_vec.at(3);
					alt_seq = str_vec.at(4);

					sv_type_str = sv_type_vec.at(0);
					sv_len = sv_len_vec.at(0);
					//if(sv_type_str.compare("DEL")==0 and sv_len_vec.at(0)>0) sv_len = -sv_len;
					//if(sv_type_str.compare("TRA")==0 or sv_type_str.compare("BND")==0){
					if(sv_type_str.compare("DEL")==0 and sv_len_vec.at(0)>0) sv_len = -sv_len;
					if(sv_type_str.compare("TRA")==0 or sv_type_str.compare("BND")==0){
						chrname2 = str_vec.at(3);
						start_pos2 = stoi(str_vec.at(4));
						endpos2 = stoi(str_vec.at(5));
					}else{
						chrname2 = "";
						start_pos2 = endpos2 = 0;
					}

					if(ref_seq.compare("N")==0 or ref_seq.compare(".")==0 or ref_seq.compare("0")==0){
						reg_str = chrname + ":" + to_string(start_pos) + "-" + to_string(endpos);
						seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
						if(strlen(seq)>=Max_SeqLen) ref_seq="-";
						else ref_seq = seq;
						free(seq);
					}
					if(alt_seq.compare("N")==0 or alt_seq.compare(".")==0 or alt_seq.compare("0")==0 or alt_seq.size()>Max_SeqLen)
						alt_seq = "-";

					//sv_item = allocateSVItem(chrname, start_pos, endpos, chrname2, start_pos2, endpos2, sv_type_str, sv_len);
					sv_item = allocateSVItem(chrname, start_pos, endpos, chrname2, start_pos2, endpos2, sv_type_str, sv_len, ref_seq, alt_seq);
					sv_item_vec.push_back(sv_item);
				}
			}
		}
	}
	infile.close();

	// remove mate sv items
	if(mate_filename.size()) removemateSVItems(mate_filename, sv_item_vec);

	// remove snv items
	if(snv_filename.size()) removeSNVItems(snv_filename, sv_item_vec);

	// save to file
	output2File(outfilename, sv_item_vec, outConvertScreenFile);

	// release memory
	releaseSVItemVec(sv_item_vec);
}

// convert data
void convertNm(const string &infilename, const string &outfilename, const string &reffilename, string &mate_filename, string &snv_filename){
	ifstream infile;
	string line, chrname, chrname2, sv_type_str, ref_seq, alt_seq, reg_str;
	vector<string> str_vec, sv_type_vec;
	int32_t sv_len, ref_seq_len, start_pos, endpos, start_pos2, endpos2, refseq_len_tmp;
	vector<int32_t> sv_len_vec;
	vector<SV_item*> sv_item_vec;
	SV_item *sv_item;
	char *seq;
	faidx_t *fai;

	fai = fai_load(reffilename.c_str());

	sv_len = 0;

	infile.open(infilename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << infilename << endl;
		exit(1);
	}

	// convert
	while(getline(infile, line)){
		if(line.size()>0){
			if(line.at(0)!='#'){
				str_vec = split(line, "\t");

				chrname = str_vec.at(0);
				start_pos = stoi(str_vec.at(1));

				ref_seq_len = str_vec.at(3).size();
				endpos = stoi(str_vec.at(1)) + ref_seq_len - 1; // endPos

				ref_seq = str_vec.at(3);
				alt_seq = str_vec.at(4);

				if(alt_seq.at(0) == '<') alt_seq = "-";

				// get sv type
				sv_type_vec = getSVType(str_vec);

				// get sv length
				//sv_len_vec = getSVLen(str_vec, sv_type_str);
				sv_len_vec = getSVLen(str_vec, sv_type_vec);

				sv_type_str = sv_type_vec.at(0);
				sv_len = sv_len_vec.at(0);
				if(sv_type_str.compare("DEL")==0 and sv_len_vec.at(0)>0) sv_len = -sv_len;
				if(sv_type_str.compare("TRA")==0 or sv_type_str.compare("BND")==0){
					chrname2 = str_vec.at(3);
					start_pos2 = stoi(str_vec.at(4));
					endpos2 = stoi(str_vec.at(5));
				}else{
					chrname2 = "";
					start_pos2 = endpos2 = 0;
				}
//				chrname2 = "";
//				start_pos2 = endpos2 = 0;

				if(sv_type_str.compare("DEL")==0) sv_len = -sv_len;

				if(ref_seq.compare("N")==0 or ref_seq.compare(".")==0 or ref_seq.compare("0")==0){
					reg_str = chrname + ":" + to_string(start_pos) + "-" + to_string(endpos);
					seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
					if(strlen(seq)>=Max_SeqLen) ref_seq="-";
					else ref_seq = seq;
					free(seq);
				}
				if(alt_seq.size()>=Max_SeqLen) alt_seq="-";

				sv_item = allocateSVItem(chrname, start_pos, endpos, chrname2, start_pos2, endpos2, sv_type_str, sv_len, ref_seq, alt_seq);
				sv_item_vec.push_back(sv_item);
			}
		}
	}
	infile.close();
	fai_destroy(fai);

	// remove mate sv items
	if(mate_filename.size()) removemateSVItems(mate_filename, sv_item_vec);

	// remove snv items
	if(snv_filename.size()) removeSNVItems(snv_filename, sv_item_vec);

	// save to file
	output2File(outfilename, sv_item_vec, outConvertScreenFile);

	// release memory
	releaseSVItemVec(sv_item_vec);
}

// get sv type
vector<string> getSVType(vector<string> &str_vec){
	string sv_type_str, sv_type_str1, sv_type_str2, str_tmp;
	size_t i, j;
	bool is_seq_flag_ref, is_seq_flag_alt, is_comma_flag;
	vector<string> comma_str_vec;
	vector<string> info_vec, sv_type_vec;
	vector<int32_t> sv_len_vec;
	int32_t sv_type_pos, SVTYPE_pos, vec_idx, str_pos, start_pos, end_pos;

	sv_type_str = "";
	for(i=3; i<str_vec.size(); i++){
		str_tmp = str_vec.at(i);
		if(str_tmp.compare("INS")==0 or str_tmp.compare("DEL")==0 or str_tmp.compare("DUP")==0 or str_tmp.compare("INV")==0 or str_tmp.compare("TRA")==0 or str_tmp.compare("BND")==0 or str_tmp.compare("INVTRA")==0 or str_tmp.compare("MIX")==0 or str_tmp.compare("MNP")==0 or str_tmp.compare("SNV")==0){
			sv_type_str = str_tmp;
		}else if(str_tmp.compare("<INS>")==0 or str_tmp.compare("<DEL>")==0 or str_tmp.compare("<DUP>")==0 or str_tmp.compare("<INV>")==0 or str_tmp.compare("<TRA>")==0 or str_tmp.compare("<BND>")==0 or str_tmp.compare("<INVTRA>")==0 or str_tmp.compare("<MIX>")==0 or str_tmp.compare("<MNP>")==0 or str_tmp.compare("<SNV>")==0){
			sv_type_str = str_tmp.substr(1, str_tmp.size()-2);
		}else if(str_tmp.compare("insertion")==0 or str_tmp.compare("deletion")==0 or str_tmp.compare("duplication")==0 or str_tmp.compare("inversion")==0 or str_tmp.compare("translocation")==0 or str_tmp.compare("snv")==0){
			sv_type_str = str_tmp;
		}
		if(sv_type_str.size()>0) {
			sv_type_vec.push_back(sv_type_str);
			break;
		}
	}

	if(sv_type_str.size()==0){
		vec_idx = -1; str_pos = -1;
		for(i=3; i<str_vec.size(); i++){
			str_tmp = str_vec.at(i);
			sv_type_pos = str_tmp.find("sv_type=");
			if(sv_type_pos!=-1){
				vec_idx = i;
				str_pos = sv_type_pos + 8;
				break;
			}else{
				SVTYPE_pos = str_tmp.find("SVTYPE=");
				if(SVTYPE_pos!=-1){
					vec_idx = i;
					str_pos = SVTYPE_pos + 7;
					break;
				}
			}
		}

		if(vec_idx!=-1 and str_pos!=-1){ // found
			str_tmp = str_vec.at(vec_idx);
			start_pos = str_pos; end_pos = -1;
			for(j=start_pos; j<str_tmp.size(); j++){ // search ';'
				if(str_tmp.at(j)==';'){
					end_pos = j - 1;
					break;
				}
			}
			if(end_pos==-1) end_pos = str_tmp.size() - 1;
			sv_type_str = str_tmp.substr(start_pos, end_pos-start_pos+1);
			if(sv_type_str.size()>0)
				sv_type_vec.push_back(sv_type_str);
		}
	}

	if(sv_type_str.size()==0){// not found
		// check sequence
		is_seq_flag_ref = isSeq(str_vec.at(3));
		if(is_seq_flag_ref){
			is_seq_flag_alt = isSeq(str_vec.at(4));
			if(is_seq_flag_alt){
				if(str_vec.at(3).size()==1 and str_vec.at(4).size()==1)sv_type_str = "SNV";
				if(str_vec.at(3).size()<str_vec.at(4).size())sv_type_str = "INS";
				if(str_vec.at(3).size()>str_vec.at(4).size())sv_type_str = "DEL";
				sv_type_vec.push_back(sv_type_str);
			}else{
				str_tmp = str_vec.at(4);
				is_comma_flag = isComma(str_vec.at(4));
				if(is_comma_flag){
					comma_str_vec = split(str_vec.at(4), ",");
					if(str_vec.at(3).size()<comma_str_vec.at(0).size())sv_type_str1 = "INS";
					if(str_vec.at(3).size()>comma_str_vec.at(0).size())sv_type_str1 = "DEL";
					if(str_vec.at(3).size()<comma_str_vec.at(1).size())sv_type_str2 = "INS";
					if(str_vec.at(3).size()>comma_str_vec.at(1).size())sv_type_str2 = "DEL";
					sv_type_vec.push_back(sv_type_str1);
					sv_type_vec.push_back(sv_type_str2);
				}
			}
		}
	}

	return sv_type_vec;
}

// get SV length
vector<int32_t> getSVLen(vector<string> &str_vec, string &sv_type){
	int32_t sv_len, sv_len1, sv_len2, vec_idx, str_pos, sv_length_pos, SVLEN_pos, start_pos, end_pos, sv_END_pos;
	size_t i, j;
	string str_tmp, sv_len_str;
	bool flag, is_seq_flag_ref, is_seq_flag_alt, is_comma_flag;
	vector<string> comma_str_vec;
	vector<int32_t> sv_len_vec;

	flag = false;
	sv_len = 0;
	if(sv_type.compare("TRA")!=0 and sv_type.compare("BND")!=0){ // not TRA
		// sv_length, or SVLEN
		vec_idx = str_pos = -1;
		for(i=3; i<str_vec.size(); i++){
			str_tmp = str_vec.at(i);
			sv_length_pos = str_tmp.find("sv_length=");
			if(sv_length_pos!=-1){
				vec_idx = i;
				str_pos = sv_length_pos + 10;
				break;
			}else{
				SVLEN_pos = str_tmp.find("SVLEN=");
				if(SVLEN_pos!=-1){
					vec_idx = i;
					str_pos = SVLEN_pos + 6;
					break;
				}
			}
		}

		if(vec_idx!=-1 and str_pos!=-1){ // found
			str_tmp = str_vec.at(vec_idx);
			start_pos = str_pos; end_pos = -1;
			for(j=start_pos; j<str_tmp.size(); j++){ // search ';'
				if(str_tmp.at(j)==';'){
					end_pos = j - 1;
					break;
				}
			}
			if(end_pos==-1) end_pos = str_tmp.size() - 1;
			sv_len_str = str_tmp.substr(start_pos, end_pos-start_pos+1);
			if(sv_len_str.find_first_not_of("-0123456789")==string::npos) {
				sv_len = stoi(sv_len_str);
				sv_len_vec.push_back(sv_len);
				flag = true;
				if(sv_len>=MAX_VALID_SVLEN) { sv_len_str = ""; sv_len = 0; flag = false; }
			}else { sv_len_str = ""; flag = false; }
		}

		// 'SVLEN' was not detected, then further detect 'END'
		if(flag==false){
			vec_idx = str_pos = -1;
			for(i=3; i<str_vec.size(); i++){
				str_tmp = str_vec.at(i);
				sv_END_pos = str_tmp.find("END=");
				if(sv_END_pos!=-1){
					vec_idx = i;
					str_pos = sv_END_pos + 4;
					break;
				}
			}

			if(vec_idx!=-1 and str_pos!=-1){ // found
				str_tmp = str_vec.at(vec_idx);
				start_pos = str_pos; end_pos = -1;
				for(j=start_pos; j<str_tmp.size(); j++){ // search ';'
					if(str_tmp.at(j)==';'){
						end_pos = j - 1;
						break;
					}
				}
				if(end_pos==-1) end_pos = str_tmp.size() - 1;
				sv_len_str = str_tmp.substr(start_pos, end_pos-start_pos+1);
				if(sv_len_str.find_first_not_of("-0123456789")==string::npos) { sv_len = stoi(sv_len_str) - stoi(str_vec.at(1)) + 1; sv_len_vec.push_back(sv_len); flag = true; }
			}else{ // not found
				// check sequence
				is_seq_flag_ref = isSeq(str_vec.at(3));
				if(is_seq_flag_ref){
					is_seq_flag_alt = isSeq(str_vec.at(4));
					if(is_seq_flag_alt){
						sv_len = (int32_t)str_vec.at(4).size() - (int32_t)str_vec.at(3).size();
						sv_len_vec.push_back(sv_len);
						flag = true;
					}else{
						str_tmp = str_vec.at(4);
						is_comma_flag = isComma(str_vec.at(4));
						if(is_comma_flag){
							comma_str_vec = split(str_vec.at(4), ",");
							sv_len1 = (int32_t)comma_str_vec.at(0).size() - (int32_t)str_vec.at(3).size();
							sv_len2 = (int32_t)comma_str_vec.at(1).size() - (int32_t)str_vec.at(3).size();
							sv_len_vec.push_back(sv_len1);
							sv_len_vec.push_back(sv_len2);
							flag = true;
						}else{
							if(str_tmp.at(0)=='<' and str_tmp.at(str_tmp.size()-1)=='>'){
								sv_len = 0;
								flag = true;
							}
						}
					}
				}

				// no sequence
				if(flag==false){
					start_pos = stoi(str_vec.at(1));
					end_pos = stoi(str_vec.at(2));
					sv_len = end_pos - start_pos + 1;
					sv_len_vec.push_back(sv_len);
					flag = true;
				}
			}
		}
		if(flag==false) cout << "missing SVLEN information" << endl;
	}else sv_len_vec.push_back(sv_len);

	return sv_len_vec;
}

// get SV length
vector<int32_t> getSVLen(vector<string> &str_vec, vector<string> &sv_type_vec){
	int32_t sv_len, sv_len1, sv_len2, vec_idx, str_pos, sv_length_pos, SVLEN_pos, start_pos, end_pos, sv_END_pos;
	size_t i, j, n;
	string str_tmp, sv_len_str, sv_type_str;
	bool flag, is_seq_flag_ref, is_seq_flag_alt, is_comma_flag;
	vector<string> comma_str_vec;
	vector<int32_t> sv_len_vec;

	for(n=0; n<sv_type_vec.size(); n++){
		flag = false;
		sv_len = 0;
		sv_type_str = sv_type_vec.at(n);
		if(sv_type_str.compare("TRA")!=0 and sv_type_str.compare("BND")!=0){ // not TRA
			// sv_length, or SVLEN
			vec_idx = str_pos = -1;
			for(i=3; i<str_vec.size(); i++){
				str_tmp = str_vec.at(i);
				sv_length_pos = str_tmp.find("sv_length=");
				if(sv_length_pos!=-1){
					vec_idx = i;
					str_pos = sv_length_pos + 10;
					break;
				}else{
					SVLEN_pos = str_tmp.find("SVLEN=");
					if(SVLEN_pos!=-1){
						vec_idx = i;
						str_pos = SVLEN_pos + 6;
						break;
					}
				}
			}

			if(vec_idx!=-1 and str_pos!=-1){ // found
				str_tmp = str_vec.at(vec_idx);
				start_pos = str_pos; end_pos = -1;
				for(j=start_pos; j<str_tmp.size(); j++){ // search ';'
					if(str_tmp.at(j)==';'){
						end_pos = j - 1;
						break;
					}
				}
				if(end_pos==-1) end_pos = str_tmp.size() - 1;
				sv_len_str = str_tmp.substr(start_pos, end_pos-start_pos+1);
				if(sv_len_str.find_first_not_of("-0123456789")==string::npos) {
					sv_len = stoi(sv_len_str);
					sv_len_vec.push_back(sv_len);
					flag = true;
					if(sv_len>=MAX_VALID_SVLEN) { sv_len_str = ""; sv_len = 0; flag = false; }
				}else { sv_len_str = ""; flag = false; }
			}

			// 'SVLEN' was not detected, then further detect 'END'
			if(flag==false){
				vec_idx = str_pos = -1;
				for(i=3; i<str_vec.size(); i++){
					str_tmp = str_vec.at(i);
					sv_END_pos = str_tmp.find("END=");
					if(sv_END_pos!=-1){
						vec_idx = i;
						str_pos = sv_END_pos + 4;
						break;
					}
				}

				if(vec_idx!=-1 and str_pos!=-1){ // found
					str_tmp = str_vec.at(vec_idx);
					start_pos = str_pos; end_pos = -1;
					for(j=start_pos; j<str_tmp.size(); j++){ // search ';'
						if(str_tmp.at(j)==';'){
							end_pos = j - 1;
							break;
						}
					}
					if(end_pos==-1) end_pos = str_tmp.size() - 1;
					sv_len_str = str_tmp.substr(start_pos, end_pos-start_pos+1);
					if(sv_len_str.find_first_not_of("-0123456789")==string::npos) { sv_len = stoi(sv_len_str) - stoi(str_vec.at(1)) + 1; sv_len_vec.push_back(sv_len); flag = true; }
				}else{ // not found
					// check sequence
					is_seq_flag_ref = isSeq(str_vec.at(3));
					if(is_seq_flag_ref){
						is_seq_flag_alt = isSeq(str_vec.at(4));
						if(is_seq_flag_alt){
							sv_len = (int32_t)str_vec.at(4).size() - (int32_t)str_vec.at(3).size();
							sv_len_vec.push_back(sv_len);
							flag = true;
						}else{
							str_tmp = str_vec.at(4);
							is_comma_flag = isComma(str_vec.at(4));
							if(is_comma_flag){
								comma_str_vec = split(str_vec.at(4), ",");
								sv_len1 = (int32_t)comma_str_vec.at(0).size() - (int32_t)str_vec.at(3).size();
								sv_len2 = (int32_t)comma_str_vec.at(1).size() - (int32_t)str_vec.at(3).size();
								sv_len_vec.push_back(sv_len1);
								sv_len_vec.push_back(sv_len2);
								flag = true;
							}else{
								if(str_tmp.at(0)=='<' and str_tmp.at(str_tmp.size()-1)=='>'){
									sv_len = 0;
									sv_len_vec.push_back(sv_len);
									flag = true;
								}
							}
						}
					}

					// no sequence
					if(flag==false){
						start_pos = stoi(str_vec.at(1));
						end_pos = stoi(str_vec.at(2));
						sv_len = end_pos - start_pos + 1;
						sv_len_vec.push_back(sv_len);
						flag = true;
					}
				}
			}
			if(flag==false) cout << "missing SVLEN information" << endl;
		}else sv_len_vec.push_back(sv_len);
	}

	return sv_len_vec;
}

bool isSeq(string &seq){
	bool flag = true;
	bool comma_flag = false;
	char ch;
	for(size_t i=0; i<seq.size(); i++){
		ch = seq.at(i);
		switch(ch){
			case 'A':
			case 'a':
			case 'C':
			case 'c':
			case 'G':
			case 'g':
			case 'T':
			case 't':
			case 'N':
			case 'n':
			case 'M':
			case 'm':
			case 'R':
			case 'r':
			case 'S':
			case 's':
			case 'V':
			case 'v':
			case 'W':
			case 'w':
			case 'Y':
			case 'y':
			case 'H':
			case 'h':
			case 'K':
			case 'k':
			case 'D':
			case 'd':
			case 'B':
			case 'b':
				//flag = true;
				break;
			case ',':
				flag = false;
				comma_flag = true;
				break;
			default:
				flag = false;
				break;
			//default: cerr << __func__ << ": unknown base: " << ch << endl; exit(1);
		}
		if(comma_flag) flag = false;
		if(!flag) break;
	}
	//if(comma_flag) flag = false;
	return flag;
}

bool isComma(string &seq){
	bool flag = false;
	char ch;
	for(size_t i=0; i<seq.size(); i++){
		ch = seq.at(i);
		if(ch==','){
			flag = true;
			break;
		}
	}
	return flag;
}

// allocate SV item
SV_item *allocateSVItem(string &chrname, size_t startPos, size_t endPos, string &chrname2, size_t startPos2, size_t endPos2, string &sv_type_str, int32_t sv_len, string &ref_seq, string &alt_seq){
	size_t sv_type;

	if(sv_type_str.compare("UNC")==0){
		sv_type = VAR_UNC;
	}else if(sv_type_str.compare("INS")==0 or sv_type_str.compare("insertion")==0){
		sv_type = VAR_INS;
	}else if(sv_type_str.compare("DEL")==0 or sv_type_str.compare("deletion")==0){
		sv_type = VAR_DEL;
	}else if(sv_type_str.compare("DUP")==0 or sv_type_str.compare("duplication")==0 or sv_type_str.compare("DUP:TANDEM")==0){
		sv_type = VAR_DUP;
	}else if(sv_type_str.compare("INV")==0 or sv_type_str.compare("inversion")==0){
		sv_type = VAR_INV;
	}else if(sv_type_str.compare("TRA")==0 or sv_type_str.compare("translocation")==0){
		sv_type = VAR_TRA;
	}else if(sv_type_str.compare("BND")==0){
		sv_type = VAR_BND;
	}else if(sv_type_str.compare("INVTRA")==0){
		sv_type = VAR_INV_TRA;
	}else if(sv_type_str.compare("MIX")==0){
		sv_type = VAR_MIX;
	}else if(sv_type_str.compare("MNP")==0){
		sv_type = VAR_MNP;
	}else if(sv_type_str.compare("SNV")==0){
		sv_type = VAR_SNV;
	}else if(sv_type_str.compare("CNV")==0 or sv_type_str.compare("cnv")==0){
		sv_type = VAR_CNV;
	}else{
		sv_type = VAR_UNC;
	}

	SV_item *item = new SV_item();
	item->chrname = chrname;
	item->startPos = startPos;
	item->endPos = endPos;
	item->chrname2 = chrname2;
	item->startPos2 = startPos2;
	item->endPos2 = endPos2;
	item->sv_type = sv_type;
	item->sv_len = sv_len;
	item->ref_seq = ref_seq;
	item->alt_seq = alt_seq;
	item->overlapped = false;
	item->validFlag = true;
//	item->seqcons = "-";
	return item;
}

// remove mate SV items
void removemateSVItems(string &mate_filename, vector<SV_item*> &sv_item_vec){
	size_t i, j, mate_num;
	SV_item *item, *item2;
	ofstream redudant_file;
	string line, sv_type_str;

	cout << ">>>>>>>>> Remove duplicated mate variant items <<<<<<<<<" << endl;
	cout << "Before removing duplicated mate variant items, data size: " << sv_item_vec.size() << endl;
	outConvertScreenFile << ">>>>>>>>> Remove duplicated mate variant items <<<<<<<<<" << endl;
	outConvertScreenFile << "Before removing duplicated mate variant items, data size: " << sv_item_vec.size() << endl;

	for(i=0; i<sv_item_vec.size(); i++){
		item = sv_item_vec.at(i);
		if(item->validFlag and item->chrname.size()>0 and item->chrname2.size()>0){
			for(j=i+1; j<sv_item_vec.size(); j++){
				item2 = sv_item_vec.at(j);
				if(item2->validFlag and item2->chrname.size()>0 and item2->chrname2.size()>0)
					if(item->chrname.compare(item2->chrname2)==0 and item->startPos==item2->startPos2 and item->endPos==item2->endPos2
						and item->chrname2.compare(item2->chrname)==0 and item->startPos2==item2->startPos and item->endPos2==item2->endPos)
						item2->validFlag = false;
			}
		}
	}
	//outputPathname
	mate_filename = suboutputDirnamePath + "/" + mate_filename;
	redudant_file.open(mate_filename);
	if(!redudant_file.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << mate_filename << endl;
		exit(1);
	}

	line = "#chr\tstartPos\tendPos\tSVType\tSVLen\tRef\tAlt";
	redudant_file << line << endl;

	// remove items
	mate_num = 0;
	for(i=0; i<sv_item_vec.size(); ){
		item = sv_item_vec.at(i);
		if(item->validFlag==false){
			switch(item->sv_type){
				case VAR_UNC: sv_type_str = "UNC"; break;
				case VAR_INS: sv_type_str = "INS"; break;
				case VAR_DEL: sv_type_str = "DEL"; break;
				case VAR_DUP: sv_type_str = "DUP"; break;
				case VAR_INV: sv_type_str = "INV"; break;
				case VAR_TRA: sv_type_str = "TRA"; break;
				case VAR_BND: sv_type_str = "BND"; break;
				case VAR_INV_TRA: sv_type_str = "INVTRA"; break;
				case VAR_MIX: sv_type_str = "MIX"; break;
				case VAR_SNV: sv_type_str = "SNV"; break;
				default:
					cerr << "line=" << __LINE__ << ", invalid sv type: " << item->sv_type << endl;
					exit(1);
			}
			line = item->chrname + "\t" + to_string(item->startPos) + "\t" + to_string(item->endPos) + "\t" + sv_type_str + "\t" + to_string(item->sv_len) + "\t" + item->ref_seq + "\t" + item->alt_seq;
			redudant_file << line << endl;
			mate_num ++;

			// delete mate items
			sv_item_vec.erase(sv_item_vec.begin()+i);
			delete item;
		}else i++;
	}
	redudant_file.close();

	cout << "After removing duplicated mate variant items, data size: " << sv_item_vec.size() << endl;
	outConvertScreenFile << "After removing duplicated mate variant items, data size: " << sv_item_vec.size() << endl;
	cout << mate_num << " duplicated mate items were saved to " << mate_filename << endl;
	outConvertScreenFile << mate_num << " duplicated mate items were saved to " << mate_filename << endl;
}

//remove SNV
void removeSNVItems(string &snv_filename, vector<SV_item*> &sv_item_vec){
	size_t i, snv_num;
	SV_item *item;
	ofstream snv_file;
	string line, sv_type_str;

	cout << ">>>>>>>>> Remove snv items <<<<<<<<<" << endl;
	cout << "Before removing snv items, data size: " << sv_item_vec.size() << endl;
	outConvertScreenFile << ">>>>>>>>> Remove snv items <<<<<<<<<" << endl;
	outConvertScreenFile << "Before removing snv items, data size: " << sv_item_vec.size() << endl;

	snv_filename = outputPathname + snv_filename;
	snv_file.open(snv_filename);
	if(!snv_file.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << snv_filename << endl;
		exit(1);
	}

	line = "#chr\tstartPos\tendPos\tSVType\tSVLen\tRef\tAlt";
	snv_file << line << endl;

	// remove items
	snv_num = 0;
	for(i=0; i<sv_item_vec.size(); ){
		item = sv_item_vec.at(i);
		if(item->sv_len == 0){
			sv_type_str = "SNV";
			line = item->chrname + "\t" + to_string(item->startPos) + "\t" + to_string(item->endPos) + "\t" + sv_type_str + "\t" + to_string(item->sv_len) + "\t" + item->ref_seq + "\t" + item->alt_seq;

			snv_file << line << endl;
			snv_num ++;

			// delete snv items
			sv_item_vec.erase(sv_item_vec.begin()+i);
			delete item;
		}else i++;
	}
	snv_file.close();

	cout << "After removing snv items, data size: " << sv_item_vec.size() << endl;
	outConvertScreenFile << "After removing snv items, data size: " << sv_item_vec.size() << endl;
	cout << snv_num << " snv items were saved to " << snv_filename << endl;
	outConvertScreenFile << snv_num << " snv items were saved to " << snv_filename << endl;
}

// release sv item vector
void releaseSVItemVec(vector<SV_item*> &sv_item_vec){
	SV_item *item;
	for(size_t i=0; i<sv_item_vec.size(); i++){
		item = sv_item_vec.at(i);
		delete item;
	}
	vector<SV_item*>().swap(sv_item_vec);
}

