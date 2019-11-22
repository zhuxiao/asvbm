#include "convert.h"
#include "util.h"

// convert data
void convertBed(const string &infilename, const string &outfilename, bool removeDuplicatedItemFlag){
	ifstream infile;
	string line, chrname, chrname2, sv_type_str;
	vector<string> str_vec;
	size_t start_pos, endpos, start_pos2, endpos2;
	int32_t sv_len;
	vector<SV_item*> sv_item_vec;
	SV_item *sv_item;

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
				sv_type_str = getSVType(str_vec);
				sv_len = getSVLen(str_vec, sv_type_str);

				if(sv_type_str.compare("DEL")==0 and sv_len>0) sv_len = -sv_len;
				if(sv_type_str.compare("TRA")==0 or sv_type_str.compare("BND")==0){
					chrname2 = str_vec.at(3);
					start_pos2 = stoi(str_vec.at(4));
					endpos2 = stoi(str_vec.at(5));
				}else{
					chrname2 = "";
					start_pos2 = endpos2 = 0;
				}

				sv_item = allocateSVItem(chrname, start_pos, endpos, chrname2, start_pos2, endpos2, sv_type_str, sv_len);
				sv_item_vec.push_back(sv_item);
			}
		}
	}
	infile.close();

	// remove duplicated sv items
	if(removeDuplicatedItemFlag) removeDuplicatedSVItems(sv_item_vec);

	// save to file
	output2File(outfilename, sv_item_vec, outConvertScreenFile);

	// release memory
	releaseSVItemVec(sv_item_vec);
}

// convert vcf result
void convertVcf(const string &infilename, const string &outfilename, bool removeDuplicatedItemFlag){
	ifstream infile;
	string line, chrname, start_pos_str, endpos_str, chrname2, start_pos_str2, endpos_str2, sv_type_str, sv_len_str, bnd_str, bnd_pos_str;
	size_t start_pos, endpos, start_pos2, endpos2;
	size_t i, end_pos_ref;
	int32_t sv_len;
	vector<string> str_vec, info_vec, sub_info_vec, bnd_pos_vec;
	bool is_seq_flag_ref;
	vector<SV_item*> sv_item_vec;
	SV_item *sv_item;

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
				start_pos_str = str_vec.at(1);

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
					if(sub_info_vec.at(0).compare("SVLEN")==0)
						sv_len_str = sub_info_vec.at(1);
				}

				// get sv type and length
				if(sv_type_str.size()==0) sv_type_str = getSVType(str_vec);
				if(sv_len_str.size()==0) { sv_len = getSVLen(str_vec, sv_type_str); sv_len_str = to_string(sv_len); }
				if(endpos_str.size()==0) {
					is_seq_flag_ref = isSeq(str_vec.at(3));
					if(is_seq_flag_ref){
						end_pos_ref = stoi(start_pos_str) + str_vec.at(3).size() - 1;
						endpos_str = to_string(end_pos_ref);
					}
				}

				if(sv_type_str.size()>0){
					if(sv_type_str.compare("TRA")==0 or sv_type_str.compare("BND")==0){ // TRA or BND
						if(sv_type_str.compare("BND")==0){ // BND
							bnd_str = str_vec.at(4);
							if(bnd_str.at(0)==']' or bnd_str.at(0)=='[')
								bnd_pos_str = bnd_str.substr(1, bnd_str.size()-3);
							else if(bnd_str.at(0)=='N')
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

					start_pos = stoi(start_pos_str);
					endpos = stoi(endpos_str);
					start_pos2 = stoi(start_pos_str2);
					endpos2 = stoi(endpos_str2);
					sv_len = stoi(sv_len_str);

					sv_item = allocateSVItem(chrname, start_pos, endpos, chrname2, start_pos2, endpos2, sv_type_str, sv_len);
					sv_item_vec.push_back(sv_item);

				}else{
					cout << line << endl;
					cout << "missing SVTYPE information" << endl;
				}
			}
		}
	}
	infile.close();

	// remove duplicated sv items
	if(removeDuplicatedItemFlag) removeDuplicatedSVItems(sv_item_vec);

	// save to file
	output2File(outfilename, sv_item_vec, outConvertScreenFile);

	// release memory
	releaseSVItemVec(sv_item_vec);
}

// convert data
void convertCsv(const string &infilename, const string &outfilename, bool removeDuplicatedItemFlag){
	ifstream infile;
	string line, chrname, chrname2, sv_type_str;
	vector<string> str_vec;
	size_t start_pos, endpos, start_pos2, endpos2;
	int32_t sv_len;
	vector<SV_item*> sv_item_vec;
	SV_item *sv_item;

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
					start_pos = stoi(str_vec.at(1));
					endpos = stoi(str_vec.at(2));

					// get sv type and length
					sv_type_str = getSVType(str_vec);
					sv_len = getSVLen(str_vec, sv_type_str);
					if(sv_type_str.compare("DEL")==0 and sv_len>0) sv_len = -sv_len;
					if(sv_type_str.compare("TRA")==0 or sv_type_str.compare("BND")==0){
						chrname2 = str_vec.at(3);
						start_pos2 = stoi(str_vec.at(4));
						endpos2 = stoi(str_vec.at(5));
					}else{
						chrname2 = "";
						start_pos2 = endpos2 = 0;
					}

					sv_item = allocateSVItem(chrname, start_pos, endpos, chrname2, start_pos2, endpos2, sv_type_str, sv_len);
					sv_item_vec.push_back(sv_item);
				}
			}
		}
	}
	infile.close();

	// remove duplicated sv items
	if(removeDuplicatedItemFlag) removeDuplicatedSVItems(sv_item_vec);

	// save to file
	output2File(outfilename, sv_item_vec, outConvertScreenFile);

	// release memory
	releaseSVItemVec(sv_item_vec);
}

// convert data
void convertNm(const string &infilename, const string &outfilename, bool removeDuplicatedItemFlag){
	ifstream infile;
	string line, chrname, chrname2, sv_type_str;
	vector<string> str_vec;
	int32_t sv_len, ref_seq_len, start_pos, endpos, start_pos2, endpos2;
	vector<SV_item*> sv_item_vec;
	SV_item *sv_item;

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

				// get sv type
				sv_type_str = getSVType(str_vec);

				// get sv length
				sv_len = getSVLen(str_vec, sv_type_str);

				chrname2 = "";
				start_pos2 = endpos2 = 0;

				sv_item = allocateSVItem(chrname, start_pos, endpos, chrname2, start_pos2, endpos2, sv_type_str, sv_len);
				sv_item_vec.push_back(sv_item);
			}
		}
	}
	infile.close();

	// remove duplicated sv items
	if(removeDuplicatedItemFlag) removeDuplicatedSVItems(sv_item_vec);

	// save to file
	output2File(outfilename, sv_item_vec, outConvertScreenFile);

	// release memory
	releaseSVItemVec(sv_item_vec);
}

// get sv type
string getSVType(vector<string> &str_vec){
	string sv_type_str, str_tmp;
	size_t i, j;
	vector<string> info_vec;
	int32_t sv_type_pos, SVTYPE_pos, vec_idx, str_pos, start_pos, end_pos;

	sv_type_str = "";
	for(i=3; i<str_vec.size(); i++){
		str_tmp = str_vec.at(i);
		if(str_tmp.compare("INS")==0 or str_tmp.compare("DEL")==0 or str_tmp.compare("DUP")==0 or str_tmp.compare("INV")==0 or str_tmp.compare("TRA")==0 or str_tmp.compare("BND")==0 or str_tmp.compare("INVTRA")==0 or str_tmp.compare("MIX")==0 or str_tmp.compare("MNP")==0){
			sv_type_str = str_tmp;
		}else if(str_tmp.compare("<INS>")==0 or str_tmp.compare("<DEL>")==0 or str_tmp.compare("<DUP>")==0 or str_tmp.compare("<INV>")==0 or str_tmp.compare("<TRA>")==0 or str_tmp.compare("<BND>")==0 or str_tmp.compare("<INVTRA>")==0 or str_tmp.compare("<MIX>")==0 or str_tmp.compare("<MNP>")==0){
			sv_type_str = str_tmp.substr(1, str_tmp.size()-2);
		}else if(str_tmp.compare("insertion")==0 or str_tmp.compare("deletion")==0 or str_tmp.compare("duplication")==0 or str_tmp.compare("inversion")==0 or str_tmp.compare("translocation")==0){
			sv_type_str = str_tmp;
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
		}
	}

	if(sv_type_str.size()==0) sv_type_str = "UNC";

	return sv_type_str;
}

// get SV length
int32_t getSVLen(vector<string> &str_vec, string &sv_type){
	int32_t sv_len, vec_idx, str_pos, sv_length_pos, SVLEN_pos, start_pos, end_pos, sv_END_pos;
	size_t i, j;
	string str_tmp, sv_len_str;
	bool flag, is_seq_flag_ref, is_seq_flag_alt;

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
			if(sv_len_str.find_first_not_of("0123456789")==string::npos) { sv_len = stoi(sv_len_str); flag = true; }
		}
//		else{ // not found
//			start_pos = stoi(str_vec.at(1));
//			end_pos = stoi(str_vec.at(2));
//			sv_len = end_pos - start_pos + 1;
//			flag = true;
//			cout << "line=" << __LINE__ << ", please check." << endl;
//			exit(1);
//		}

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
				if(sv_len_str.find_first_not_of("0123456789")==string::npos) { sv_len = stoi(sv_len_str); flag = true; }
			}else{ // not found
				// check sequence
				is_seq_flag_ref = isSeq(str_vec.at(3));
				if(is_seq_flag_ref){
					is_seq_flag_alt = isSeq(str_vec.at(4));
					if(is_seq_flag_alt){
						sv_len = (int32_t)str_vec.at(4).size() - (int32_t)str_vec.at(3).size();
						flag = true;
					}else{
						str_tmp = str_vec.at(4);
						if(str_tmp.at(0)=='<' and str_tmp.at(str_tmp.size()-1)=='>'){
							sv_len = 0;
							flag = true;
						}
					}
				}

				// no sequence
				if(flag==false){
					start_pos = stoi(str_vec.at(1));
					end_pos = stoi(str_vec.at(2));
					sv_len = end_pos - start_pos + 1;
					flag = true;
				}
			}
		}
		if(flag==false) cout << "missing SVLEN information" << endl;
	}

	return sv_len;
}

bool isSeq(string &seq){
	bool flag = true;
	char ch;
	for(size_t i=0; i<seq.size(); i++){
		ch = seq.at(i);
		if(ch!='A' and ch!='C' and ch!='G' and ch!='T' and ch!='a' and ch!='c' and ch!='g' and ch!='t'){
			flag = false;
			break;
		}
	}

	return flag;
}

// allocate SV item
SV_item *allocateSVItem(string &chrname, size_t startPos, size_t endPos, string &chrname2, size_t startPos2, size_t endPos2, string &sv_type_str, int32_t sv_len){
	size_t sv_type;

	if(sv_type_str.compare("UNC")==0){
		sv_type = VAR_UNC;
	}else if(sv_type_str.compare("INS")==0 or sv_type_str.compare("insertion")==0){
		sv_type = VAR_INS;
	}else if(sv_type_str.compare("DEL")==0 or sv_type_str.compare("deletion")==0){
		sv_type = VAR_DEL;
	}else if(sv_type_str.compare("DUP")==0 or sv_type_str.compare("duplication")==0){
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
	item->overlapped = false;
	item->validFlag = true;
	return item;
}

// remove duplicate SV items
void removeDuplicatedSVItems(vector<SV_item*> &sv_item_vec){
	size_t i, j;
	SV_item *item, *item2;

	cout << ">>>>>>>>> Remove duplicated variant items <<<<<<<<<" << endl;
	cout << "Before removing duplicated variant items, data size: " << sv_item_vec.size() << endl;
	outConvertScreenFile << ">>>>>>>>> Remove duplicated variant items <<<<<<<<<" << endl;
	outConvertScreenFile << "Before removing duplicated variant items, data size: " << sv_item_vec.size() << endl;

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

	// remove items
	for(i=0; i<sv_item_vec.size(); ){
		item = sv_item_vec.at(i);
		if(item->validFlag==false){
			sv_item_vec.erase(sv_item_vec.begin()+i);
			delete item;
		}else i++;
	}

	cout << "After removing duplicated variant items, data size: " << sv_item_vec.size() << endl;
	outConvertScreenFile << "After removing duplicated variant items, data size: " << sv_item_vec.size() << endl;
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

