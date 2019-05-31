#include "size_dif_stat.h"
#include "extvab.h"
#include "util.h"


void SVSizeDifStat(string &user_file, string &benchmark_file, int32_t max_valid_reg_thres){
	sizeDifStatDirname = outputPathname + sizeDifStatDirname;
	mkdir(sizeDifStatDirname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	if(max_valid_reg_thres>0){
		cout << ">>>>>>>>> Before filtering long SV regions: <<<<<<<<<" << endl;
		outStatScreenFile << ">>>>>>>>> Before filtering long SV regions: <<<<<<<<<" << endl;
	}
	SVSizeDifStatOp(user_file, benchmark_file, 0, sizeDifStatDirname);
	if(max_valid_reg_thres>0){
		cout << "\n>>>>>>>>> After filtering long SV regions: <<<<<<<<<" << endl;
		outStatScreenFile << "\n>>>>>>>>> After filtering long SV regions: <<<<<<<<<" << endl;
		SVSizeDifStatOp(user_file, benchmark_file, max_valid_reg_thres, sizeDifStatDirname);
	}
}

void SVSizeDifStatOp(string &user_file, string &benchmark_file, int32_t max_valid_reg_thres, string &dirname){
	vector<SV_item*> user_data, benchmark_data, long_sv_data;
	vector<SV_pair*> sv_pair_vec;
	vector< vector<int32_t> > dif_stat_vec;
	vector<size_t> ratio_stat_vec;
	string svSizeDifRatioFilename_tmp, svSizeDifStatFilename_tmp, svSizeRatioStatFilename_tmp;

	// load data
	user_data = loadData(user_file);
	benchmark_data = loadData(benchmark_file);

	if(max_valid_reg_thres>0)
		long_sv_data = getLongSVReg(user_data, max_valid_reg_thres);
	cout << "user data size: " << user_data.size() << endl;
	cout << "benchmark data size: " << benchmark_data.size() << endl;

	outStatScreenFile << "user data size: " << user_data.size() << endl;
	outStatScreenFile << "benchmark data size: " << benchmark_data.size() << endl;

	// compute overlapped SV pairs
	sv_pair_vec = computeOverlapSVPair(user_data, benchmark_data);

	// output pair statistics data to file
	svSizeDifRatioFilename_tmp = dirname + svSizeDifRatioFilename;
	if(max_valid_reg_thres>0) svSizeDifRatioFilename_tmp += "_long_filtered";
	outputPairDataToFile(svSizeDifRatioFilename_tmp, sv_pair_vec);

	// compute size difference statistics
	dif_stat_vec = computeDifStatVec(sv_pair_vec);

	// save to file
	svSizeDifStatFilename_tmp = dirname + svSizeDifStatFilename;
	if(max_valid_reg_thres>0) svSizeDifStatFilename_tmp += "_long_filtered";
	outputDifStatToFile(svSizeDifStatFilename_tmp, dif_stat_vec);

	// compute size ratio statistics
	ratio_stat_vec = computeRatioStatVec(sv_pair_vec, ratio_div_vec);

	// save to file
	svSizeRatioStatFilename_tmp = dirname + svSizeRatioStatFilename;
	if(max_valid_reg_thres>0) svSizeRatioStatFilename_tmp += "_long_filtered";
	outputRatioStatToFile(svSizeRatioStatFilename_tmp, ratio_stat_vec, ratio_div_vec);

	// release memory
	destroyPairData(sv_pair_vec);
	destroyData(user_data);
	destroyData(benchmark_data);
	if(max_valid_reg_thres>0) destroyData(long_sv_data);
}

vector<SV_pair*> computeOverlapSVPair(vector<SV_item*> &data1, vector<SV_item*> &data2){
	vector<SV_pair*> sv_pair_vec;
	SV_item *item1, *item2;
	SV_pair *pair_item;
	vector<size_t> overlap_type_vec;
	size_t i, j;
	int32_t mid_loc1, mid_loc2, reg_size1, reg_size2, dif_tmp;
	double dif_size, size_ratio;

	for(i=0; i<data1.size(); i++){
		item1 = data1.at(i);
		for(j=0; j<data2.size(); j++){
			item2 = data2.at(j);
			overlap_type_vec = computeOverlapType(item1, item2);
			if(overlap_type_vec.size()>0 and overlap_type_vec.at(0)!=NO_OVERLAP){
				if(item1->sv_type!=VAR_TRA and item2->sv_type!=VAR_TRA){
					pair_item = new SV_pair();
					pair_item->sv_item1 = item1;
					pair_item->sv_item2 = item2;
					sv_pair_vec.push_back(pair_item);
					break;
				}
			}
		}
	}

	for(i=0; i<sv_pair_vec.size(); i++){
		pair_item = sv_pair_vec.at(i);
		item1 = pair_item->sv_item1;
		item2 = pair_item->sv_item2;
		mid_loc1 = (double)(item1->startPos + item1->endPos) / 2;
		mid_loc2 = (double)(item2->startPos + item2->endPos) / 2;
		reg_size1 = item1->endPos - item1->startPos + 1;
		reg_size2 = item2->endPos - item2->startPos + 1;

		dif_size = mid_loc1 - mid_loc2;
		dif_tmp = reg_size1 - reg_size2;
		if(dif_tmp<0) dif_tmp = -dif_tmp;
		if(dif_tmp<MIN_DIF_SIZE) size_ratio = 1;
		else size_ratio = (double)reg_size1 / reg_size2;
		pair_item->dif_size = dif_size;
		pair_item->size_ratio = size_ratio;
	}

	return sv_pair_vec;
}

// output pair data to file
void outputPairDataToFile(string &filename, vector<SV_pair*> &sv_pair_vec){
	ofstream outfile;
	string line, sv_type_str1, sv_type_str2, balancedTraFlag1_str, balancedTraFlag2_str;
	SV_pair *pair_item;
	SV_item *item1, *item2;


	outfile.open(filename);
	if(!outfile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << filename << endl;
		exit(1);
	}

	for(size_t i=0; i<sv_pair_vec.size(); i++){
		pair_item = sv_pair_vec.at(i);
		item1 = pair_item->sv_item1;
		item2 = pair_item->sv_item2;

		switch(item1->sv_type){
			case VAR_UNC: sv_type_str1 = "UNC"; break;
			case VAR_INS: sv_type_str1 = "INS"; break;
			case VAR_DEL: sv_type_str1 = "DEL"; break;
			case VAR_DUP: sv_type_str1 = "DUP"; break;
			case VAR_INV: sv_type_str1 = "INV"; break;
			case VAR_TRA: sv_type_str1 = "TRA"; break;
			case VAR_BND: sv_type_str1 = "BND"; break;
			case VAR_INV_TRA: sv_type_str1 = "INVTRA"; break;
			case VAR_MIX: sv_type_str1 = "MIX"; break;
			default: cerr << "line=" << __LINE__ << ", invalid SV type: " << item1->sv_type << endl; exit(1);
		}
		switch(item2->sv_type){
			case VAR_UNC: sv_type_str2 = "UNC"; break;
			case VAR_INS: sv_type_str2 = "INS"; break;
			case VAR_DEL: sv_type_str2 = "DEL"; break;
			case VAR_DUP: sv_type_str2 = "DUP"; break;
			case VAR_INV: sv_type_str2 = "INV"; break;
			case VAR_TRA: sv_type_str2 = "TRA"; break;
			case VAR_BND: sv_type_str2 = "BND"; break;
			case VAR_INV_TRA: sv_type_str2 = "INVTRA"; break;
			case VAR_MIX: sv_type_str2 = "MIX"; break;
			default: cerr << "line=" << __LINE__ << ", invalid SV type: " << item2->sv_type << endl; exit(1);
		}

		if(item1->sv_type==VAR_TRA or item1->sv_type==VAR_BND){
			if(item1->balancedTraFlag) balancedTraFlag1_str = BALANCED_TRA_STR;
			else balancedTraFlag1_str = UNBALANCED_TRA_STR;
		}else balancedTraFlag1_str = "-";
		if(item2->sv_type==VAR_TRA or item2->sv_type==VAR_BND){
			if(item2->balancedTraFlag) balancedTraFlag2_str = BALANCED_TRA_STR;
			else balancedTraFlag2_str = UNBALANCED_TRA_STR;
		}else balancedTraFlag2_str = "-";

		line = to_string(pair_item->dif_size) + "\t" + to_string(pair_item->size_ratio);
		line += "\t" + item1->chrname + "\t" + to_string(item1->startPos) + "\t" + to_string(item1->endPos) + "\t" + sv_type_str1 + "\t" + to_string(item1->sv_len) + "\t" + balancedTraFlag1_str;
		line += "\t" + item2->chrname + "\t" + to_string(item2->startPos) + "\t" + to_string(item2->endPos) + "\t" + sv_type_str2 + "\t" + to_string(item2->sv_len) + "\t" + balancedTraFlag2_str;
		outfile << line << endl;
	}

	outfile.close();
}

vector< vector<int32_t> > computeDifStatVec(vector<SV_pair*> &sv_pair_vec){
	vector< vector<int32_t> > dif_stat_vec;
	vector<int32_t> positive_dif_stat_vec, negative_dif_stat_vec;
	size_t i, positive_vec_size, negative_vec_size;
	int32_t minValue, maxValue, tmp, idx;
	SV_pair *item;

	if(sv_pair_vec.size()>0){
		// get the minimal and the maximal difference size
		minValue = INT_MAX; maxValue = INT_MIN;
		for(i=0; i<sv_pair_vec.size(); i++){
			item = sv_pair_vec.at(i);
			if(minValue>item->dif_size) minValue = item->dif_size;
			if(maxValue<item->dif_size) maxValue = item->dif_size;
		}

		if(minValue<0){
			tmp = -minValue;
			if(tmp<=MAX_DIF_SIZE) negative_vec_size = (tmp - 1) / SUB_DIF_REG_SIZE + 1;
			else negative_vec_size = (MAX_DIF_SIZE - 1) / SUB_DIF_REG_SIZE + 2;
		}else negative_vec_size = 0;
		if(maxValue>=0){
			tmp = maxValue;
			if(tmp<=MAX_DIF_SIZE) positive_vec_size = (tmp - 1) / SUB_DIF_REG_SIZE + 1;
			else positive_vec_size = (MAX_DIF_SIZE - 1) / SUB_DIF_REG_SIZE + 2;
		}else positive_vec_size = 0;

		// push back items
		for(i=0; i<negative_vec_size; i++) negative_dif_stat_vec.push_back(0);
		for(i=0; i<positive_vec_size; i++) positive_dif_stat_vec.push_back(0);

		// update vector
		for(i=0; i<sv_pair_vec.size(); i++){
			item = sv_pair_vec.at(i);
			if(item->dif_size<0){
				tmp = -item->dif_size;
				if(tmp<=MAX_DIF_SIZE) idx = (tmp - 1) / SUB_DIF_REG_SIZE;
				else idx = negative_dif_stat_vec.size() - 1;
				negative_dif_stat_vec.at(idx) ++;
			}else{
				tmp = item->dif_size;
				if(tmp==0) idx = 0;
				else if(tmp<=MAX_DIF_SIZE) idx = (tmp - 1) / SUB_DIF_REG_SIZE;
				else idx = positive_dif_stat_vec.size() - 1;
				positive_dif_stat_vec.at(idx) ++;
			}
		}
	}
	dif_stat_vec.push_back(negative_dif_stat_vec);
	dif_stat_vec.push_back(positive_dif_stat_vec);

	return dif_stat_vec;
}

void outputDifStatToFile(string &svSizeDifStatFilename, vector< vector<int32_t> > &dif_stat_vec){
	vector<int32_t> positive_dif_stat_vec, negative_dif_stat_vec;
	ofstream out_file;
	int32_t i, start_dif, end_dif, num, tmp, max_vec_size, start_idx;
	string line;

	out_file.open(svSizeDifStatFilename);
	if(!out_file.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << svSizeDifStatFilename << endl;
		exit(1);
	}

	negative_dif_stat_vec = dif_stat_vec.at(0);
	positive_dif_stat_vec = dif_stat_vec.at(1);

	max_vec_size = (MAX_DIF_SIZE - 1) / SUB_DIF_REG_SIZE + 2;
	num = 0;
	if(negative_dif_stat_vec.size()>0){
		if(max_vec_size==(int32_t)negative_dif_stat_vec.size()){
			tmp = negative_dif_stat_vec.at(negative_dif_stat_vec.size()-1);
			end_dif = -MAX_DIF_SIZE - 1;
			line = "<" + to_string(end_dif) + "\t" + to_string(tmp);
			out_file << line << endl;
			num ++;
			start_idx = negative_dif_stat_vec.size() - 2;
		}else start_idx = negative_dif_stat_vec.size() - 1;

		for(i=start_idx; i>=0; i--){
			tmp = negative_dif_stat_vec.at(i);
			start_dif = -SUB_DIF_REG_SIZE * (i + 1);
			end_dif = -SUB_DIF_REG_SIZE * i - 1;
			line = to_string(start_dif) + "-" + to_string(end_dif) + "\t" + to_string(tmp);
			out_file << line << endl;
			num ++;
		}
	}
	if(positive_dif_stat_vec.size()>0){
		for(i=0; i<(int32_t)positive_dif_stat_vec.size()-1; i++){
			tmp = positive_dif_stat_vec.at(i);
			if(i==0) start_dif = 0;
			else start_dif = SUB_DIF_REG_SIZE * i + 1;
			end_dif = SUB_DIF_REG_SIZE * (i + 1);
			line = to_string(start_dif) + "-" + to_string(end_dif) + "\t" + to_string(tmp);
			out_file << line << endl;
			num ++;
		}
		if(max_vec_size==(int32_t)positive_dif_stat_vec.size()){
			tmp = positive_dif_stat_vec.at(positive_dif_stat_vec.size()-1);
			start_dif = SUB_DIF_REG_SIZE * i + 1;
			line = ">" + to_string(start_dif) + "\t" + to_string(tmp);
			out_file << line << endl;
			num ++;
		}
	}

	out_file.close();

	cout << num << " size difference items have been saved to file: " << svSizeDifStatFilename << endl;
	outStatScreenFile << num << " size difference items have been saved to file: " << svSizeDifStatFilename << endl;
}

vector<size_t> computeRatioStatVec(vector<SV_pair*> &sv_pair_vec, vector<double> &ratio_div_vec){
	vector<size_t> ratio_stat_vec;
	size_t i, j;
	SV_pair *item;
	int32_t idx;

	// initialize
	for(i=0; i<ratio_div_vec.size()+1; i++) ratio_stat_vec.push_back(0);

	// fill data
	for(i=0; i<sv_pair_vec.size(); i++){
		item = sv_pair_vec.at(i);
		idx = -1;
		for(j=0; j<ratio_div_vec.size(); j++){
			if(item->size_ratio<ratio_div_vec.at(j)){
				idx = j;
				break;
			}
		}
		if(idx==-1) idx = ratio_stat_vec.size() - 1;

		if(idx!=-1) ratio_stat_vec.at(idx) ++;
		else{
			cerr << "line=" << __LINE__ << ", invalid idx=" << idx << ", error!" << endl;
			exit(1);
		}
	}

	return ratio_stat_vec;
}

void outputRatioStatToFile(string &svSizeRatioStatFilename, vector<size_t> &ratio_stat_vec, vector<double> &ratio_div_vec){
	ofstream out_file;
	size_t i, num, ratio_num;
	string line, start_ratio_str, end_ratio_str;
	double start_ratio, end_ratio;
	int32_t fixed_precision;

	out_file.open(svSizeRatioStatFilename);
	if(!out_file.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << svSizeRatioStatFilename << endl;
		exit(1);
	}

	fixed_precision = 1;
	num = 0;
	for(i=0; i<ratio_stat_vec.size(); i++){
		ratio_num = ratio_stat_vec.at(i);
		if(i==0){
			start_ratio = 0;
			end_ratio = ratio_div_vec.at(i);
			start_ratio_str = double2Str(start_ratio, fixed_precision);
			end_ratio_str = double2Str(end_ratio, fixed_precision);
			line = start_ratio_str + "-" + end_ratio_str + "\t" + to_string(ratio_num);
		}else if(i<ratio_stat_vec.size()-1){
			start_ratio = ratio_div_vec.at(i-1);
			end_ratio = ratio_div_vec.at(i);
			start_ratio_str = double2Str(start_ratio, fixed_precision);
			end_ratio_str = double2Str(end_ratio, fixed_precision);
			line = start_ratio_str + "-" + end_ratio_str + "\t" + to_string(ratio_num);
		}else{
			start_ratio = ratio_div_vec.at(i-1);
			start_ratio_str = double2Str(start_ratio, fixed_precision);
			line = ">" + start_ratio_str + "\t" + to_string(ratio_num);
		}
		out_file << line << endl;
		num ++;
	}

	out_file.close();

	cout << num << " size ratio items have been saved to file: " << svSizeRatioStatFilename << endl;
	outStatScreenFile << num << " size ratio items have been saved to file: " << svSizeRatioStatFilename << endl;
}

void destroyPairData(vector<SV_pair*> &sv_vec){
	vector<SV_pair*>::iterator it;
	for(it=sv_vec.begin(); it!=sv_vec.end(); it++)
		delete *it;
	vector<SV_pair*>().swap(sv_vec);
}


