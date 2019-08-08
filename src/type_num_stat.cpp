#include "type_num_stat.h"
#include "extvab.h"

void SVTypeNumStat(string &user_file, string &benchmark_file, int32_t max_valid_reg_thres){
	typeNumStatDirname = outputPathname + typeNumStatDirname;
	mkdir(typeNumStatDirname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	if(max_valid_reg_thres>0){
		cout << ">>>>>>>>> Before filtering long SV regions: <<<<<<<<<" << endl;
		outStatScreenFile << ">>>>>>>>> Before filtering long SV regions: <<<<<<<<<" << endl;
	}
	SVTypeNumStatOp(user_file, benchmark_file, size_div_vec, 0, typeNumStatDirname);
	if(max_valid_reg_thres>0){
		cout << "\n>>>>>>>>> After filtering long SV regions: <<<<<<<<<" << endl;
		outStatScreenFile << "\n>>>>>>>>> After filtering long SV regions: <<<<<<<<<" << endl;
		SVTypeNumStatOp(user_file, benchmark_file, size_div_vec, max_valid_reg_thres, typeNumStatDirname);
	}
}

void SVTypeNumStatOp(string &user_file, string &benchmark_file, vector<size_t> &size_div_vec, int32_t max_valid_reg_thres, string &dirname){
	vector<SV_item*> user_data, benchmark_data, long_sv_data;
	vector< vector<SV_item*> > divided_vec1, divided_vec2;

	user_data = loadData(user_file);
	benchmark_data = loadData(benchmark_file);

	if(max_valid_reg_thres>0){
		long_sv_data = getLongSVReg(user_data, max_valid_reg_thres);
		cout << "Total long SV data size: " << long_sv_data.size() << endl;
		outStatScreenFile << "Total long SV data size: " << long_sv_data.size() << endl;
	}

	cout << "Total user data size: " << user_data.size() << endl;
	cout << "Total benchmark data size: " << benchmark_data.size() << endl;

	outStatScreenFile << "Total user data size: " << user_data.size() << endl;
	outStatScreenFile << "Total benchmark data size: " << benchmark_data.size() << endl;

	// prepare user called result
	divided_vec1 = sizeDivideSV(user_data, size_div_vec);

	// prepare benchmark sv result
	divided_vec2 = sizeDivideSV(benchmark_data, size_div_vec);

	// compute Num statistics for each sub-set
	computeTypeNumStat(divided_vec1, divided_vec2, dirname);

	// free memory
	destroyData(user_data);
	destroyData(benchmark_data);
	if(max_valid_reg_thres>0) destroyData(long_sv_data);
	destroySizeDividedData(divided_vec1);
	destroySizeDividedData(divided_vec2);
}

vector< vector<SV_item*> > sizeDivideSV(vector<SV_item*> &sv_data, vector<size_t> &size_div_vec){
	vector< vector<SV_item*> > divided_vec;
	vector<SV_item*> sub_vec[size_div_vec.size()+2];
	SV_item *item, *item_tmp;
	size_t i, j, div_pos;
	int32_t dist, idx;

	for(i=0; i<sv_data.size(); i++){
		item = sv_data.at(i);
		idx = -1;
		if(item->sv_type!=VAR_TRA and item->sv_type!=VAR_BND){ // INS, DEL, DUP, INV
			dist = item->endPos - item->startPos + 1;
			for(j=0; j<size_div_vec.size(); j++){
				div_pos = size_div_vec.at(j);
				if(dist<(int32_t)div_pos){ idx = j; break; }
			}
			if(idx==-1) idx = size_div_vec.size();
		}else // TRA or BND
			idx = size_div_vec.size() + 1;

		item_tmp = itemdup(item);
		sub_vec[idx].push_back(item_tmp);
	}

	for(i=0; i<size_div_vec.size()+2; i++) divided_vec.push_back(sub_vec[i]);

	return divided_vec;
}

void computeTypeNumStat(vector< vector<SV_item*> > &divided_vec1, vector< vector<SV_item*> > &divided_vec2, string &dirname){
	vector<SV_item*> sv_vec1, sv_vec2;
	string file_prefix, file_prefix_tmp;
	size_t start_pos, end_pos;

	file_prefix = dirname + "type_num_stat";

	start_pos = 1;
	for(size_t i=0; i<divided_vec1.size(); i++){
		sv_vec1 = divided_vec1.at(i);
		sv_vec2 = divided_vec2.at(i);

		if(i<size_div_vec.size()){
			end_pos = size_div_vec.at(i);
			cout << "\n>>>>>>>>> SV size: " << start_pos << "-" << end_pos << ":" << " <<<<<<<<<" << endl;
			outStatScreenFile << "\n>>>>>>>>> SV size: " << start_pos << "-" << end_pos << ":" << " <<<<<<<<<" << endl;
			file_prefix_tmp = file_prefix + "_" + to_string(start_pos) + "_" + to_string(end_pos);
			start_pos = end_pos + 1;
		}else if(i==size_div_vec.size()){
			cout << "\n>>>>>>>>> SV size: " << ">=" << start_pos << ":" << " <<<<<<<<<" << endl;
			outStatScreenFile << "\n>>>>>>>>> SV size: " << ">=" << start_pos << ":" << " <<<<<<<<<" << endl;
			file_prefix_tmp = file_prefix + "_" + to_string(start_pos) + "_larger";
		}else if(i==size_div_vec.size()+1){
			cout << "\n>>>>>>>>> TRA breakpoints:" << " <<<<<<<<<" << endl;
			outStatScreenFile << "\n>>>>>>>>> TRA breakpoints:" << " <<<<<<<<<" << endl;
			file_prefix_tmp = file_prefix + "_TRA_BND_BP";
		}
		if(i<size_div_vec.size()+1) computeNumStat(sv_vec1, sv_vec2, file_prefix_tmp);
		else computeNumStatTra(sv_vec1, sv_vec2, file_prefix_tmp);
	}
}

void destroySizeDividedData(vector< vector<SV_item*> > &divided_vec){
	destroyResultData(divided_vec);
}

