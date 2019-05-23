#include "type_num_stat.h"
#include "extvab.h"

void SVTypeNumStat(string &sv_file1, string &sv_file2, int32_t max_valid_reg_thres){
	if(max_valid_reg_thres>0)
		cout << ">>>>>>>>> Before filtering long SV regions: <<<<<<<<<" << endl;
	SVTypeNumStatOp(sv_file1, sv_file2, size_div_vec, 0);
	if(max_valid_reg_thres>0){
		cout << "\n>>>>>>>>> After filtering long SV regions: <<<<<<<<<" << endl;
		SVTypeNumStatOp(sv_file1, sv_file2, size_div_vec, max_valid_reg_thres);
	}
}

void SVTypeNumStatOp(string &sv_file1, string &sv_file2, vector<size_t> &size_div_vec, int32_t max_valid_reg_thres){
	vector<SV_item*> sv_data1, sv_data2, long_sv_data;
	vector< vector<SV_item*> > divided_vec1, divided_vec2;

	sv_data1 = loadData(sv_file1);
	sv_data2 = loadData(sv_file2);

	if(max_valid_reg_thres>0){
		long_sv_data = getLongSVReg(sv_data1, max_valid_reg_thres);
		cout << "long_sv_data.size: " << long_sv_data.size() << endl;
	}

	cout << "data1.size: " << sv_data1.size() << endl;
	cout << "data2.size: " << sv_data2.size() << endl;

	// prepare user called result
	divided_vec1 = sizeDivideSV(sv_data1, size_div_vec);

	// prepare benchmark sv result
	divided_vec2 = sizeDivideSV(sv_data2, size_div_vec);

	// compute Num statistics for each sub-set
	computeTypeNumStat(divided_vec1, divided_vec2);

	// free memory
	destroyData(sv_data1);
	destroyData(sv_data2);
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

void computeTypeNumStat(vector< vector<SV_item*> > &divided_vec1, vector< vector<SV_item*> > &divided_vec2){
	vector<SV_item*> sv_vec1, sv_vec2;
	string type_num_dirname, file_prefix, file_prefix_tmp;
	size_t start_pos, end_pos;

	type_num_dirname = outputPathname + "type_num_stat";
	mkdir(type_num_dirname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	file_prefix = type_num_dirname + "/type_num_stat";

	start_pos = 1;
	for(size_t i=0; i<divided_vec1.size(); i++){
		sv_vec1 = divided_vec1.at(i);
		sv_vec2 = divided_vec2.at(i);

		if(i<size_div_vec.size()){
			end_pos = size_div_vec.at(i);
			cout << "><><><><><>< SV size: " << start_pos << "-" << end_pos << ":" << endl;
			file_prefix_tmp = file_prefix + "_" + to_string(start_pos) + "_" + to_string(end_pos);
			start_pos = end_pos + 1;
		}else if(i==size_div_vec.size()){
			cout << "><><><><><>< SV size: " << ">=" << start_pos << ":" << endl;
			file_prefix_tmp = file_prefix + "_" + to_string(start_pos) + "_larger";
		}else if(i==size_div_vec.size()+1){
			cout << "><><><><><>< TRA:" << endl;
			file_prefix_tmp = file_prefix + "_TRA_BND";
		}
		computeNumStat(sv_vec1, sv_vec2, file_prefix_tmp);
	}
}

void destroySizeDividedData(vector< vector<SV_item*> > &divided_vec){
	destroyResultData(divided_vec);
}

