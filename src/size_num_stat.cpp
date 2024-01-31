#include "extvab.h"
#include "size_num_stat.h"
#include "gnuplotcall.h"

void FolderInit(){
	data.clear();
	data1.clear();
	data_4.clear();
	data1_4.clear();
	MeticsValues_4.clear();
	MeticsValues1_4.clear();
	allmetric.clear();
	centerDistance.clear();
	sizeratio.clear();
	outputInsideToolDirname = "";
	refRegSizeStatDirname = "1_ref_reg_size_stat/";
	numStatDirname = "2_num_stat/";
	sizeDifStatDirname = "3_size_dif_stat/";
	sizeNumStatDirname = "4_size_num_stat/";
	statScreenFilename = "stat_screen";
}

void SVSizeNumStat(string &user_file, string &benchmark_file, string &ref_file, int32_t max_valid_reg_thres, vector<string> &sv_files1){
//	sizeNumStatDirname = outputPathname + sizeNumStatDirname;
	if(sv_files1.size()>1) sizeNumStatDirname = outputInsideToolDirname + '/' + sizeNumStatDirname;
	else sizeNumStatDirname = outputPathname + sizeNumStatDirname;
//	sizeNumStatDirname = outputInsideToolDirname + '/' + sizeNumStatDirname;
	mkdir(sizeNumStatDirname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

//	if(max_valid_reg_thres>0){
//		cout << ">>>>>>>>> Before filtering long SV regions: <<<<<<<<<" << endl;
//		outStatScreenFile << ">>>>>>>>> Before filtering long SV regions: <<<<<<<<<" << endl;
//	}
//	SVSizeNumStatOp(user_file, benchmark_file, ref_file, size_div_vec, 0, sizeNumStatDirname);
	if(max_valid_reg_thres>0){
		cout << "\n>>>>>>>>> After filtering long SV regions: <<<<<<<<<" << endl;
		outStatScreenFile << "\n>>>>>>>>> After filtering long SV regions: <<<<<<<<<" << endl;
		SVSizeNumStatOp(user_file, benchmark_file, ref_file, size_div_vec, max_valid_reg_thres, sizeNumStatDirname);
	}
	//Call gnuplot to plot
	SVsizeAndNumstatistics(sizeNumStatDirname, MeticsValues_4);
	SVsizeAndNumstatistics(sizeNumStatDirname, MeticsValues1_4);
}

void SVSizeNumStatOp(string &user_file, string &benchmark_file, string &ref_file, vector<size_t> &size_div_vec, int32_t max_valid_reg_thres, string &dirname){
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
	computeSizeNumStat(divided_vec1, divided_vec2, dirname, ref_file);

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

void computeSizeNumStat(vector< vector<SV_item*> > &divided_vec1, vector< vector<SV_item*> > &divided_vec2, string &dirname, string &ref_file){
	vector<SV_item*> sv_vec1, sv_vec2;
	string file_prefix, file_prefix_tmp;
	size_t start_pos, end_pos;
	faidx_t *fai;

	fai = fai_load(ref_file.c_str());

	file_prefix = dirname + "size_num_stat";

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
		if(i<size_div_vec.size()+1) {
			regionmetric.push_back(region[i]);
			computeNumStat(sv_vec1, sv_vec2, file_prefix_tmp, fai, 4);
			MeticsValues_4.push_back(data_4);
			MeticsValues1_4.push_back(data1_4);
			InformationStorage(i, data_4);
			data_4.clear();
			data1_4.clear();
		}
		else computeNumStatTra(sv_vec1, sv_vec2, file_prefix_tmp);
	}
	allregionmetrics.push_back(regionmetrics);
	if(regionmetrics.size()>1) regionmetrics.erase(regionmetrics.begin()+1,regionmetrics.end());
	fai_destroy(fai);
}

void destroySizeDividedData(vector< vector<SV_item*> > &divided_vec){
	destroyResultData(divided_vec);
}

void InformationStorage(size_t num, vector<float> &Data){
	if(num == 0)	MeticsValues4_0.push_back(Data);
	else if (num == 1) MeticsValues4_1.push_back(Data);
	else if (num == 2) MeticsValues4_2.push_back(Data);
	else if (num == 3) MeticsValues4_3.push_back(Data);
	else if (num == 4) MeticsValues4_4.push_back(Data);
	else if (num == 5) MeticsValues4_5.push_back(Data);
	else if (num == 6) MeticsValues4_6.push_back(Data);
}
