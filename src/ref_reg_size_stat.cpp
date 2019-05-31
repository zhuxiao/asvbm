#include "ref_reg_size_stat.h"
#include "extvab.h"

void refRegSizeStat(string &user_file, string &benchmark_file, int32_t max_valid_reg_thres){
	string refRegSizeFilename_benchmark, refRegSizeFilename_user;
	refRegSizeFilename_benchmark = "ref_reg_size_benchmark";
	refRegSizeFilename_user = "ref_reg_size_user";

	refRegSizeStatDirname = outputPathname + refRegSizeStatDirname;
	mkdir(refRegSizeStatDirname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	cout << ">>>>>>>>> The SV reference region size statistics for benchmark set: <<<<<<<<<" << endl;
	refRegSizeStatOp(refRegSizeFilename_benchmark, benchmark_file, 0, refRegSizeStatDirname);

	if(max_valid_reg_thres>0) cout << ">>>>>>>>> The SV reference region size statistics before filtering for user-called set: <<<<<<<<<" << endl;
	else cout << ">>>>>>>>> The SV reference region size statistics for user-called set: <<<<<<<<<" << endl;
	refRegSizeStatOp(refRegSizeFilename_user, user_file, 0, refRegSizeStatDirname);

	if(max_valid_reg_thres>0){
		cout << "\n>>>>>>>>> The SV reference region size statistics after filtering for user-called set: <<<<<<<<<" << endl;
		refRegSizeStatOp(refRegSizeFilename_user, user_file, max_valid_reg_thres, refRegSizeStatDirname);
	}
}

void refRegSizeStatOp(string &refRegSizeFinename, string &sv_file, int32_t max_valid_reg_thres, string &dirname){
	vector<SV_item*> sv_data, long_sv_data;
	SV_item *item;
	size_t i, count_array[SV_SIZE_ARR_SIZE+2], reg_size, num;
	string refRegSizeFinename_tmp, line;
	ofstream outfile;

	// load data
	sv_data = loadData(sv_file);

	if(max_valid_reg_thres>0) long_sv_data = getLongSVReg(sv_data, max_valid_reg_thres);
	cout << "data.size: " << sv_data.size() << endl;
	if(max_valid_reg_thres>0) cout << "data_long.size: " << long_sv_data.size() << endl;

	for(i=0; i<SV_SIZE_ARR_SIZE+2; i++) count_array[i] = 0;

	for(i=0; i<sv_data.size(); i++){
		item = sv_data.at(i);
		if(item->sv_type!=VAR_TRA and item->sv_type!=VAR_BND and item->sv_type!=VAR_INV_TRA){
			if(item->chrname.size()>0) reg_size = item->endPos - item->startPos + 1;
			else reg_size = item->endPos2 - item->startPos2 + 1;
			if(reg_size<=SV_SIZE_ARR_SIZE) count_array[reg_size] ++;
			else count_array[SV_SIZE_ARR_SIZE+1] ++;
		}
	}
	destroyData(sv_data); // release data

	// compute statistics
	refRegSizeFinename_tmp = dirname + refRegSizeFinename;
	if(max_valid_reg_thres>0) refRegSizeFinename_tmp += "_long_filtered";

	outfile.open(refRegSizeFinename_tmp);
	if(!outfile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << refRegSizeFinename_tmp << endl;
		exit(1);
	}

	num = 0;
	for(i=0; i<SV_SIZE_ARR_SIZE+1; i++){
		line = to_string(i) + "\t" + to_string(count_array[i]);
		outfile << line << endl;
		num ++;
	}
	line = ">=" + to_string(SV_SIZE_ARR_SIZE+1) + "\t" + to_string(count_array[SV_SIZE_ARR_SIZE+1]);
	outfile << line << endl;
	num ++;

	outfile.close();

	cout << num << " reference region size items have been saved to file: " << refRegSizeFinename_tmp << endl;
}


