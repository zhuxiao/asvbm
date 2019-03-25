#include "sv_stat.h"

// global options
int32_t maxValidRegThres;
string outputPathname = "output";
string longSVFilename = "long_sv_reg.bed";
string svSizeDifRatioFilename = "sv_size_dif_ratio";
string svSizeDifStatFilename = "sv_size_dif_stat";
string svSizeRatioStatFilename = "sv_size_ratio_stat";

vector<size_t> size_div_vec = {100, 250, 500, 1000, 5000, 10000};
vector<double> ratio_div_vec = {0.5, 0.7, 1.2, 2, 5, 10, 50, 100};

int main(int argc, char *argv[]){

	string command;

	if (argc < 2) { showUsage(); return 1; }

    if (strcmp(argv[1], "-h") == 0 or strcmp(argv[1], "help") == 0 or strcmp(argv[1], "--help") == 0) {
        if (argc == 2) { showUsage(); return 0; }
        argv++;
        argc = 2;
    }

    if (strcmp(argv[1], "convert")==0){
    	command = "convert";
    	return parseConvert(argc-1, argv+1);
    }else if(strcmp(argv[1], "stat")==0){
    	command = "stat";
    	return parseStat(argc-1, argv+1);
    }else{
    	cerr << "invalid command: " << argv[1] << endl << endl;;
    	showUsage(); return 1;
    }

	return 0;
}

// show usage
void showUsage(){

	cout << "Usage:  sv_stat  <command> [options]" << endl << endl;

	cout << "Commands:" << endl;
	cout << "     convert      convert SV results to 5-column bed file format" << endl;
	cout << "     stat         compute the SV statistics" << endl;
}

// show usage for convert command
void showUsageConvert(){

	cout << "Usage:  sv_stat convert [options] <infile> <outfile>" << endl << endl;

	cout << "Options:" << endl;
	cout << "     -f STR       SV input file format (required):" << endl;
	cout << "                  bed: BED format" << endl;
	cout << "                  vcf: VCF format" << endl;
	cout << "                  csv: CSV format" << endl;
	cout << "     -d INT       remove duplicated variant items [1]: 1 for yes, 0 for no" << endl;
	cout << "     -h           show this help message and exit" << endl;
}

// show usage for stat command
void showUsageStat(){

	cout << "Usage:  sv_stat stat [options] <USER_SV_FILE> <FULL_SV_FILE>" << endl << endl;

	cout << "Description:" << endl;
	cout << "     USER_SV_FILE     User called SV result file." << endl;
	cout << "     FULL_SV_FILE     Full SV set file." << endl << endl;

	cout << "Options:" << endl;
	cout << "     -m           valid maximal region size for statistics: [0]" << endl;
	cout << "                  0 is for all variant size are valid, and while positive" << endl;
	cout << "                  values are for the valid maximal region size, then longer" << endl;
	cout << "                  regions are omitted and saved to the file specified with -l" << endl;
	cout << "     -o FILE      output path name for SV statistics: [output]" << endl;
	cout << "     -l FILE      file name for long SV regions: [long_sv_reg.bed]" << endl;
	cout << "     -h           show this help message and exit" << endl;
}

// parse the parameters for convert command
int parseConvert(int argc, char **argv)
{
	int opt;
	string sv_format, in_file, out_file, remove_dup_str = "1";
	bool remove_dup_falg;

	while( (opt = getopt(argc, argv, ":f:d:h")) != -1 ){
		switch(opt){
			case 'f': sv_format = optarg; break;
			case 'd': remove_dup_str = optarg; break;
			case 'h': showUsageConvert(); exit(0);
			case '?': cout << "unknown option -" << (char)optopt << endl; exit(1);
			case ':': cout << "the option -" << (char)optopt << " needs a value" << endl; exit(1);
		}
	}

	if(remove_dup_str.compare("1")==0) remove_dup_falg = true;
	else if(remove_dup_str.compare("0")==0) remove_dup_falg = false;
	else{
		cout << "Error: Please specify the correct value for '-d' option" << endl << endl;
		showUsageConvert();
		return 1;
	}

	opt = argc - optind; // the number of SAMs on the command line
	if(opt==2) {
		in_file = argv[optind];
		out_file = argv[optind+1];
	}else { showUsageConvert(); return 1; }

	if(sv_format.compare("bed")==0)
		convertBed(in_file, out_file, remove_dup_falg);
	else if(sv_format.compare("vcf")==0)
		convertVcf(in_file, out_file, remove_dup_falg);
	else if(sv_format.compare("csv")==0)
		convertCsv(in_file, out_file, remove_dup_falg);
	else if(sv_format.compare("nm")==0)  // private usage: nm
		convertNm(in_file, out_file, remove_dup_falg);
	else{
		cout << "Error: Please specify the correct SV file format" << endl << endl;
		showUsageConvert();
		return 1;
	}

	return 0;
}

// parse the parameters for stat command
int parseStat(int argc, char **argv){
	int32_t opt;
	string sv_file1, sv_file2;

	if (argc < 2) { showUsageStat(); return 1; }

	maxValidRegThres = 0;
	while( (opt = getopt(argc, argv, ":m:o:l:h")) != -1 ){
		switch(opt){
			case 'm': maxValidRegThres = stoi(optarg); break;
			case 'o': outputPathname = optarg; break;
			case 'l': longSVFilename = optarg; break;
			case 'h': showUsageStat(); exit(0);
			case '?': cout << "unknown option -" << (char)optopt << endl; exit(1);
			case ':': cout << "the option -" << (char)optopt << " needs a value" << endl; exit(1);
		}
	}

	if(maxValidRegThres<0){
		cout << "Error: Please specify the correct positive value for option: -m" << endl;
		showUsageStat();
		return 1;
	}
	if(outputPathname.at(outputPathname.size()-1)!='/') outputPathname += "/";

	opt = argc - optind; // the number of SAMs on the command line
	if(opt==2) {
		sv_file1 = argv[optind];
		sv_file2 = argv[optind+1];
	}else { showUsageStat(); return 1; }

	SVStat(sv_file1, sv_file2);

	return 0;
}

// SV stat
void SVStat(string &sv_file1, string &sv_file2){

	mkdir(outputPathname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	cout << "############# Phage 1: SV size statistics: #############" << endl;
	refRegSizeStat(sv_file1, sv_file2, maxValidRegThres);

	cout << "############# Phage 2: Num statistics: #############" << endl;
	SVNumStat(sv_file1, sv_file2, maxValidRegThres);

	// compute SV size difference statistics
	cout << "\n\n############# Phage 3: SV size difference statistics: #############" << endl;
	SVSizeDifStat(sv_file1, sv_file2, maxValidRegThres);

	cout << "\n\n############# Phage 4: SV type and num statistics: #############" << endl;
	SVTypeNumStat(sv_file1, sv_file2, maxValidRegThres);
}

void refRegSizeStat(string &user_file, string &standard_file, int32_t max_valid_reg_thres){
	string refRegSizeFilename_standard, refRegSizeFilename_user;
	refRegSizeFilename_standard = "ref_reg_size_standard";
	refRegSizeFilename_user = "ref_reg_size_user";

	cout << ">>>>>>>>> The SV reference region size statistics for gold-standard set: <<<<<<<<<" << endl;
	refRegSizeStatOp(refRegSizeFilename_standard, standard_file, 0);

	if(max_valid_reg_thres>0) cout << ">>>>>>>>> The SV reference region size statistics before filtering for user-called set: <<<<<<<<<" << endl;
	else cout << ">>>>>>>>> The SV reference region size statistics for user-called set: <<<<<<<<<" << endl;
	refRegSizeStatOp(refRegSizeFilename_user, user_file, 0);

	if(max_valid_reg_thres>0){
		cout << "\n>>>>>>>>> The SV reference region size statistics after filtering for user-called set: <<<<<<<<<" << endl;
		refRegSizeStatOp(refRegSizeFilename_user, user_file, max_valid_reg_thres);
	}
}

void refRegSizeStatOp(string &refRegSizeFinename, string &sv_file, int32_t max_valid_reg_thres){
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
		if(item->sv_type!=VAR_TRA or item->sv_type!=VAR_TRA or item->sv_type!=VAR_INV_TRA){
			if(item->chrname.size()>0) reg_size = item->endPos - item->startPos + 1;
			else reg_size = item->endPos2 - item->startPos2 + 1;
			if(reg_size<=SV_SIZE_ARR_SIZE) count_array[reg_size] ++;
			else count_array[SV_SIZE_ARR_SIZE+1] ++;
		}
	}

	refRegSizeFinename_tmp = outputPathname + refRegSizeFinename;
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

void SVSizeDifStat(string &sv_file1, string &sv_file2, int32_t max_valid_reg_thres){
	if(max_valid_reg_thres>0)
		cout << ">>>>>>>>> Before filtering long SV regions: <<<<<<<<<" << endl;
	SVSizeDifStatOp(sv_file1, sv_file2, 0);
	if(max_valid_reg_thres>0){
		cout << "\n>>>>>>>>> After filtering long SV regions: <<<<<<<<<" << endl;
		SVSizeDifStatOp(sv_file1, sv_file2, max_valid_reg_thres);
	}
}

void SVSizeDifStatOp(string &sv_file1, string &sv_file2, int32_t max_valid_reg_thres){
	vector<SV_item*> sv_data1, sv_data2, long_sv_data;
	vector<SV_pair*> sv_pair_vec;
	vector< vector<int32_t> > dif_stat_vec;
	vector<size_t> ratio_stat_vec;
	string svSizeDifRatioFilename_tmp, svSizeDifStatFilename_tmp, svSizeRatioStatFilename_tmp;

	// load data
	sv_data1 = loadData(sv_file1);
	sv_data2 = loadData(sv_file2);

	if(max_valid_reg_thres>0)
		long_sv_data = getLongSVReg(sv_data1, max_valid_reg_thres);
	cout << "data1.size: " << sv_data1.size() << endl;
	cout << "data2.size: " << sv_data2.size() << endl;

	// compute overlapped SV pairs
	sv_pair_vec = computeOverlapSVPair(sv_data1, sv_data2);

	// output pair statistics data to file
	svSizeDifRatioFilename_tmp = outputPathname + svSizeDifRatioFilename;
	if(max_valid_reg_thres>0) svSizeDifRatioFilename_tmp += "_long_filtered";
	outputPairDataToFile(svSizeDifRatioFilename_tmp, sv_pair_vec);

	// compute size difference statistics
	dif_stat_vec = computeDifStatVec(sv_pair_vec);

	// save to file
	svSizeDifStatFilename_tmp = outputPathname + svSizeDifStatFilename;
	if(max_valid_reg_thres>0) svSizeDifStatFilename_tmp += "_long_filtered";
	outputDifStatToFile(svSizeDifStatFilename_tmp, dif_stat_vec);

	// compute size ratio statistics
	ratio_stat_vec = computeRatioStatVec(sv_pair_vec, ratio_div_vec);

	// save to file
	svSizeRatioStatFilename_tmp = outputPathname + svSizeRatioStatFilename;
	if(max_valid_reg_thres>0) svSizeRatioStatFilename_tmp += "_long_filtered";
	outputRatioStatToFile(svSizeRatioStatFilename_tmp, ratio_stat_vec, ratio_div_vec);

	// release memory
	destroyPairData(sv_pair_vec);
	destroyData(sv_data1);
	destroyData(sv_data2);
	if(max_valid_reg_thres>0) destroyData(long_sv_data);
}

vector<SV_pair*> computeOverlapSVPair(vector<SV_item*> &data1, vector<SV_item*> &data2){
	vector<SV_pair*> sv_pair_vec;
	SV_item *item1, *item2;
	SV_pair *pair_item;
	size_t i, j;
	int32_t mid_loc1, mid_loc2, reg_size1, reg_size2, dif_tmp;
	double dif_size, size_ratio;

	for(i=0; i<data1.size(); i++){
		item1 = data1.at(i);
		for(j=0; j<data2.size(); j++){
			item2 = data2.at(j);
			if(haveOverlap(item1, item2)){
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
	string line, sv_type_str1, sv_type_str2;
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

		line = to_string(pair_item->dif_size) + "\t" + to_string(pair_item->size_ratio);
		line += "\t" + item1->chrname + "\t" + to_string(item1->startPos) + "\t" + to_string(item1->endPos) + "\t" + sv_type_str1 + "\t" + to_string(item1->sv_len);
		line += "\t" + item2->chrname + "\t" + to_string(item2->startPos) + "\t" + to_string(item2->endPos) + "\t" + sv_type_str2 + "\t" + to_string(item2->sv_len);
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
}

void SVNumStat(string &sv_file1, string &sv_file2, int32_t max_valid_reg_thres){
	if(max_valid_reg_thres>0)
		cout << ">>>>>>>>> Before filtering long SV regions: <<<<<<<<<" << endl;
	SVNumStatOp(sv_file1, sv_file2, 0, outputPathname);

	if(max_valid_reg_thres>0){
		cout << "\n>>>>>>>>> After filtering long SV regions: <<<<<<<<<" << endl;
		SVNumStatOp(sv_file1, sv_file2, max_valid_reg_thres, outputPathname);
	}
}

// SV number stat operation
void SVNumStatOp(string &sv_file1, string &sv_file2, int32_t max_valid_reg_thres, string &dirname){
	vector<SV_item*> sv_data1, sv_data2, long_sv_data;
	string description_str_long, file_prefix, longFilename_tmp;

	//mkdir(dirname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	longFilename_tmp = dirname + longSVFilename;

	sv_data1 = loadData(sv_file1);
	sv_data2 = loadData(sv_file2);

	if(max_valid_reg_thres>0){
		long_sv_data = getLongSVReg(sv_data1, max_valid_reg_thres);
		output2File(longFilename_tmp, long_sv_data);
		cout << "data1.size: " << sv_data1.size() << endl;
		cout << "data2.size: " << sv_data2.size() << endl;
		cout << "long_sv_data.size: " << long_sv_data.size() << endl;
		description_str_long = "Length statistics for long SV regions: ";
		computeLenStat(long_sv_data, description_str_long);
	}else{
		cout << "data1.size: " << sv_data1.size() << endl;
		cout << "data2.size: " << sv_data2.size() << endl;
	}

	file_prefix = dirname + "num_stat";
	computeNumStat(sv_data1, sv_data2, file_prefix);

	destroyData(sv_data1);
	destroyData(sv_data2);
	if(max_valid_reg_thres>0) destroyData(long_sv_data);
}

void computeNumStat(vector<SV_item*> &sv_data1, vector<SV_item*> &sv_data2, string &file_prefix){
	vector<SV_item*> long_sv_data;
	vector< vector<SV_item*> > result;
	int32_t TP, FP, FN, positive_num_called;
	float recall, precision, F1_score, sv_num_per_reg, percent;
	string filename_intersect, filename_private1, filename_private2;

	filename_intersect = file_prefix + "_intersect";
	filename_private1 = file_prefix + "_private1";
	filename_private2 = file_prefix + "_private2";

	// compute intersection
	result = intersect(sv_data1, sv_data2);

	if(sv_data2.size()>0) percent = (double)result.at(0).size() / sv_data2.size();
	else percent = 0;
	cout << "intersection.size: " << result.at(0).size() << ", percent: " << percent << endl;

	if(sv_data1.size()>0) percent = (double)result.at(1).size()/sv_data1.size();
	else percent = 0;
	cout << "private1.size: " << result.at(1).size() << ", percent: " << percent << endl;

	if(sv_data2.size()>0) percent = (double)result.at(2).size()/sv_data2.size();
	else percent = 0;
	cout << "private2.size: " << result.at(2).size() << ", percent: " << percent << endl;

	// output to file
	output2File(filename_intersect, result.at(0));
	output2File(filename_private1, result.at(1));
	output2File(filename_private2, result.at(2));

	//TP = sv_data1.size() - result.at(1).size();
	TP = result.at(0).size();
	FP = result.at(1).size();
	FN = result.at(2).size();
	positive_num_called = sv_data1.size() - result.at(1).size();
	if(sv_data2.size()>0) recall = (float)TP / sv_data2.size();
	else recall = 0;
	//precision = (float)TP / sv_data1.size();
	if(TP+result.at(1).size()>0) precision = (float)TP / (TP + result.at(1).size());
	else precision = 0;

	F1_score = 2.0 * (recall * precision) / (recall + precision);

	if(positive_num_called>0) sv_num_per_reg = (float)TP / positive_num_called;
	else sv_num_per_reg = 0;

	cout << "TP=" << TP << ", FP=" << FP << ", FN=" << FN << endl;
	cout << "Identified SV regions: " << positive_num_called << endl;
	cout << "Recall=" << recall << ", precision=" << precision << ", F1 score=" << F1_score << ", sv_num_per_reg=" << sv_num_per_reg << endl;

	destroyResultData(result);
}

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

	// prepare full sv result
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
		cout << "><><><><><>< i=" << i << endl;
		sv_vec1 = divided_vec1.at(i);
		sv_vec2 = divided_vec2.at(i);

		if(i<size_div_vec.size()){
			end_pos = size_div_vec.at(i);
			file_prefix_tmp = file_prefix + "_" + to_string(start_pos) + "_" + to_string(end_pos);
			start_pos = end_pos + 1;
		}else if(i==size_div_vec.size()){
			file_prefix_tmp = file_prefix + "_" + to_string(start_pos) + "_larger";
		}else if(i==size_div_vec.size()+1){
			file_prefix_tmp = file_prefix + "_TRA_BND";
		}
		computeNumStat(sv_vec1, sv_vec2, file_prefix_tmp);
	}
}

// load data
vector<SV_item*> loadData(const string &filename){
	vector<SV_item*> sv_vec;
	ifstream infile;
	string line, str_tmp;
	vector<string> str_vec;
	SV_item *item;

	infile.open(filename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << filename << endl;
		exit(1);
	}

	// load data
	while(getline(infile, line)){
		if(line.size()){
			str_vec = split(line, "\t");

			item = new SV_item;
			item->chrname = str_vec.at(0);
			if(str_vec.at(1).compare("-")==0) item->startPos = 0;
			else item->startPos = stoi(str_vec.at(1));
			if(str_vec.at(2).compare("-")==0) item->endPos = 0;
			else item->endPos = stoi(str_vec.at(2));
			item->chrname2 = "";
			item->startPos2 = 0;
			item->endPos2 = 0;
			str_tmp = str_vec.at(3);
			if(str_tmp.compare("INS")==0){
				item->sv_type = VAR_INS;
			}else if(str_tmp.compare("DEL")==0){
				item->sv_type = VAR_DEL;
			}else if(str_tmp.compare("DUP")==0){
				item->sv_type = VAR_DUP;
			}else if(str_tmp.compare("INV")==0){
				item->sv_type = VAR_INV;
			}else{
				if(str_vec.size()>=7 and (str_vec.at(6).compare("TRA")==0 or str_vec.at(6).compare("BND")==0)){
					item->chrname2 = str_vec.at(3);
					if(str_vec.at(4).compare("-")==0) item->startPos2 = 0;
					else item->startPos2 = stoi(str_vec.at(4));
					if(str_vec.at(5).compare("-")==0) item->endPos2 = 0;
					else item->endPos2 = stoi(str_vec.at(5));
					if(str_vec.at(6).compare("TRA")==0) item->sv_type = VAR_TRA;
					else item->sv_type = VAR_BND;
				}else
					item->sv_type = VAR_MIX;
			}

			if(item->sv_type==VAR_TRA or item->sv_type==VAR_BND) item->sv_len = 0;
			else item->sv_len = stoi(str_vec.at(4));
			item->overlapped = false;

			sv_vec.push_back(item);
		}
	}
	infile.close();

	return sv_vec;
}

void destroySizeDividedData(vector< vector<SV_item*> > &divided_vec){
	vector<SV_item*> sub_vec;
	for(size_t i=0; i<divided_vec.size(); i++){
		sub_vec = divided_vec.at(i);
		destroyData(sub_vec);
	}
	vector<vector<SV_item*>>().swap(divided_vec);
}

void destroyData(vector<SV_item*> &sv_vec){
	vector<SV_item*>::iterator it;
	for(it=sv_vec.begin(); it!=sv_vec.end(); it++)
		delete *it;
	vector<SV_item*>().swap(sv_vec);
}

void destroyResultData(vector<vector<SV_item*>> &result){
	vector<vector<SV_item*>>::iterator it;
	for(it=result.begin(); it!=result.end(); it++)
		destroyData(*it);
	vector<vector<SV_item*>>().swap(result);
}

void destroyPairData(vector<SV_pair*> &sv_vec){
	vector<SV_pair*>::iterator it;
	for(it=sv_vec.begin(); it!=sv_vec.end(); it++)
		delete *it;
	vector<SV_pair*>().swap(sv_vec);
}

vector<SV_item*> getLongSVReg(vector<SV_item*> &dataset, int32_t thres){
	vector<SV_item*> sv_vec;
	size_t i;
	int32_t reg_len;
	SV_item *sv_item;

	for(i=0; i<dataset.size(); ){
		sv_item = dataset.at(i);
		reg_len = 0;
		if(sv_item->sv_type!=VAR_TRA and sv_item->sv_type!=VAR_BND)
			reg_len = sv_item->endPos - sv_item->startPos + 1;

		if(reg_len>thres){
			sv_vec.push_back(sv_item);
			dataset.erase(dataset.begin()+i);
		}else i++;
	}

	return sv_vec;
}

vector<vector<SV_item*>> intersect(vector<SV_item*> &data1, vector<SV_item*> &data2){
	vector< vector<SV_item*> > result;
	vector<SV_item*> intersect_vec, private_vec1, private_vec2;
	SV_item *item1, *item2, *item;
	size_t i, j;

	for(i=0; i<data1.size(); i++) data1.at(i)->overlapped = false;
	for(i=0; i<data2.size(); i++) data2.at(i)->overlapped = false;

	for(i=0; i<data1.size(); i++){
		item1 = data1.at(i);
		for(j=0; j<data2.size(); j++){
			item2 = data2.at(j);
			if(haveOverlap(item1, item2)){
				item1->overlapped = true;
				item2->overlapped = true;
			}
		}
	}

	for(i=0; i<data2.size(); i++){
		item2 = data2.at(i);
		item = itemdup(item2);
		if(item2->overlapped)
			intersect_vec.push_back(item);
		else
			private_vec2.push_back(item);
	}

	for(i=0; i<data1.size(); i++){
		item1 = data1.at(i);
		if(!item1->overlapped){
			item = itemdup(item1);
			private_vec1.push_back(item);
		}
	}

	result.push_back(intersect_vec);
	result.push_back(private_vec1);
	result.push_back(private_vec2);

	return result;
}

SV_item* itemdup(SV_item* item){
	SV_item *item_new = new SV_item;
	item_new->chrname = item->chrname;
	item_new->startPos = item->startPos;
	item_new->endPos = item->endPos;
	item_new->chrname2 = item->chrname2;
	item_new->startPos2 = item->startPos2;
	item_new->endPos2 = item->endPos2;
	item_new->sv_type = item->sv_type;
	item_new->sv_len = item->sv_len;
	item_new->overlapped = item->overlapped;
	return item_new;
}

bool haveOverlap(SV_item* item1, SV_item* item2){
	bool flag = false, flag_tmp;
	int32_t startPos1, endPos1, startPos2, endPos2, tra_end_idx, tra_end_idx2;

	if(item1->sv_type!=VAR_TRA and item1->sv_type!=VAR_BND){ // INS, DEL, DUP, INV
		if(item1->chrname.compare(item2->chrname)==0){
			startPos1 = item1->startPos - ADJUST_SIZE;
			if(startPos1<1) startPos1 = 1;
			endPos1 = item1->endPos + ADJUST_SIZE;
			startPos2 = item2->startPos - ADJUST_SIZE;
			if(startPos2<1) startPos2 = 1;
			endPos2 = item2->endPos + ADJUST_SIZE;

			flag = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
		}

		if(flag==false and (item2->sv_type==VAR_TRA or item2->sv_type==VAR_BND)){
			if(item1->chrname.compare(item2->chrname2)==0){
				startPos1 = item1->startPos - ADJUST_SIZE;
				if(startPos1<1) startPos1 = 1;
				endPos1 = item1->endPos + ADJUST_SIZE;
				startPos2 = item2->startPos2 - ADJUST_SIZE;
				if(startPos2<1) startPos2 = 1;
				endPos2 = item2->endPos2 + ADJUST_SIZE;

				flag = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
			}
		}
	}else{ // TRA or BND
		tra_end_idx = -1;
		if(item1->chrname.compare(item2->chrname)==0){
			startPos1 = item1->startPos - ADJUST_SIZE;
			if(startPos1<1) startPos1 = 1;
			endPos1 = item1->endPos + ADJUST_SIZE;
			startPos2 = item2->startPos - ADJUST_SIZE;
			if(startPos2<1) startPos2 = 1;
			endPos2 = item2->endPos + ADJUST_SIZE;

			flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
			if(flag_tmp) tra_end_idx = 1;
		}
		if(flag_tmp==false and item1->chrname.compare(item2->chrname2)==0){
			startPos1 = item1->startPos - ADJUST_SIZE;
			if(startPos1<1) startPos1 = 1;
			endPos1 = item1->endPos + ADJUST_SIZE;
			startPos2 = item2->startPos2 - ADJUST_SIZE;
			if(startPos2<1) startPos2 = 1;
			endPos2 = item2->endPos2 + ADJUST_SIZE;

			flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
			if(flag_tmp) tra_end_idx = 2;
		}

		tra_end_idx2 = -1;
		if(item1->chrname2.compare(item2->chrname)==0){
			startPos1 = item1->startPos2 - ADJUST_SIZE;
			if(startPos1<1) startPos1 = 1;
			endPos1 = item1->endPos2 + ADJUST_SIZE;
			startPos2 = item2->startPos - ADJUST_SIZE;
			if(startPos2<1) startPos2 = 1;
			endPos2 = item2->endPos + ADJUST_SIZE;

			flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
			if(flag_tmp) tra_end_idx2 = 1;
		}
		if(flag_tmp==false and item1->chrname2.compare(item2->chrname2)==0){
			startPos1 = item1->startPos2 - ADJUST_SIZE;
			if(startPos1<1) startPos1 = 1;
			endPos1 = item1->endPos2 + ADJUST_SIZE;
			startPos2 = item2->startPos2 - ADJUST_SIZE;
			if(startPos2<1) startPos2 = 1;
			endPos2 = item2->endPos2 + ADJUST_SIZE;

			flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
			if(flag_tmp) tra_end_idx2 = 2;
		}

		if(tra_end_idx!=-1 or tra_end_idx2!=-1)
			flag = true;
	}

	return flag;
}

// determine whether the given positions of two regions are overlapped
bool isOverlappedPos(size_t startPos1, size_t endPos1, size_t startPos2, size_t endPos2){
	bool flag = false;
	if((startPos1>=startPos2 and startPos1<=endPos2)
		or (endPos2>=startPos1 and endPos2<=endPos1)
		or (startPos2>=startPos1 and startPos2<=endPos1)
		or (endPos1>=startPos2 and endPos1<=endPos2))
			flag = true;
	return flag;
}

void output2File(const string &filename, vector<SV_item*> &data){
	ofstream outfile;
	string data_out, sv_type_str;
	SV_item *item;

	outfile.open(filename);
	if(!outfile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << filename << endl;
		exit(1);
	}

	for(size_t i=0; i<data.size(); i++){
		item = data.at(i);

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
			default:
				cerr << "line=" << __LINE__ << ", invalid sv type: " << item->sv_type << endl;
				exit(1);
		}

		if(item->sv_type!=VAR_TRA and item->sv_type!=VAR_BND)
			data_out = item->chrname + "\t" + to_string(item->startPos) + "\t" + to_string(item->endPos) + "\t" + sv_type_str + "\t" + to_string(item->sv_len);
		else
			data_out = item->chrname + "\t" + to_string(item->startPos) + "\t" + to_string(item->endPos) + "\t" + item->chrname2 + "\t" + to_string(item->startPos2) + "\t" + to_string(item->endPos2) + "\t" + sv_type_str + "\t" + to_string(item->sv_len);
		outfile << data_out << endl;
	}

	outfile.close();
}

void computeLenStat(vector<SV_item*> &data, string &description_str){
	int64_t len, sum, max_len, min_len, aver_len;
	SV_item *item;
	sum = 0; max_len = 0; min_len = INT_MAX;
	for(size_t i=0; i<data.size(); i++){
		item = data.at(i);
		len = item->endPos - item->startPos + 1;
		sum += len;
		if(max_len<len) max_len = len;
		if(min_len>len) min_len = len;
	}
	if(data.size()>0) aver_len = sum / data.size();
	else { aver_len = 0; min_len = 0; }
	cout << description_str << " reg_num=" << data.size() << ", total_len=" << sum << ", max_len=" << max_len << ", min_len=" << min_len << ", aver_len=" << aver_len << endl;
}

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
	output2File(outfilename, sv_item_vec);

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
	output2File(outfilename, sv_item_vec);

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

					chrname2 = "";
					start_pos2 = endpos2 = 0;

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
	output2File(outfilename, sv_item_vec);

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
	infile.close();

	// remove duplicated sv items
	if(removeDuplicatedItemFlag) removeDuplicatedSVItems(sv_item_vec);

	// save to file
	output2File(outfilename, sv_item_vec);

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
	}else if(sv_type_str.compare("INS")==0){
		sv_type = VAR_INS;
	}else if(sv_type_str.compare("DEL")==0){
		sv_type = VAR_DEL;
	}else if(sv_type_str.compare("DUP")==0){
		sv_type = VAR_DUP;
	}else if(sv_type_str.compare("INV")==0){
		sv_type = VAR_INV;
	}else if(sv_type_str.compare("TRA")==0){
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

// string split function
vector<string> split(const string& s, const string& delim)
{
    vector<string> elems;
    size_t pos = 0;
    size_t len = s.length();
    size_t delim_len = delim.length();
    if (delim_len == 0) return elems;
    while (pos < len)
    {
        int find_pos = s.find(delim, pos);
        if (find_pos < 0)
        {
            elems.push_back(s.substr(pos, len - pos));
            break;
        }
        elems.push_back(s.substr(pos, find_pos - pos));
        pos = find_pos + delim_len;
    }
    return elems;
}

void printSV(vector<SV_item*> &data){
	string data_out, sv_type_str;
	SV_item *item;

	for(size_t i=0; i<data.size(); i++){
		item = data.at(i);

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
			default:
				cerr << "line=" << __LINE__ << ", invalid sv type: " << item->sv_type << endl;
				exit(1);
		}

		if(item->sv_type!=VAR_TRA and item->sv_type!=VAR_BND)
			data_out = item->chrname + "\t" + to_string(item->startPos) + "\t" + to_string(item->endPos) + "\t" + sv_type_str + "\t" + to_string(item->sv_len);
		else
			data_out = item->chrname + "\t" + to_string(item->startPos) + "\t" + to_string(item->endPos) + "\t" + item->chrname2 + "\t" + to_string(item->startPos2) + "\t" + to_string(item->endPos2) + "\t" + sv_type_str + "\t" + to_string(item->sv_len);
		cout << data_out << endl;
	}
}

// transform double to string with specified fixed precision
string double2Str(double num, int32_t fixedPrecision){
	string num_str, tmp_str;
	int32_t i, dot_idx, float_bit_num, len;
	char ch;

	tmp_str = to_string(num);
	dot_idx = tmp_str.find(".");
	if(dot_idx>0){
		float_bit_num = (int32_t)tmp_str.size() - 1 - dot_idx;
		if(float_bit_num>fixedPrecision){
			ch = tmp_str.at(dot_idx+fixedPrecision+1);
			len = dot_idx + fixedPrecision + 1;
			if(ch<'5') num_str = tmp_str.substr(0, len);
			else tmp_str.at(dot_idx+fixedPrecision) ++;
		}else{
			num_str = tmp_str;
			len = fixedPrecision - float_bit_num;
			for(i=0; i<len; i++) num_str += "0";
		}
	}else{
		num_str = tmp_str;
		if(fixedPrecision>0){
			num_str += ".";
			for(i=0; i<fixedPrecision; i++) num_str += "0";
		}
	}

	return num_str;
}
