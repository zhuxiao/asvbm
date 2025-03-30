#include "dataLoader.h"
#include "num_stat.h"
#include "extvab.h"
#include <htslib/faidx.h>
#include "convert.h"


void SVNumStat(string &user_file, string &benchmark_file, string &ref_file, int32_t max_valid_reg_thres, string &outputPathname, vector<string> &sv_files1){
	if(sv_files1.size()>=1) numStatDirname = outputInsideToolDirname + '/' + numStatDirname;
	else numStatDirname = outputPathname + numStatDirname;
	mkdir(numStatDirname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	 //non-TRA
//	if(max_valid_reg_thres>0){
//		cout << ">>>>>>>>> Before filtering long SV regions: <<<<<<<<<" << endl;
//		outStatScreenFile << ">>>>>>>>> Before filtering long SV regions: <<<<<<<<<" << endl;
//	}
//	SVNumStatOp(user_file, benchmark_file, ref_file, 0, numStatDirname);

	if(max_valid_reg_thres>0){
		cout << "\n>>>>>>>>> After filtering long SV regions: <<<<<<<<<" << endl;
		outStatScreenFile << "\n>>>>>>>>> After filtering long SV regions: <<<<<<<<<" << endl;
		SVNumStatOp(user_file, benchmark_file, ref_file, max_valid_reg_thres, numStatDirname);
	}
	MeticsValues.push_back(data);
	MeticsValues1.push_back(data1);

	// TRA
//	cout << "\n>>>>>>>>> TRA breakpoint statistics: <<<<<<<<<" << endl;
//	outStatScreenFile << "\n>>>>>>>>> TRA breakpoint statistics: <<<<<<<<<" << endl;
//	SVNumStatTraOp(user_file, benchmark_file, numStatDirname);
}

// SV number stat operation
void SVNumStatOp(string &user_file, string &benchmark_file, string &ref_file, int32_t max_valid_reg_thres, string &dirname){
	vector<SV_item*> user_data, benchmark_data, long_sv_data;
	string description_str_long, file_prefix, longFilename_tmp;
	faidx_t *fai;

	fai = fai_load(ref_file.c_str());

	longFilename_tmp = dirname + longSVFilename;

	// non-TRA variants
	user_data = loadDataWithoutTra(user_file);
	benchmark_data = loadDataWithoutTra(benchmark_file);

	if(max_valid_reg_thres>0){
		long_sv_data = getLongSVReg(user_data, max_valid_reg_thres);
		output2File(longFilename_tmp, long_sv_data, outStatScreenFile);
		cout << "Non-TRA: user data size: " << user_data.size() << endl;
		cout << "Non-TRA: benchmark data size: " << benchmark_data.size() << endl;
		cout << "Non-TRA: long_sv_data.size: " << long_sv_data.size() << endl;
		outStatScreenFile << "Non-TRA: user data size: " << user_data.size() << endl;
		outStatScreenFile << "Non-TRA: benchmark data size: " << benchmark_data.size() << endl;
		outStatScreenFile << "Non-TRA: long_sv_data.size: " << long_sv_data.size() << endl;
		description_str_long = "Size statistics for long SV regions: ";
		computeLenStat(long_sv_data, description_str_long);
	}else{
		cout << "Non-TRA: user data size: " << user_data.size() << endl;
		cout << "Non-TRA: benchmark data size: " << benchmark_data.size() << endl;
		outStatScreenFile << "Non-TRA: user data size: " << user_data.size() << endl;
		outStatScreenFile << "Non-TRA: benchmark data size: " << benchmark_data.size() << endl;
	}

	file_prefix = dirname + "num_stat";
	computeNumStat(user_data, benchmark_data, file_prefix, fai, 2);

	destroyData(user_data);
	destroyData(benchmark_data);
	fai_destroy(fai);
	if(max_valid_reg_thres>0) destroyData(long_sv_data);
}

void computeNumStat(vector<SV_item*> &user_data, vector<SV_item*> &benchmark_data, string &file_prefix, faidx_t *fai, int Markers){
	vector< vector<SV_item*> > result;
	int32_t TP_benchmark, TP_user, FP, FN ,LP;
	float recall, precision_user, F1_score_user, sv_num_per_reg, percent;
	string filename_intersect_user, filename_latentpositive_user, filename_intersect_benchmark, filename_private_user, filename_private_benchmark, out_str;
	string filename_TP_user, filename_TP_bench, filename_FP, filename_FN, filename_LP;
	SeqconsSum = SeqconsNum = Seqcons = LpNum = 0;

	filename_intersect_user = file_prefix + "_intersect_user";
	filename_intersect_benchmark = file_prefix + "_intersect_benchmark";
	filename_private_user = file_prefix + "_private_user";
	filename_private_benchmark = file_prefix + "_private_benchmark";
	filename_latentpositive_user = file_prefix + "_latentpositive_user";

	filename_TP_user = file_prefix + "_TP_user.vcf";
	filename_TP_bench = file_prefix + "_TP_bench.vcf";
	filename_FP = file_prefix + "_FP.vcf";
	filename_FN = file_prefix + "_FN.vcf";
	filename_LP = file_prefix + "_LP.vcf";

	sharedFPFilenames.push_back(filename_FP);
	benchmarkAddTP_user.push_back(filename_intersect_user);
	// compute intersection
	result = intersect(user_data, benchmark_data, fai);

	if(user_data.size()>0) percent = (double)result.at(0).size() / user_data.size();
	else percent = 0;
	out_str = "User intersection data size: " + to_string(result.at(0).size()) + ", percent: " + to_string(percent);
	cout << out_str << endl;
	outStatScreenFile << out_str << endl;

	if(benchmark_data.size()>0) percent = (double)result.at(1).size() / benchmark_data.size();
	else percent = 0;
	out_str = "Benchmark intersection data size: " + to_string(result.at(1).size()) + ", percent: " + to_string(percent);
	cout << out_str << endl;
	outStatScreenFile << out_str << endl;

//	if(user_data.size()>0) percent = (double)result.at(1).size()/user_data.size();
//	else percent = 0;
//	out_str = "User private data size: " + to_string(result.at(1).size()) + ", percent: " + to_string(percent);
//	cout << out_str << endl;
//	outStatScreenFile << out_str << endl;

	if(user_data.size()>0) percent = (double)result.at(2).size()/user_data.size();
		else percent = 0;
		out_str = "User private data size: " + to_string(result.at(2).size()) + ", percent: " + to_string(percent);
		cout << out_str << endl;
		outStatScreenFile << out_str << endl;

//	if(benchmark_data.size()>0) percent = (double)result.at(2).size()/benchmark_data.size();
//	else percent = 0;
//	out_str = "Benchmark private data size: " + to_string(result.at(2).size()) + ", percent: " + to_string(percent);
//	cout << out_str << endl;
//	outStatScreenFile << out_str << endl;

	if(benchmark_data.size()>0) percent = (double)result.at(3).size()/benchmark_data.size();
	else percent = 0;
	out_str = "Benchmark private data size: " + to_string(result.at(3).size()) + ", percent: " + to_string(percent);
	cout << out_str << endl;
	outStatScreenFile << out_str << endl;

	// output to file
	output3File(filename_intersect_benchmark, result.at(1), outStatScreenFile);
	output3File(filename_intersect_user, result.at(0), outStatScreenFile);
	output3File(filename_private_user, result.at(2), outStatScreenFile);
	output3File(filename_private_benchmark, result.at(3), outStatScreenFile);
	output3File(filename_latentpositive_user, result.at(4), outStatScreenFile);

	if(usersets_num <= (int32_t)SVcallernames.size()){
		outputvcfFile(filename_TP_bench, result.at(1));
		outputvcfFile(filename_TP_user, result.at(0));
		outputvcfFile(filename_FP, result.at(2));
		outputvcfFile(filename_FN, result.at(3));
		outputvcfFile(filename_LP, result.at(4));
		usersets_num +=1;
	}

	if(Markers==2){
		TPfilesPath.push_back(filename_intersect_benchmark);
		FNfilesPath.push_back(filename_private_benchmark);
		ResultdataPath.push_back(filename_intersect_benchmark);
		ResultdataPath.push_back(filename_intersect_user);
		ResultdataPath.push_back(filename_private_user);
		ResultdataPath.push_back(filename_private_benchmark);
	}
	//TP_user = user_data.size() - result.at(0).size();
	//TP_benchmark = user_data.size() - result.at(1).size();
	TP_user = result.at(0).size();
	TP_benchmark = result.at(1).size();
	FP = result.at(2).size();
	FN = result.at(3).size();
	LP = result.at(4).size();
//	cout<<"11111111111:  "<< LpNum <<endl;
	if(benchmark_data.size()>0) recall = (float)TP_benchmark / benchmark_data.size();
	else recall = 0;
	//precision_user = (float)TP_user / user_data.size();
	if(TP_benchmark+result.at(2).size()>0) precision_user = (float)(TP_benchmark) / (TP_benchmark + result.at(2).size());
	else precision_user = 0;
	//precision_benchmark = (float)TP_benchmark / user_data.size();
	/*if(TP_benchmark+result.at(2).size()>0) precision_benchmark = (float)TP_benchmark / (TP_benchmark + result.at(2).size());
	else precision_benchmark = 0;*/

	if(recall+precision_user>0) F1_score_user = 2.0 * (recall * precision_user) / (recall + precision_user);
	else F1_score_user = 0;
	/*if(recall+precision_benchmark>0) F1_score_benchmark = 2.0 * (recall * precision_benchmark) / (recall + precision_benchmark);
	else F1_score_benchmark = 0;*/

	if(TP_user>0) sv_num_per_reg = (float)TP_benchmark / TP_user;
	else sv_num_per_reg = 0;

	if(SeqconsNum != 0)	Seqcons = SeqconsSum/(double)SeqconsNum;
	else	Seqcons = 0;
//	cout << "11111111:" << Seqcons << "	num:" << SeqconsNum << " SeqconsSum:" << SeqconsSum << endl;

	cout << "TP_benchmark=" << TP_benchmark << ", TP_user=" << TP_user << ", FP=" << FP << ", FN=" << FN << ", LP=" << LP << endl;
	cout << "Recall=" << recall << ", Precision=" << precision_user << ", F1 score=" << F1_score_user << endl; //", precision_benchmark=" << precision_benchmark <<
//	cout << "precision_user=" << precision_user << ", F1 score_user=" << F1_score_user << endl;
	cout <<"Identity=" << Seqcons << ", sv_num_per_reg=" << sv_num_per_reg << endl; // ", F1 score_user=" << F1_score_user <<

	outStatScreenFile << "TP_benchmark=" << TP_benchmark << ", TP_user=" << TP_user << ", FP=" << FP << ", FN=" << FN << ", LP=" << LP << endl;
	outStatScreenFile << "Recall=" << recall << ", Precision=" << precision_user << ", F1 score=" << F1_score_user << endl;
//	outStatScreenFile << "precision_user=" << precision_user << ", F1 score_user=" << F1_score_user << endl;
	outStatScreenFile << "Identity=" << Seqcons << ", sv_num_per_reg=" << sv_num_per_reg << endl;
	if(Markers == 2){
		CollectData(LP, TP_benchmark, FP, FN, data1, 4, 2);
		CollectData(recall, precision_user, F1_score_user, Seqcons, data, 3, 2);
	}
	allmetric.clear();
	if(Markers == 4){
		CollectData(LP, TP_benchmark, FP, FN, data1_4, 4, 4);
		CollectData(recall, precision_user, F1_score_user, Seqcons, data_4, 3, 4);
	}
	regionmetric.clear();
	destroyResultData(result);
}

void computeNumStatFromFile(vector<SV_item*> TPbench_data, vector<SV_item*> TPuser_data, vector<SV_item*> FP_data, vector<SV_item*> FN_data, string &file_prefix, int32_t endpos, int Markers){
	vector<SV_item*> sub_TP_bench, sub_TP_user,sub_FN, sub_FP;
	int32_t TP_benchmark, TP_user, FP, FN;
	float recall, precision_benchmark, precision_user, F1_score_benchmark, F1_score_user, sv_num_per_reg, percent;
	string filename_intersect_user, filename_intersect_benchmark, filename_private_user, filename_private_benchmark, out_str;
	int32_t startpos;
	SV_item *item;
	SeqconsSum = SeqconsNum = Seqcons = 0;

	filename_intersect_user = file_prefix + "_intersect_user";
	filename_intersect_benchmark = file_prefix + "_intersect_benchmark";
	filename_private_user = file_prefix + "_private_user";
	filename_private_benchmark = file_prefix + "_private_benchmark";

	// compute intersection
//	result = intersect(user_data, benchmark_data, fai);
	for(size_t i=0; i<=size_div_vec.size(); ++i){
		if(size_div_vec[i] == (size_t)endpos){
			if(i==0)	startpos = 1;
			else if(sign) startpos = size_div_vec[i] + 1;
			else if(i==size_div_vec.size()-1){
				if((size_t)endpos == size_div_vec[size_div_vec.size()-1]) sign = true;
				startpos = size_div_vec[i-1]+1;
			}else	startpos = size_div_vec[i-1] + 1;
			break;
		}
	}
	//TP_bench
	for(size_t i=0;i<TPbench_data.size();i++){
		if(startpos < endpos){
			item = TPbench_data[i];
			if(item->sv_len >= startpos and item->sv_len <= endpos){
				sub_TP_bench.push_back(item);
				if(item->seqcons.compare("-")!=0){
					double value = stod((item->seqcons).substr((item->seqcons).find('=') + 1));
					SeqconsSum += value;
					SeqconsNum += 1;
				}
			}
		}else{
			item = TPbench_data[i];
			if(item->sv_len >= startpos){
				sub_TP_bench.push_back(item);
				if(item->seqcons.compare("-")!=0){
					double value = stod((item->seqcons).substr((item->seqcons).find('=') + 1));
					SeqconsSum += value;
					SeqconsNum += 1;
				}
			}
		}
	}
	//TP_user
	for(size_t i=0;i<TPuser_data.size();i++){
		if(startpos < endpos){
			item = TPuser_data[i];
			if(item->sv_len >= startpos and item->sv_len <= endpos)	sub_TP_user.push_back(item);
		}else{
			item = TPuser_data[i];
			if(item->sv_len >= startpos)	sub_TP_user.push_back(item);
		}
	}
	//FP
	for(size_t i=0;i<FP_data.size();i++){
		if(startpos < endpos){
			item = FP_data[i];
			if(item->sv_len >= startpos and item->sv_len <= endpos)	sub_FP.push_back(item);
		}else{
			item = FP_data[i];
			if(item->sv_len >= startpos)	sub_FP.push_back(item);
		}
	}
	//FN
	for(size_t i=0;i<FN_data.size();i++){
		if(startpos < endpos){
			item = FN_data[i];
			if(item->sv_len >= startpos and item->sv_len <= endpos)	sub_FN.push_back(item);
		}else{
			item = FN_data[i];
			if(item->sv_len >= startpos)	sub_FN.push_back(item);
		}
	}

	if(sub_TP_user.size() + sub_FP.size()>0) percent = (double)sub_TP_user.size() / (sub_TP_user.size() + sub_FP.size());
	else percent = 0;
	out_str = "User intersection data size: " + to_string(sub_TP_user.size()) + ", percent: " + to_string(percent);
	cout << out_str << endl;
	outStatScreenFile << out_str << endl;

	if(sub_TP_bench.size() + sub_FN.size()>0) percent = (double)sub_TP_bench.size() / (sub_TP_bench.size() + sub_FN.size());
	else percent = 0;
	out_str = "Benchmark intersection data size: " + to_string(sub_TP_bench.size()) + ", percent: " + to_string(percent);
	cout << out_str << endl;
	outStatScreenFile << out_str << endl;

	if(sub_TP_user.size() + sub_FP.size()>0) percent = (double)sub_FP.size()/(sub_TP_user.size() + sub_FP.size());
	else percent = 0;
	out_str = "User private data size: " + to_string(sub_FP.size()) + ", percent: " + to_string(percent);
	cout << out_str << endl;
	outStatScreenFile << out_str << endl;

	if(sub_TP_bench.size() + sub_FN.size()>0) percent = (double)sub_FN.size()/(sub_TP_bench.size() + sub_FN.size());
	else percent = 0;
	out_str = "Benchmark private data size: " + to_string(sub_FN.size()) + ", percent: " + to_string(percent);
	cout << out_str << endl;
	outStatScreenFile << out_str << endl;

	// output to file
	output2File(filename_intersect_benchmark, sub_TP_bench, outStatScreenFile);
	output2File(filename_intersect_user, sub_TP_user, outStatScreenFile);
	output2File(filename_private_user, sub_FP, outStatScreenFile);
	output2File(filename_private_benchmark, sub_FN, outStatScreenFile);
//	if(Markers==2){
//		TPfilesPath.push_back(filename_intersect_benchmark);
//		FNfilesPath.push_back(filename_private_benchmark);
//	}
	//TP_user = user_data.size() - result.at(0).size();
	//TP_benchmark = user_data.size() - result.at(1).size();
	TP_user = sub_TP_user.size();
	TP_benchmark = sub_TP_bench.size();
	FP = sub_FP.size();
	FN = sub_FN.size();

	if(TP_benchmark+FN>0) recall = (float)TP_benchmark / (TP_benchmark+FN);
	else recall = 0;
	//precision_user = (float)TP_user / user_data.size();
	if(TP_user+FP>0) precision_user = (float)TP_user / (TP_user + FP);
	else precision_user = 0;
	//precision_benchmark = (float)TP_benchmark / user_data.size();
	if(TP_benchmark+FP>0) precision_benchmark = (float)TP_benchmark / (TP_benchmark + FP);
	else precision_benchmark = 0;

	if(recall+precision_user>0) F1_score_user = 2.0 * (recall * precision_user) / (recall + precision_user);
	else F1_score_user = 0;
	if(recall+precision_benchmark>0) F1_score_benchmark = 2.0 * (recall * precision_benchmark) / (recall + precision_benchmark);
	else F1_score_benchmark = 0;

	if(TP_user>0) sv_num_per_reg = (float)TP_benchmark / TP_user;
	else sv_num_per_reg = 0;

	if(SeqconsNum != 0)	Seqcons = SeqconsSum/(double)SeqconsNum;
	else	Seqcons = 0;

	cout << "TP_benchmark=" << TP_benchmark << ", TP_user=" << TP_user << ", FP=" << FP << ", FN=" << FN << endl;
	cout << "Recall=" << recall << ", precision_benchmark=" << precision_benchmark << ", precision_user=" << precision_user << endl;
	cout << "F1 score_benchmark=" << F1_score_benchmark << ", F1 score_user=" << F1_score_user << ", Identity=" << Seqcons << ", sv_num_per_reg=" << sv_num_per_reg << endl;

	outStatScreenFile << "TP_benchmark=" << TP_benchmark << ", TP_user=" << TP_user << ", FP=" << FP << ", FN=" << FN << endl;
	outStatScreenFile << "Recall=" << recall << ", precision_benchmark=" << precision_benchmark << ", precision_user=" << precision_user << endl;
	outStatScreenFile << "F1 score_benchmark=" << F1_score_benchmark << ", F1 score_user=" << F1_score_user << ", Identity=" << Seqcons << ", sv_num_per_reg=" << sv_num_per_reg << endl;
	if(Markers == 2){
		CollectData(TP_user, TP_benchmark, FP, FN, data1, 4, 2);
		CollectData(recall, precision_user, F1_score_user, Seqcons, data, 3, 2);
	}
	allmetric.clear();
	if(Markers == 4){
		CollectData(TP_user, TP_benchmark, FP, FN, data1_4, 4, 4);
		CollectData(recall, precision_user, F1_score_user, Seqcons, data_4, 3, 4);
	}
	regionmetric.clear();
}

void CollectData(float recall, float precision_user, float F1_score_user, double seqcons, vector<float> &Data, size_t num, int Markers){
	Data.push_back(seqcons);
	Data.push_back(recall);
	Data.push_back(precision_user);
	Data.push_back(F1_score_user);
//	Data.push_back(seqcons);
	if(Markers==2){
		allmetric.push_back(to_string(seqcons));
		allmetric.push_back(to_string(recall));
		allmetric.push_back(to_string(precision_user));
		allmetric.push_back(to_string(F1_score_user));
//		allmetric.push_back(to_string(seqcons));
		allmetrics.push_back(allmetric);
	}else if(Markers==4){
		regionmetric.push_back(to_string(seqcons));
		regionmetric.push_back(to_string(recall));
		regionmetric.push_back(to_string(precision_user));
		regionmetric.push_back(to_string(F1_score_user));
//		regionmetric.push_back(to_string(seqcons));
		regionmetrics.push_back(regionmetric);
	}
//	if(Data.size()>num){
//		Data.clear();
//		Data = {recall, precision_user, F1_score_user};
//	}
}
void CollectData(int LP, int TP_benchmark, int FP, int FN, vector<int> &Data, size_t num, int Markers){
	Data.push_back(TP_benchmark);
	Data.push_back(FP);
	Data.push_back(FN);
	Data.push_back(LP);
	if(Markers==2){
		allmetric.push_back(to_string(TP_benchmark));
		allmetric.push_back(to_string(FP));
		allmetric.push_back(to_string(FN));
		allmetric.push_back(to_string(LP));
	}else if (Markers==4){
		regionmetric.push_back(to_string(TP_benchmark));
		regionmetric.push_back(to_string(FP));
		regionmetric.push_back(to_string(FN));
		regionmetric.push_back(to_string(LP));
	}
//	if(Data.size()>num){
//		data1.clear();
//		data1 = {TP_user, TP_benchmark, FP, FN};
//		}
}


void computeLenStat(vector<SV_item*> &data, string &description_str){
	int64_t len, sum, max_len, min_len, aver_len;
	SV_item *item;
	string out_str;

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

	out_str = description_str + " reg_num=" + to_string(data.size()) + ", total_len=" + to_string(sum) + ", max_len=" + to_string(max_len) + ", min_len=" + to_string(min_len) + ", aver_len=" + to_string(aver_len);
	cout << out_str << endl;
	outStatScreenFile << out_str << endl;
}

vector<vector<SV_item*>> intersect(vector<SV_item*> &user_data, vector<SV_item*> &benchmark_data, faidx_t *fai){
	vector< vector<SV_item*> > result, subsets;

	subsets = constructSubsetByChr(user_data, benchmark_data);
	sortSubsets(subsets); // sort
	checkOrder(subsets); // check order

	result = intersectOp(subsets, fai);

	return result;
}

//generate two datasets
vector<vector<SV_item*>> constructSubsetByChr(vector<SV_item*> &user_data, vector<SV_item*> &benchmark_data){
	vector< vector<SV_item*> > result;
	set<string> chrname_set1, chrname_set2;
	set<string> chrname_set_tmp;
	vector<string> chrname_vec_sorted;

	chrname_set1 = getChrnames(user_data);
	chrname_set2 = getChrnames(benchmark_data);
	chrname_set_tmp = getChrUnion(chrname_set1, chrname_set2);
	chrname_vec_sorted = sortChrnames(chrname_set_tmp); // sort chromosome names

	result = constructSubsetByChrOp(user_data, benchmark_data, chrname_vec_sorted);

	return result;
}

//get chrnames
set<string> getChrnames(vector<SV_item*> &dataset){
	set<string> chrname_set;
	size_t i;
	SV_item *item;
	set<string>::iterator iter;

	for(i=0;i<dataset.size();i++){
		item = dataset.at(i);
		if(chrname_set.find(item->chrname)==chrname_set.end()) chrname_set.insert(item->chrname);
	}

	return chrname_set;
}

//union
set<string> getChrUnion(set<string> &chrname_set1, set<string> &chrname_set2){
	set<string> chrname_union;
	set<string>::iterator iter, iter1;
	string str1, str2, str3, tmp_str, chr = "chr";
	size_t len = chr.size();
	int32_t cond_idx;

	for(iter=chrname_set1.begin(); iter!=chrname_set1.end(); ++iter) chrname_union.insert(*iter);

	for(iter=chrname_set2.begin(); iter!=chrname_set2.end(); ++iter){
		str1 = *iter;
		str2 = "";
		str3 = "";
		if((*iter).size()>len){
			tmp_str = (*iter).substr(0,len);
			if(tmp_str.compare(chr)!=0) str2 = chr + *iter;
			else str3 = (*iter).substr(len);
		}else{
			str2 = chr + *iter;
		}

		cond_idx = -1;
		if(chrname_union.find(str1)!=chrname_union.end()) cond_idx = 0;
		if(cond_idx==-1 and str2.size()>0 and chrname_union.find(str2)!=chrname_union.end()) cond_idx = 1;
		if(cond_idx==-1 and str3.size()>0 and chrname_union.find(str3)!=chrname_union.end()) cond_idx = 2;
		if(cond_idx==-1) chrname_union.insert(*iter);//not found
	}

	return chrname_union;
}

vector<string> sortChrnames(set<string> &chrname_set){
	int8_t *selected_flag_array;
	size_t i, j;
	int32_t idx_chr;
	set<string>::iterator iter;
	string chr_str1, chr_str2, head_str, head_str_chr, chrname;
	vector<string> chrname_vec_sorted, chr_str_vec, chrname_vec;
	bool have_chr_prefix_flag;

	chr_str1 = "chr1_chr2_chr3_chr4_chr5_chr6_chr7_chr8_chr9_chr10_chr11_chr12_chr13_chr14_chr15_chr16_chr17_chr18_chr19_chr20_chr21_chr22_chrX_chrY_chrM";
	chr_str2 = "1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19_20_21_22_X_Y_MT";

	for(iter=chrname_set.begin(); iter!=chrname_set.end(); ++iter) chrname_vec.push_back(*iter);
	selected_flag_array = (int8_t*) calloc(chrname_vec.size(), sizeof(int8_t));

	// determine the have_chr_prefix_flag
	have_chr_prefix_flag = false;
	for(i=0; i<chrname_vec.size(); i++){
		if(chrname_vec.at(i).compare("chr1")==0){
			have_chr_prefix_flag = true;
			break;
		}
	}

	if(have_chr_prefix_flag) chr_str_vec = split(chr_str1, "_");
	else chr_str_vec = split(chr_str2, "_");

	// select chromosomes by 'chr'
	for(i=0; i<chr_str_vec.size(); i++){
		idx_chr = -1;
		for(j=0; j<chrname_vec.size(); j++){
			chrname = chrname_vec.at(j);
			if(chrname.compare(chr_str_vec.at(i))==0){
				idx_chr = j;
				break;
			}
		}
		if(idx_chr!=-1){
			chrname = chrname_vec.at(idx_chr);
			chrname_vec_sorted.push_back(chrname);
			selected_flag_array[idx_chr] = 1;
		}
	}

	// select chromosomes by 'chrA_*'
	for(i=0; i<chr_str_vec.size(); i++){
		head_str = chr_str_vec.at(i) + "_";
		for(j=0; j<chrname_vec.size(); j++){
			chrname = chrname_vec.at(j);
			head_str_chr = chrname.substr(0, head_str.size());
			if(head_str_chr.compare(head_str)==0){
				chrname_vec_sorted.push_back(chrname);
				selected_flag_array[j] = 1;
			}
		}
	}

	// add unselected items
	for(i=0; i<chrname_vec.size(); i++){
		if(selected_flag_array[i]==0){
			chrname = chrname_vec.at(i);
			chrname_vec_sorted.push_back(chrname);
		}
	}

	vector<string>().swap(chrname_vec);	// free memory

	free(selected_flag_array);

	return chrname_vec_sorted;
}

vector<vector<SV_item*>> constructSubsetByChrOp(vector<SV_item*> &user_data, vector<SV_item*> &benchmark_data, vector<string> &chrname_vec){
	vector< vector<SV_item*> > result;
	string chr;
	vector<SV_item*> subset1, subset2;

	for(size_t i=0; i<chrname_vec.size(); i++){
		chr = chrname_vec.at(i);
		subset1 = getItemsByChr(chr, user_data);
		subset2 = getItemsByChr(chr, benchmark_data);
		result.push_back(subset1);
		result.push_back(subset2);
	}

	return result;
}

vector<SV_item*> getItemsByChr(string &chrname, vector<SV_item*> &dataset){
	vector<SV_item*> result;
	size_t i;
	SV_item *item;

	for(i=0; i<dataset.size(); i++){
		item = dataset.at(i);
		if(IsSameChrname(item->chrname, chrname)) result.push_back(item);
	}

	return result;
}

void sortSubsets(vector<vector<SV_item*>> &subsets){
	sortWork_opt *sort_opt;
	size_t i;

	hts_tpool *p = hts_tpool_init(num_threads);
	hts_tpool_process *q = hts_tpool_process_init(p, num_threads*2, 1);

	pthread_mutex_init(&mtx_overlap, NULL);

	cout << "Sort SV items ..." << endl;

	for(i=0; i<subsets.size(); i++){
		sort_opt =  new sortWork_opt();
		sort_opt->dataset = &(subsets.at(i));

		//if(sort_opt->dataset->size()>0) cout << "\t" << sort_opt->dataset->at(0)->chrname << ", data size: " << sort_opt->dataset->size() << endl;
		hts_tpool_dispatch(p, q, sortSubsetOp, sort_opt);
	}

    hts_tpool_process_flush(q);
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);
}

void* sortSubsetOp(void *arg){
	sortWork_opt *sort_opt = (sortWork_opt *)arg;
	sort(sort_opt->dataset->begin(), sort_opt->dataset->end(), sortFunSameChr);
	delete sort_opt;
	return NULL;
}

// sort ascending
bool sortFunSameChr(const SV_item *item1, const SV_item *item2){
	return item1->startPos < item2->startPos;
}

// check order
void checkOrder(vector<vector<SV_item*>> &subsets){
	vector<SV_item*> subset;
	size_t i, j;
	SV_item *item1, *item2;

	for(i=0; i<subsets.size(); i++){
		subset = subsets.at(i);
		for(j=1; j<subset.size(); j++){
			item1 = subset.at(j-1);
			item2 = subset.at(j);
			if(item1->chrname.compare(item2->chrname)!=0){
				cerr << "different chrname, error!" << endl;
				exit(1);
			}else if(item1->startPos>item2->startPos){
				cerr << "item1.startPos=" << item1->startPos << " > item2.startPos=" << item2->startPos << ", error!" << endl;
				exit(1);
			}
		}
	}

	//cout << "Check sort finished." << endl;
}

vector<vector<SV_item*>> intersectOp(vector<vector<SV_item*>> &subsets, faidx_t *fai){
	vector<vector<SV_item*>> result;
	size_t i, j, subset_num = subsets.size() >> 1;
	vector<SV_item*> vec_tmp, intersect_vec_user[subset_num+1], intersect_vec_benchmark[subset_num+1], private_vec_user[subset_num+1], private_vec_benchmark[subset_num+1], intersect_vec_latentpositive_user[subset_num+1];
	overlapWork_opt *overlap_opt;

	hts_tpool *p = hts_tpool_init(num_threads);
	hts_tpool_process *q = hts_tpool_process_init(p, num_threads*2, 1);

	pthread_mutex_init(&mtx_overlap, NULL);
	pthread_mutex_init(&mutex_mem, NULL);

	cout << "Computing intersect information between user data and benchmark data ..." << endl;

	for(i=0; i<subset_num; i++){
		overlap_opt = new overlapWork_opt();
		overlap_opt->subset1 = subsets.at(i*2);
		overlap_opt->subset2 = subsets.at(i*2+1);
		overlap_opt->intersect_vec_user = intersect_vec_user + i;
		overlap_opt->intersect_vec_benchmark = intersect_vec_benchmark + i;
		overlap_opt->private_vec_user = private_vec_user + i;
		overlap_opt->private_vec_benchmark = private_vec_benchmark + i;
		overlap_opt->intersect_vec_latentpositive_user = intersect_vec_latentpositive_user + i;
		overlap_opt->fai = fai;

//		if(overlap_opt->subset1.size()>0) cout << "\t" << overlap_opt->subset1.at(0)->chrname << ", user_data: " << overlap_opt->subset1.size() << ", benchmark_data: " << overlap_opt->subset2.size() << endl;
		hts_tpool_dispatch(p, q, intersectSubset, overlap_opt);
	}

    hts_tpool_process_flush(q);
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);

    // merge sub-results
    for(i=0; i<subset_num; i++){
//    	if(subsets.at(i*2).size()>0) cout << "\t" << subsets.at(i*2).at(0)->chrname << ", user_data: " << subsets.at(i*2).size() << ", benchmark_data: " << subsets.at(i*2+1).size() << endl;
//    	cout << "\t\t" << intersect_vec_user[i].size() << ", " << intersect_vec_benchmark[i].size() << ", " << private_vec_user[i].size() << ", " << private_vec_benchmark[i].size() << endl;

    	for(j=0; j<intersect_vec_user[i].size(); j++) intersect_vec_user[subset_num].push_back(intersect_vec_user[i].at(j));
    	vector<SV_item*>().swap(intersect_vec_user[i]);
    	for(j=0; j<intersect_vec_benchmark[i].size(); j++) intersect_vec_benchmark[subset_num].push_back(intersect_vec_benchmark[i].at(j));
    	vector<SV_item*>().swap(intersect_vec_benchmark[i]);
    	for(j=0; j<private_vec_user[i].size(); j++) private_vec_user[subset_num].push_back(private_vec_user[i].at(j));
    	vector<SV_item*>().swap(private_vec_user[i]);
    	for(j=0; j<private_vec_benchmark[i].size(); j++) private_vec_benchmark[subset_num].push_back(private_vec_benchmark[i].at(j));
    	vector<SV_item*>().swap(private_vec_benchmark[i]);
    	for(j=0; j<intersect_vec_latentpositive_user[i].size(); j++) intersect_vec_latentpositive_user[subset_num].push_back(intersect_vec_latentpositive_user[i].at(j));
    	vector<SV_item*>().swap(intersect_vec_latentpositive_user[i]);
    }

	result.push_back(intersect_vec_user[subset_num]);
	result.push_back(intersect_vec_benchmark[subset_num]);
	result.push_back(private_vec_user[subset_num]);
	result.push_back(private_vec_benchmark[subset_num]);
	result.push_back(intersect_vec_latentpositive_user[subset_num]);

	return result;
}

void* intersectSubset(void *arg){
	overlapWork_opt *overlap_opt = (overlapWork_opt *)arg;
	vector<SV_item*> intersect_vec_user, intersect_vec_benchmark, private_vec_user, private_vec_benchmark, intersect_vec_latentpositive_user;
	SV_item *item1, *item2, *item, *item3;
	size_t i, j, k;
	vector<size_t> overlap_type_vec;
	vector<SV_item*> subset1, subset2;
	int64_t h, start_idx, end_idx, new_start_idx, start_search_pos, end_search_pos, extendsize_tmp, startpos_tmp, endpos_tmp;
	int64_t start_idx_tmp0, end_idx_tmp0, start_idx_tmp, end_idx_tmp, svlen_1, svlen_2, sv_len1, sv_len2;
	float svlenRatio_tmp, svlenRatio_tmp1;
	bool flag;
	double consistency;
	string seq1, seq2, aln_seq1, aln_seq2;

	//if(overlap_opt->subset1.size()>0) cout << "\t" << overlap_opt->subset1.at(0)->chrname << ", user_data: " << overlap_opt->subset1.size() << ", benchmark_data: " << overlap_opt->subset2.size() << endl;

	subset1 = overlap_opt->subset1;
	subset2 = overlap_opt->subset2;
	for(i=0; i<subset1.size(); i++) {subset1.at(i)->overlapped = false; subset1.at(i)->seqcons = "-"; subset1.at(i)->mergeMark = false; subset1.at(i)-> LpFlag = false;}
	for(i=0; i<subset2.size(); i++) {subset2.at(i)->overlapped = false;	subset2.at(i)->seqcons = "-"; subset2.at(i)->mergeMark = false; subset2.at(i)-> LpFlag = false;}

	start_idx = end_idx = -1;
	for(i=0; i<subset1.size(); i++){

		item1 = subset1.at(i);

//		cout << "i=" << i << ", " << item1->chrname << ":" << item1->startPos << "-" << item1->endPos << ", svtype=" << item1->sv_type << ", svlen=" << item1->sv_len << endl;

		start_idx_tmp0 = start_idx;
		end_idx_tmp0 = end_idx;

		for(k=0; k<2; k++){

			if(k==0) extendsize_tmp = 2 * extendSize;	//2 *
			else{
				svlen_1 = item1->sv_len;
				if(svlen_1<0) svlen_1 = -svlen_1;
				if(svlen_1>=minSizeLargeSV)
					extendsize_tmp = extendSizeLargeSV;
				else
					continue;
			}

			//compute limit search regions
			start_search_pos = item1->startPos - extendsize_tmp;
			end_search_pos = item1->endPos + extendsize_tmp;
			if(start_search_pos<1) start_search_pos = 1;

			//compute start element index
			if(start_idx_tmp0==-1) start_idx = 0;
			else{
				new_start_idx = -1;
				if(k==0){
					for(j=start_idx_tmp0; j<subset2.size(); j++){
						item2 = subset2.at(j);

						overlap_type_vec = computeOverlapType(item1, item2);
						if(overlap_type_vec.size()>0 and overlap_type_vec.at(0)!=NO_OVERLAP){
							new_start_idx = j;
							break;
						}else if(item2->endPos>=start_search_pos){
							new_start_idx = j;
							break;
						}
					}
				}else{
					if(start_idx_tmp0>=0 and start_idx_tmp0<(int64_t)subset2.size()){
						for(h=start_idx_tmp0; h>=0; h--){
							item2 = subset2.at(h);
							startpos_tmp = item2->startPos - extendSize;
							endpos_tmp = item2->endPos + extendSize;
							if(startpos_tmp<1) startpos_tmp = 1;
							flag = isOverlappedPos(startpos_tmp, endpos_tmp, start_search_pos, end_search_pos);
							if(flag){
								new_start_idx = h;
							}else
								break;
						}
					}
				}

				if(new_start_idx!=-1)
					start_idx = new_start_idx;
				else
					start_idx = end_idx_tmp0;
			}

			//begin overlap
			if(start_idx!=-1){ // valid start index
				end_idx = -1;
				if(k==0){
					for(j=start_idx; j<subset2.size(); j++){
						item2 = subset2.at(j);
						if(item2->startPos<start_search_pos) continue;
						overlap_type_vec = computeOverlapType(item1, item2);
						if(overlap_type_vec.size()>0 and overlap_type_vec.at(0)!=NO_OVERLAP){
							//merge_judge
							if(item1->sv_type == item2->sv_type)
								IsmergeJudge(i, item1, item2, subset1, overlap_opt->fai);

							sv_len1 = item1->sv_len; sv_len2 = item2->sv_len;
							if(sv_len1<0) sv_len1 = -sv_len1;
							if(sv_len2<0) sv_len2 = -sv_len2;

							if(sv_len1<=sv_len2) svlenRatio_tmp1 = (float)sv_len1 / sv_len2;
							else svlenRatio_tmp1 = (float)sv_len2 / sv_len1;

							flag = false;
							bool alleleflag = false;
							if(svlenRatio_tmp1>=svlenRatio){
								if(item1->sv_type==item2->sv_type or (item1->sv_type==VAR_INS and item2->sv_type==VAR_DUP) or (item1->sv_type==VAR_DUP and item2->sv_type==VAR_INS)){
									if(i!=subset1.size()-1 and j!=subset2.size()-1){

										if((subset1.at(i+1)->startPos - item1->startPos <= ALLELIC_DISTANCE and subset1.at(i+1)->sv_type == item1->sv_type and !item1->mergeMark) or (i>0 and item1->startPos - subset1.at(i-1)->startPos <= ALLELIC_DISTANCE and subset1.at(i-1)->sv_type == item1->sv_type and !item1->mergeMark))	alleleflag = true;

										if((subset2.at(j+1)->startPos - item2->startPos <= ALLELIC_DISTANCE and subset2.at(j+1)->sv_type == item2->sv_type  and !item1->mergeMark) or (j>0 and item2->startPos - subset2.at(j-1)->startPos <= ALLELIC_DISTANCE and subset2.at(j-1)->sv_type == item2->sv_type and !item1->mergeMark))	alleleflag = true; //or item2->startPos - subset2.at(j-1)->startPos <=10
									}else{
										if(i>0 and item1->startPos - subset1.at(i-1)->startPos <= ALLELIC_DISTANCE and subset1.at(i-1)->sv_type == item1->sv_type and !item1->mergeMark)	alleleflag = true;
										if(j>0 and item2->startPos - subset2.at(j-1)->startPos <= ALLELIC_DISTANCE and subset2.at(j-1)->sv_type == item2->sv_type and !item1->mergeMark)	alleleflag = true;
									}
									if(alleleflag){
										if(((item1->sv_type == VAR_INS and item2->sv_type == VAR_INS) or (item1->sv_type == VAR_DUP and item2->sv_type== VAR_DUP) or (item1->sv_type == VAR_INV and item2->sv_type== VAR_INV)) and typeMatchLevel == MATCHLEVEL_S){ //or (item1->sv_type == VAR_INS and item2->sv_type== VAR_DUP) or (item1->sv_type == VAR_DUP and item2->sv_type== VAR_INS)
//											if (item1->alt_seq.compare("-")==0 or item2->alt_seq.compare("-")==0 or item1->alt_seq.length()<(size_t)(sv_len1-1) or item1->alt_seq.length()>(size_t)(sv_len1+1)
//												or item2->alt_seq.length()<(size_t)(sv_len2-1) or item2->alt_seq.length()>(size_t)(sv_len2+1)) {
											if (item1->alt_seq.compare("-")==0 or item2->alt_seq.compare("-")==0) {
												flag = true;
											}else{
												consistency = computeVarseqConsistency(item1, item2, overlap_opt->fai, percentAlleleSeqIdentity);
												if(item1->seqcons.compare("-") == 0 and item2->seqcons.compare("-") == 0){
													item1->seqcons = item2->seqcons = to_string(consistency);
													if (consistency >= percentAlleleSeqIdentity) {
														SeqConsNumStat(consistency);
														flag = true;
													}
												}else if(item2->seqcons.compare("-") != 0 and stod(item2->seqcons) >= percentAlleleSeqIdentity){
													if(stod(item2->seqcons) < consistency){
														SeqconsSum -= stod(item2->seqcons);	SeqconsNum = SeqconsNum - 1;
														if(item1->seqcons.compare("-") == 0)
															item1->seqcons = item2->seqcons = to_string(consistency);
														else{
															if(consistency >= stod(item1->seqcons))	item1->seqcons = to_string(consistency);
															item2->seqcons = to_string(consistency);
														}
														if (consistency >= percentAlleleSeqIdentity) {
															SeqConsNumStat(consistency);
															flag = true;
														}
													}else{
														if (consistency >= percentAlleleSeqIdentity) {
															if(item1->seqcons.compare("-") ==0)
																item1->seqcons = to_string(consistency);
															else{
																if(consistency >= stod(item1->seqcons))	item1->seqcons = to_string(consistency);
															}
															flag = true;
														}
													}
												}else{
													if(item1->seqcons.compare("-") == 0)
														item1->seqcons = item2->seqcons = to_string(consistency);
													else{
														if(consistency >= stod(item1->seqcons))	item1->seqcons = to_string(consistency);
														item2->seqcons = to_string(consistency);
													}
													if (consistency >= percentAlleleSeqIdentity) {
														flag = true;
													}
												}
											}
										}else if(((item1->sv_type == VAR_INS and item2->sv_type == VAR_INS) or (item1->sv_type == VAR_DUP and item2->sv_type== VAR_DUP) or (item1->sv_type == VAR_INV and item2->sv_type== VAR_INV) or (item1->sv_type == VAR_INS and item2->sv_type== VAR_DUP) or (item1->sv_type == VAR_DUP and item2->sv_type== VAR_INS)) and typeMatchLevel == MATCHLEVEL_L){
											if (item1->alt_seq.compare("-")==0 or item2->alt_seq.compare("-")==0) {
													flag = true;
												}else{
													consistency = computeVarseqConsistency(item1, item2, overlap_opt->fai, percentAlleleSeqIdentity);
//													item1->seqcons = item2->seqcons = to_string(consistency);
//													if (consistency >= percentSeqIdentity) {
//														SeqConsNumStat(consistency);
//														flag = true;
//													}
													if(item1->seqcons.compare("-") == 0 and item2->seqcons.compare("-") == 0){
														item1->seqcons = item2->seqcons = to_string(consistency);
														if (consistency >= percentAlleleSeqIdentity) {
															SeqConsNumStat(consistency);
															flag = true;
														}
													}else if(item2->seqcons.compare("-") != 0 and stod(item2->seqcons) >= percentAlleleSeqIdentity){
														if(stod(item2->seqcons) < consistency){
															SeqconsSum -= stod(item2->seqcons);	SeqconsNum = SeqconsNum - 1;
//															item1->seqcons = item2->seqcons = to_string(consistency);
															if(item1->seqcons.compare("-") == 0)
																item1->seqcons = item2->seqcons = to_string(consistency);
															else{
																if(consistency >= stod(item1->seqcons))	item1->seqcons = to_string(consistency);
																item2->seqcons = to_string(consistency);
															}
															if (consistency >= percentAlleleSeqIdentity) {
																SeqConsNumStat(consistency);
																flag = true;
															}
														}else{
															if (consistency >= percentAlleleSeqIdentity) {
//																item1->seqcons = to_string(consistency);
																if(item1->seqcons.compare("-") ==0)
																	item1->seqcons = to_string(consistency);
																else{
																	if(consistency >= stod(item1->seqcons))	item1->seqcons = to_string(consistency);
																}
																flag = true;
															}
														}
													}else{
//														item1->seqcons = item2->seqcons = to_string(consistency);
														if(item1->seqcons.compare("-") == 0)
															item1->seqcons = item2->seqcons = to_string(consistency);
														else{
															if(consistency >= stod(item1->seqcons))	item1->seqcons = to_string(consistency);
															item2->seqcons = to_string(consistency);
														}
														if (consistency >= percentAlleleSeqIdentity) {
//															item1->seqcons = item2->seqcons = to_string(consistency);
															flag = true;
														}
													}
												}
										}else if(item1->sv_type == VAR_DEL and item2->sv_type == VAR_DEL){
											if (item1->ref_seq.compare("-") == 0 or item2->ref_seq.compare("-")== 0) {
												flag = true;
											}else{
												consistency = computeVarseqConsistency(item1, item2, overlap_opt->fai, percentAlleleSeqIdentity);
//												item1->seqcons = item2->seqcons = to_string(consistency);
//												if (consistency >= percentSeqIdentity) {
//													SeqConsNumStat(consistency);
//													flag = true;
//												}
												if(item1->seqcons.compare("-") == 0 and item2->seqcons.compare("-") == 0){
													item1->seqcons = item2->seqcons = to_string(consistency);
													if (consistency >= percentAlleleSeqIdentity) {
														SeqConsNumStat(consistency);
														flag = true;
													}
												}else if(item2->seqcons.compare("-") != 0 and stod(item2->seqcons) >= percentAlleleSeqIdentity){
													if(stod(item2->seqcons) < consistency){
														SeqconsSum -= stod(item2->seqcons);	SeqconsNum = SeqconsNum - 1;
//														item1->seqcons = item2->seqcons = to_string(consistency);
														if(item1->seqcons.compare("-") == 0)
															item1->seqcons = item2->seqcons = to_string(consistency);
														else{
															if(consistency >= stod(item1->seqcons))	item1->seqcons = to_string(consistency);
															item2->seqcons = to_string(consistency);
														}
														if (consistency >= percentAlleleSeqIdentity) {
															SeqConsNumStat(consistency);
															flag = true;
														}
													}else{
														if (consistency >= percentAlleleSeqIdentity) {
//															item1->seqcons = to_string(consistency);
															if(item1->seqcons.compare("-") ==0)
																item1->seqcons = to_string(consistency);
															else{
																if(consistency >= stod(item1->seqcons))	item1->seqcons = to_string(consistency);
															}
															flag = true;
														}
													}
												}else{
													item1->seqcons = item2->seqcons = to_string(consistency);
													if (consistency >= percentAlleleSeqIdentity) {
//														item1->seqcons = item2->seqcons = to_string(consistency);
														flag = true;
													}
												}
											}
										}else if((item1->sv_type == VAR_TRA and item2->sv_type == VAR_TRA) or(item1->sv_type == VAR_BND and item2->sv_type == VAR_BND) or(item1->sv_type == VAR_INV_TRA and item2->sv_type == VAR_INV_TRA)){	//TRA,BND,TRA_BND
											flag = true;
										}
									}else{
										if(((item1->sv_type == VAR_INS and item2->sv_type == VAR_INS) or (item1->sv_type == VAR_DUP and item2->sv_type== VAR_DUP) or (item1->sv_type == VAR_INV and item2->sv_type== VAR_INV)) and typeMatchLevel == MATCHLEVEL_S){ //or (item1->sv_type == VAR_INS and item2->sv_type== VAR_DUP) or (item1->sv_type == VAR_DUP and item2->sv_type== VAR_INS)
//											if (item1->alt_seq.compare("-")==0 or item2->alt_seq.compare("-")==0 or item1->alt_seq.length()<(size_t)(sv_len1-1) or item1->alt_seq.length()>(size_t)(sv_len1+1)
//												or item2->alt_seq.length()<(size_t)(sv_len2-1) or item2->alt_seq.length()>(size_t)(sv_len2+1)) {
											if (item1->alt_seq.compare("-")==0 or item2->alt_seq.compare("-")==0) {
												flag = true;
											}else{
												consistency = computeVarseqConsistency(item1, item2, overlap_opt->fai, percentSeqIdentity);
												if(item1->seqcons.compare("-") == 0 and item2->seqcons.compare("-") == 0){
													item1->seqcons = item2->seqcons = to_string(consistency);
													if (consistency >= percentSeqIdentity) {
														SeqConsNumStat(consistency);
														flag = true;
													}
												}else if(item2->seqcons.compare("-") != 0 and stod(item2->seqcons) >= percentSeqIdentity){
													if(stod(item2->seqcons) < consistency){
														SeqconsSum -= stod(item2->seqcons);	SeqconsNum = SeqconsNum - 1;
														if(item1->seqcons.compare("-") == 0)
															item1->seqcons = item2->seqcons = to_string(consistency);
														else{
															if(consistency >= stod(item1->seqcons))	item1->seqcons = to_string(consistency);
															item2->seqcons = to_string(consistency);
														}
														if (consistency >= percentSeqIdentity) {
															SeqConsNumStat(consistency);
															flag = true;
														}
													}else{
														if (consistency >= percentSeqIdentity) {
															if(item1->seqcons.compare("-") ==0)
																item1->seqcons = to_string(consistency);
															else{
																if(consistency >= stod(item1->seqcons))	item1->seqcons = to_string(consistency);
															}
															flag = true;
														}
													}
												}else{
													if(item1->seqcons.compare("-") == 0)
														item1->seqcons = item2->seqcons = to_string(consistency);
													else{
														if(consistency >= stod(item1->seqcons))	item1->seqcons = to_string(consistency);
														item2->seqcons = to_string(consistency);
													}
													if (consistency >= percentSeqIdentity) {
														flag = true;
													}
												}
											}
										}else if(((item1->sv_type == VAR_INS and item2->sv_type == VAR_INS) or (item1->sv_type == VAR_DUP and item2->sv_type== VAR_DUP) or (item1->sv_type == VAR_INV and item2->sv_type== VAR_INV) or (item1->sv_type == VAR_INS and item2->sv_type== VAR_DUP) or (item1->sv_type == VAR_DUP and item2->sv_type== VAR_INS)) and typeMatchLevel == MATCHLEVEL_L){
											if (item1->alt_seq.compare("-")==0 or item2->alt_seq.compare("-")==0) {
													flag = true;
												}else{
													consistency = computeVarseqConsistency(item1, item2, overlap_opt->fai, percentSeqIdentity);
//													item1->seqcons = item2->seqcons = to_string(consistency);
//													if (consistency >= percentSeqIdentity) {
//														SeqConsNumStat(consistency);
//														flag = true;
//													}
													if(item1->seqcons.compare("-") == 0 and item2->seqcons.compare("-") == 0){
														item1->seqcons = item2->seqcons = to_string(consistency);
														if (consistency >= percentSeqIdentity) {
															SeqConsNumStat(consistency);
															flag = true;
														}
													}else if(item2->seqcons.compare("-") != 0 and stod(item2->seqcons) >= percentSeqIdentity){
														if(stod(item2->seqcons) < consistency){
															SeqconsSum -= stod(item2->seqcons);	SeqconsNum = SeqconsNum - 1;
//															item1->seqcons = item2->seqcons = to_string(consistency);
															if(item1->seqcons.compare("-") == 0)
																item1->seqcons = item2->seqcons = to_string(consistency);
															else{
																if(consistency >= stod(item1->seqcons))	item1->seqcons = to_string(consistency);
																item2->seqcons = to_string(consistency);
															}
															if (consistency >= percentSeqIdentity) {
																SeqConsNumStat(consistency);
																flag = true;
															}
														}else{
															if (consistency >= percentSeqIdentity) {
//																item1->seqcons = to_string(consistency);
																if(item1->seqcons.compare("-") ==0)
																	item1->seqcons = to_string(consistency);
																else{
																	if(consistency >= stod(item1->seqcons))	item1->seqcons = to_string(consistency);
																}
																flag = true;
															}
														}
													}else{
//														item1->seqcons = item2->seqcons = to_string(consistency);
														if(item1->seqcons.compare("-") == 0)
															item1->seqcons = item2->seqcons = to_string(consistency);
														else{
															if(consistency >= stod(item1->seqcons))	item1->seqcons = to_string(consistency);
															item2->seqcons = to_string(consistency);
														}
														if (consistency >= percentSeqIdentity) {
//															item1->seqcons = item2->seqcons = to_string(consistency);
															flag = true;
														}
													}
												}
										}else if(item1->sv_type == VAR_DEL and item2->sv_type == VAR_DEL){
											if (item1->ref_seq.compare("-") == 0 or item2->ref_seq.compare("-")== 0) {
												flag = true;
											}else{
//												auto start = chrono::high_resolution_clock::now();
												consistency = computeVarseqConsistency(item1, item2, overlap_opt->fai, percentSeqIdentity);
//												auto end = chrono::high_resolution_clock::now();
//												chrono::duration<double> duration = end - start;
//												cout << "Execution time: " << duration.count() << " seconds" << endl;
//												item1->seqcons = item2->seqcons = to_string(consistency);
//												if (consistency >= percentSeqIdentity) {
//													SeqConsNumStat(consistency);
//													flag = true;
//												}
												if(item1->seqcons.compare("-") == 0 and item2->seqcons.compare("-") == 0){
													item1->seqcons = item2->seqcons = to_string(consistency);
													if (consistency >= percentSeqIdentity) {
														SeqConsNumStat(consistency);
														flag = true;
													}
												}else if(item2->seqcons.compare("-") != 0 and stod(item2->seqcons) >= percentSeqIdentity){
													if(stod(item2->seqcons) < consistency){
														SeqconsSum -= stod(item2->seqcons);	SeqconsNum = SeqconsNum - 1;
//														item1->seqcons = item2->seqcons = to_string(consistency);
														if(item1->seqcons.compare("-") == 0)
															item1->seqcons = item2->seqcons = to_string(consistency);
														else{
															if(consistency >= stod(item1->seqcons))	item1->seqcons = to_string(consistency);
															item2->seqcons = to_string(consistency);
														}
														if (consistency >= percentSeqIdentity) {
															SeqConsNumStat(consistency);
															flag = true;
														}
													}else{
														if (consistency >= percentSeqIdentity) {
//															item1->seqcons = to_string(consistency);
															if(item1->seqcons.compare("-") ==0)
																item1->seqcons = to_string(consistency);
															else{
																if(consistency >= stod(item1->seqcons))	item1->seqcons = to_string(consistency);
															}
															flag = true;
														}
													}
												}else{
													item1->seqcons = item2->seqcons = to_string(consistency);
													if (consistency >= percentSeqIdentity) {
//														item1->seqcons = item2->seqcons = to_string(consistency);
														flag = true;
													}
												}
											}
										}else if((item1->sv_type == VAR_TRA and item2->sv_type == VAR_TRA) or(item1->sv_type == VAR_BND and item2->sv_type == VAR_BND) or(item1->sv_type == VAR_INV_TRA and item2->sv_type == VAR_INV_TRA)){	//TRA,BND,TRA_BND
											flag = true;
										}
									}
								}
								if(flag){
									item1->overlapped = true;
									item2->overlapped = true;
									end_idx = j;
								}
							}
						}else if(item2->startPos>end_search_pos){
							end_idx = j - 1;
							if(end_idx<0) end_idx = 0;
							if(end_idx<start_idx) end_idx = start_idx;
							break;
						}
					}
				}else{ // k=1
					for(j=start_idx; j<subset2.size(); j++){
						item2 = subset2.at(j);
//						if(item2->startPos<start_search_pos) continue;
						if(item2->overlapped==false and item2->startPos<=end_search_pos and item2->startPos>=start_search_pos){
							svlen_2 = item2->sv_len;
							if(svlen_2<0) svlen_2 = -svlen_2;

							if(svlen_1<svlen_2) svlenRatio_tmp = (float)svlen_1 / svlen_2;
							else svlenRatio_tmp = (float)svlen_2 / svlen_1;

							flag = false;
							bool alleleflag1 = false;
							if(svlenRatio_tmp>=svlenRatio){
								if(item1->sv_type==item2->sv_type or (item1->sv_type==VAR_INS and item2->sv_type==VAR_DUP) or (item1->sv_type==VAR_DUP and item2->sv_type==VAR_INS)){ // or (item1->sv_type==VAR_INS and item2->sv_type==VAR_DUP) or (item1->sv_type==VAR_DUP and item2->sv_type==VAR_INS)
									if(i!=subset1.size()-1 and j!=subset2.size()-1){

											if((subset1.at(i+1)->startPos - item1->startPos <= ALLELIC_DISTANCE and subset1.at(i+1)->sv_type == item1->sv_type) or (i>0 and item1->startPos - subset1.at(i-1)->startPos <= ALLELIC_DISTANCE and subset1.at(i-1)->sv_type == item1->sv_type))	alleleflag1 = true; //or item1->startPos - subset1.at(i-1)->startPos <=10

											if((subset2.at(j+1)->startPos - item2->startPos <= ALLELIC_DISTANCE and subset2.at(j+1)->sv_type == item2->sv_type) or (j>0 and item2->startPos - subset2.at(j-1)->startPos <= ALLELIC_DISTANCE and subset2.at(j-1)->sv_type == item2->sv_type))	alleleflag1 = true;
									}else{
										if(i>0 and item1->startPos - subset1.at(i-1)->startPos <= ALLELIC_DISTANCE and subset1.at(i-1)->sv_type == item1->sv_type)	alleleflag1 = true;
										if(j>0 and item2->startPos - subset2.at(j-1)->startPos <= ALLELIC_DISTANCE and subset2.at(j-1)->sv_type == item2->sv_type)	alleleflag1 = true;
									}
									if(alleleflag1){
										// compute the sequence consistency
										if(((item1->sv_type==VAR_INS and item2->sv_type==VAR_INS) or (item1->sv_type==VAR_DUP and item2->sv_type==VAR_DUP) or (item1->sv_type==VAR_INV and item2->sv_type==VAR_INV)) and typeMatchLevel == MATCHLEVEL_S){ //(item1->sv_type==VAR_INS and item2->sv_type==VAR_DUP) or (item1->sv_type==VAR_DUP and item2->sv_type==VAR_INS) or
											if (item1->alt_seq.compare("-") == 0 or item2->alt_seq.compare("-") == 0){
												flag = true;
											}else{
												consistency = computeVarseqConsistency(item1, item2, overlap_opt->fai, percentAlleleSeqIdentity);
	//											item1->seqcons = item2->seqcons = to_string(consistency);
	//											// determine the overlap flag according to consistency
	//											if(consistency>=percentSeqIdentity) {
	//												SeqConsNumStat(consistency);
	//												flag = true;
	//											}
												if(item1->seqcons.compare("-") == 0 and item2->seqcons.compare("-") == 0){
													item1->seqcons = item2->seqcons = to_string(consistency);
													if (consistency >= percentAlleleSeqIdentity) {
														SeqConsNumStat(consistency);
														flag = true;
													}
												}else if(item2->seqcons.compare("-") != 0 and stod(item2->seqcons) >= percentAlleleSeqIdentity){
													if(stod(item2->seqcons) < consistency){
														SeqconsSum -= stod(item2->seqcons);	SeqconsNum = SeqconsNum - 1;
														item1->seqcons = item2->seqcons = to_string(consistency);
														if (consistency >= percentAlleleSeqIdentity) {
															SeqConsNumStat(consistency);
															flag = true;
														}
													}else{
														if (consistency >= percentAlleleSeqIdentity) {
															item1->seqcons = to_string(consistency);
															flag = true;
														}
													}
												}else{
													item1->seqcons = item2->seqcons = to_string(consistency);
													if (consistency >= percentAlleleSeqIdentity) {
														flag = true;
													}
												}
											}
										}else if(((item1->sv_type == VAR_INS and item2->sv_type == VAR_INS) or (item1->sv_type == VAR_DUP and item2->sv_type== VAR_DUP) or (item1->sv_type == VAR_INV and item2->sv_type== VAR_INV) or (item1->sv_type == VAR_INS and item2->sv_type== VAR_DUP) or (item1->sv_type == VAR_DUP and item2->sv_type== VAR_INS)) and typeMatchLevel == MATCHLEVEL_L){
											if (item1->alt_seq.compare("-")==0 or item2->alt_seq.compare("-")==0) {
													flag = true;
												}else{
													consistency = computeVarseqConsistency(item1, item2, overlap_opt->fai, percentAlleleSeqIdentity);
	//												item1->seqcons = item2->seqcons = to_string(consistency);
	//												if (consistency >= percentSeqIdentity) {
	//													SeqConsNumStat(consistency);
	//													flag = true;
	//												}
													if(item1->seqcons.compare("-") == 0 and item2->seqcons.compare("-") == 0){
														item1->seqcons = item2->seqcons = to_string(consistency);
														if (consistency >= percentSeqIdentity) {
															SeqConsNumStat(consistency);
															flag = true;
														}
													}else if(item2->seqcons.compare("-") != 0 and stod(item2->seqcons) >= percentAlleleSeqIdentity){
														if(stod(item2->seqcons) < consistency){
															SeqconsSum -= stod(item2->seqcons);	SeqconsNum = SeqconsNum - 1;
															item1->seqcons = item2->seqcons = to_string(consistency);
															if (consistency >= percentAlleleSeqIdentity) {
																SeqConsNumStat(consistency);
																flag = true;
															}
														}else{
															if (consistency >= percentAlleleSeqIdentity) {
																item1->seqcons = to_string(consistency);
																flag = true;
															}
														}
													}else{
														item1->seqcons = item2->seqcons = to_string(consistency);
														if (consistency >= percentAlleleSeqIdentity) {
	//														item1->seqcons = item2->seqcons = to_string(consistency);
															flag = true;
														}
													}
												}
										}else if(item1->sv_type==VAR_DEL and item2->sv_type==VAR_DEL){
											if (item1->ref_seq.compare("-") == 0 or item2->ref_seq.compare("-")== 0) {
												flag = true;
											}else{
												consistency = computeVarseqConsistency(item1, item2, overlap_opt->fai, percentAlleleSeqIdentity);
	//											item1->seqcons = item2->seqcons = to_string(consistency);
	//											// determine the overlap flag according to consistency
	//											if(consistency>=percentSeqIdentity) {
	//												SeqConsNumStat(consistency);
	//												flag = true;
	//											}
												if(item1->seqcons.compare("-") == 0 and item2->seqcons.compare("-") == 0){
													item1->seqcons = item2->seqcons = to_string(consistency);
													if (consistency >= percentAlleleSeqIdentity) {
														SeqConsNumStat(consistency);
														flag = true;
													}
												}else if(item2->seqcons.compare("-") != 0 and stod(item2->seqcons) >= percentAlleleSeqIdentity){
													if(stod(item2->seqcons) < consistency){
														SeqconsSum -= stod(item2->seqcons);	SeqconsNum = SeqconsNum - 1;
														item1->seqcons = item2->seqcons = to_string(consistency);
														if (consistency >= percentAlleleSeqIdentity) {
															SeqConsNumStat(consistency);
															flag = true;
														}
													}else{
														if (consistency >= percentAlleleSeqIdentity) {
															item1->seqcons = to_string(consistency);
															flag = true;
														}
													}
												}else{
													item1->seqcons = item2->seqcons = to_string(consistency);
													if (consistency >= percentAlleleSeqIdentity) {
	//													item1->seqcons = item2->seqcons = to_string(consistency);
														flag = true;
													}
												}
											}
										}else if((item1->sv_type == VAR_TRA and item2->sv_type == VAR_TRA) or(item1->sv_type == VAR_BND and item2->sv_type == VAR_BND) or(item1->sv_type == VAR_INV_TRA and item2->sv_type == VAR_INV_TRA)){
											flag = true;
										}
									}else{
										// compute the sequence consistency
										if(((item1->sv_type==VAR_INS and item2->sv_type==VAR_INS) or (item1->sv_type==VAR_DUP and item2->sv_type==VAR_DUP) or (item1->sv_type==VAR_INV and item2->sv_type==VAR_INV)) and typeMatchLevel == MATCHLEVEL_S){ //(item1->sv_type==VAR_INS and item2->sv_type==VAR_DUP) or (item1->sv_type==VAR_DUP and item2->sv_type==VAR_INS) or
											if (item1->alt_seq.compare("-") == 0 or item2->alt_seq.compare("-") == 0){
												flag = true;
											}else{
												consistency = computeVarseqConsistency(item1, item2, overlap_opt->fai, percentSeqIdentity);
		//											item1->seqcons = item2->seqcons = to_string(consistency);
		//											// determine the overlap flag according to consistency
		//											if(consistency>=percentSeqIdentity) {
		//												SeqConsNumStat(consistency);
		//												flag = true;
		//											}
												if(item1->seqcons.compare("-") == 0 and item2->seqcons.compare("-") == 0){
													item1->seqcons = item2->seqcons = to_string(consistency);
													if (consistency >= percentSeqIdentity) {
														SeqConsNumStat(consistency);
														flag = true;
													}
												}else if(item2->seqcons.compare("-") != 0 and stod(item2->seqcons) >= percentSeqIdentity){
													if(stod(item2->seqcons) < consistency){
														SeqconsSum -= stod(item2->seqcons);	SeqconsNum = SeqconsNum - 1;
														item1->seqcons = item2->seqcons = to_string(consistency);
														if (consistency >= percentSeqIdentity) {
															SeqConsNumStat(consistency);
															flag = true;
														}
													}else{
														if (consistency >= percentSeqIdentity) {
															item1->seqcons = to_string(consistency);
															flag = true;
														}
													}
												}else{
													item1->seqcons = item2->seqcons = to_string(consistency);
													if (consistency >= percentSeqIdentity) {
														flag = true;
													}
												}
											}
										}else if(((item1->sv_type == VAR_INS and item2->sv_type == VAR_INS) or (item1->sv_type == VAR_DUP and item2->sv_type== VAR_DUP) or (item1->sv_type == VAR_INV and item2->sv_type== VAR_INV) or (item1->sv_type == VAR_INS and item2->sv_type== VAR_DUP) or (item1->sv_type == VAR_DUP and item2->sv_type== VAR_INS)) and typeMatchLevel == MATCHLEVEL_L){
											if (item1->alt_seq.compare("-")==0 or item2->alt_seq.compare("-")==0) {
													flag = true;
												}else{
													consistency = computeVarseqConsistency(item1, item2, overlap_opt->fai, percentSeqIdentity);
		//												item1->seqcons = item2->seqcons = to_string(consistency);
		//												if (consistency >= percentSeqIdentity) {
		//													SeqConsNumStat(consistency);
		//													flag = true;
		//												}
													if(item1->seqcons.compare("-") == 0 and item2->seqcons.compare("-") == 0){
														item1->seqcons = item2->seqcons = to_string(consistency);
														if (consistency >= percentSeqIdentity) {
															SeqConsNumStat(consistency);
															flag = true;
														}
													}else if(item2->seqcons.compare("-") != 0 and stod(item2->seqcons) >= percentSeqIdentity){
														if(stod(item2->seqcons) < consistency){
															SeqconsSum -= stod(item2->seqcons);	SeqconsNum = SeqconsNum - 1;
															item1->seqcons = item2->seqcons = to_string(consistency);
															if (consistency >= percentSeqIdentity) {
																SeqConsNumStat(consistency);
																flag = true;
															}
														}else{
															if (consistency >= percentSeqIdentity) {
																item1->seqcons = to_string(consistency);
																flag = true;
															}
														}
													}else{
														item1->seqcons = item2->seqcons = to_string(consistency);
														if (consistency >= percentSeqIdentity) {
		//														item1->seqcons = item2->seqcons = to_string(consistency);
															flag = true;
														}
													}
												}
										}else if(item1->sv_type==VAR_DEL and item2->sv_type==VAR_DEL){
											if (item1->ref_seq.compare("-") == 0 or item2->ref_seq.compare("-")== 0) {
												flag = true;
											}else{
												consistency = computeVarseqConsistency(item1, item2, overlap_opt->fai, percentSeqIdentity);
		//											item1->seqcons = item2->seqcons = to_string(consistency);
		//											// determine the overlap flag according to consistency
		//											if(consistency>=percentSeqIdentity) {
		//												SeqConsNumStat(consistency);
		//												flag = true;
		//											}
												if(item1->seqcons.compare("-") == 0 and item2->seqcons.compare("-") == 0){
													item1->seqcons = item2->seqcons = to_string(consistency);
													if (consistency >= percentSeqIdentity) {
														SeqConsNumStat(consistency);
														flag = true;
													}
												}else if(item2->seqcons.compare("-") != 0 and stod(item2->seqcons) >= percentSeqIdentity){
													if(stod(item2->seqcons) < consistency){
														SeqconsSum -= stod(item2->seqcons);	SeqconsNum = SeqconsNum - 1;
														item1->seqcons = item2->seqcons = to_string(consistency);
														if (consistency >= percentSeqIdentity) {
															SeqConsNumStat(consistency);
															flag = true;
														}
													}else{
														if (consistency >= percentSeqIdentity) {
															item1->seqcons = to_string(consistency);
															flag = true;
														}
													}
												}else{
													item1->seqcons = item2->seqcons = to_string(consistency);
													if (consistency >= percentSeqIdentity) {
		//													item1->seqcons = item2->seqcons = to_string(consistency);
														flag = true;
													}
												}
											}
										}else if((item1->sv_type == VAR_TRA and item2->sv_type == VAR_TRA) or(item1->sv_type == VAR_BND and item2->sv_type == VAR_BND) or(item1->sv_type == VAR_INV_TRA and item2->sv_type == VAR_INV_TRA)){
											flag = true;
										}
									}
							}
								if(flag){
									item1->overlapped = true;
									item2->overlapped = true;
									end_idx = j;
								}
							}
						}

						if(item2->startPos>end_search_pos){
							end_idx = j - 1;
							if(end_idx<0) end_idx = 0;
							if(end_idx<start_idx) end_idx = start_idx;
							break;
						}
	 				}
				}

				if(k==0){
					start_idx_tmp = start_idx;
					end_idx_tmp = end_idx;
				}else{
					start_idx = start_idx_tmp;
					end_idx = end_idx_tmp;
				}
			}
		}
	}

	for(i=0; i<subset2.size(); i++){
		item2 = subset2.at(i);
		item = itemdup(item2);
		if(item2->overlapped) intersect_vec_benchmark.push_back(item);
		else{ private_vec_benchmark.push_back(item);
			auto it = benchmarklineMap.find(item->lineInfo);
			    if (it != benchmarklineMap.end()) {
			        ++it->second;
			    } else {
			        cout << "Record not found: " << item->lineInfo << std::endl;
			    }
		}
	}

	for(i=0; i<subset1.size(); i++){
		item1 = subset1.at(i);
		item = itemdup(item1);
		/*if(!item1->overlapped) private_vec_user.push_back(item);
		else intersect_vec_user.push_back(item);*/
		if (item1->overlapped) {
			if(!item1->mergeMark)
				intersect_vec_user.push_back(item);
			else{
				item3 = itemdup(item1);
				intersect_vec_user.push_back(item3);
				intersect_vec_latentpositive_user.push_back(item);
			}
		} else {
			if (item1->mergeMark) {
				intersect_vec_latentpositive_user.push_back(item); //intersect_vec_latentpositive_user
			} else {
				private_vec_user.push_back(item);
			}
		}
	}

	pthread_mutex_lock(&mtx_overlap);
	for(j=0; j<intersect_vec_user.size(); j++) overlap_opt->intersect_vec_user->push_back(intersect_vec_user.at(j));
	for(j=0; j<intersect_vec_benchmark.size(); j++) overlap_opt->intersect_vec_benchmark->push_back(intersect_vec_benchmark.at(j));
	for(j=0; j<private_vec_user.size(); j++) overlap_opt->private_vec_user->push_back(private_vec_user.at(j));
	for(j=0; j<private_vec_benchmark.size(); j++) overlap_opt->private_vec_benchmark->push_back(private_vec_benchmark.at(j));
	for(j=0; j<intersect_vec_latentpositive_user.size(); j++) overlap_opt->intersect_vec_latentpositive_user->push_back(intersect_vec_latentpositive_user.at(j));
	//if(subset1.size()>0) cout << "\t" << subset1.at(0)->chrname << ", user_data: " << subset1.size() << ", benchmark_data: " << subset2.size() << endl;
	//cout << "\t\t" << intersect_vec_user.size() << ", " << intersect_vec_benchmark.size() << ", " << private_vec_user.size() << ", " << private_vec_benchmark.size() << endl;
	pthread_mutex_unlock(&mtx_overlap);

	delete overlap_opt;

	return NULL;
}

// upper the sequence
void upperSeq(string &seq){
	for(size_t i=0; i<seq.size(); i++)
		if(seq[i]>='a' and seq[i]<='z') seq[i] -= 32;
}

SV_item* itemdup(SV_item* item){
	SV_item *item_new = new SV_item();
	item_new->chrname = item->chrname;
	item_new->startPos = item->startPos;
	item_new->endPos = item->endPos;
	item_new->chrname2 = item->chrname2;
	item_new->startPos2 = item->startPos2;
	item_new->endPos2 = item->endPos2;
	item_new->sv_type = item->sv_type;
	item_new->sv_len = item->sv_len;
	item_new->overlapped = item->overlapped;
	item_new->mergeMark = item->mergeMark;
	item_new->LpFlag = item->LpFlag;
	item_new->validFlag = item->validFlag;
	for(size_t i=0; i<4; i++) item_new->traOverlappedArr[i] = item->traOverlappedArr[i];
	item_new->ref_seq = item->ref_seq;
	item_new->alt_seq = item->alt_seq;
	item_new->seqcons = item->seqcons;
	item_new->lineInfo = item->lineInfo;
	return item_new;
}

//judge chrname is same
bool IsSameChrname(string &chrname1, string &chrname2){
	string chr = "chr";
	if(chrname1.compare(chrname2)==0 or (chr + chrname1).compare(chrname2)==0 or chrname1.compare(chr + chrname2)==0){
		return true;
	}
	return false;
}

vector<size_t> computeOverlapType(SV_item* item1, SV_item* item2){
	bool flag = false, flag_tmp = false;
	int32_t startPos1, endPos1, startPos2, endPos2;
	vector<size_t> overlap_type_vec;
	size_t overlap_type;

	if(item1->sv_type!=VAR_TRA and item1->sv_type!=VAR_BND){ // INS, DEL, DUP, INV
		if(item1->chrname.size()>0 and IsSameChrname(item1->chrname, item2->chrname)){
			startPos1 = item1->startPos - extendSize;
			if(startPos1<1) startPos1 = 1;
			endPos1 = item1->endPos + extendSize;
			startPos2 = item2->startPos - extendSize;
			if(startPos2<1) startPos2 = 1;
			endPos2 = item2->endPos + extendSize;

			flag = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
		}

		if(flag==false and (item2->sv_type==VAR_TRA or item2->sv_type==VAR_BND)){
			if(item1->chrname.size()>0 and IsSameChrname(item1->chrname, item2->chrname2)){
				startPos1 = item1->startPos - extendSize;
				if(startPos1<1) startPos1 = 1;
				endPos1 = item1->endPos + extendSize;
				startPos2 = item2->startPos2 - extendSize;
				if(startPos2<1) startPos2 = 1;
				endPos2 = item2->endPos2 + extendSize;

				flag = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
			}
		}

		if(flag) overlap_type_vec.push_back(WHOLE_REG_OVERLAP);

	}else{ // TRA or BND

		if(item1->chrname.size()>0){
			if(IsSameChrname(item1->chrname, item2->chrname)){
				if(item1->startPos!=0){
					// 1: one_overlap_one_bp
					startPos1 = item1->startPos - extendSize;
					if(startPos1<1) startPos1 = 1;
					endPos1 = item1->startPos + extendSize;

					if(item2->startPos!=0){
						startPos2 = item2->startPos - extendSize;
						if(startPos2<1) startPos2 = 1;
						endPos2 = item2->startPos + extendSize;
						flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
					}

					if(flag_tmp) overlap_type = ONE_OVERLAP_ONE_BP;
					else if(item2->startPos!=item2->endPos and item2->endPos!=0){
						// one_overlap_two_bp
						startPos2 = item2->endPos - extendSize;
						if(startPos2<1) startPos2 = 1;
						endPos2 = item2->endPos + extendSize;

						flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
						if(flag_tmp) overlap_type = ONE_OVERLAP_TWO_BP;
					}
					if(flag_tmp) overlap_type_vec.push_back(overlap_type);
				}

				if(item1->startPos!=item1->endPos and item1->endPos!=0){
					// 2: two_overlap_one_bp
					startPos1 = item1->endPos - extendSize;
					if(startPos1<1) startPos1 = 1;
					endPos1 = item1->endPos + extendSize;

					if(item2->startPos!=0){
						startPos2 = item2->startPos - extendSize;
						if(startPos2<1) startPos2 = 1;
						endPos2 = item2->startPos + extendSize;
						flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
					}

					if(flag_tmp) overlap_type = TWO_OVERLAP_ONE_BP;
					else if(item2->startPos!=item2->endPos and item2->endPos!=0){
						// two_overlap_two_bp
						startPos2 = item2->endPos - extendSize;
						if(startPos2<1) startPos2 = 1;
						endPos2 = item2->endPos + extendSize;

						flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
						if(flag_tmp) overlap_type = TWO_OVERLAP_TWO_BP;
					}
					if(flag_tmp) overlap_type_vec.push_back(overlap_type);
				}
			}else if(IsSameChrname(item1->chrname, item2->chrname2)){

				if(item1->startPos!=0){
					// one_overlap_three_bp
					startPos1 = item1->startPos - extendSize;
					if(startPos1<1) startPos1 = 1;
					endPos1 = item1->startPos + extendSize;

					if(item2->startPos2!=0){
						startPos2 = item2->startPos2 - extendSize;
						if(startPos2<1) startPos2 = 1;
						endPos2 = item2->startPos2 + extendSize;
						flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
					}

					if(flag_tmp) overlap_type = ONE_OVERLAP_THREE_BP;
					else if(item2->startPos2!=item2->endPos2 and item2->endPos2!=0){
						// one_overlap_four_bp
						startPos2 = item2->endPos2 - extendSize;
						if(startPos2<1) startPos2 = 1;
						endPos2 = item2->endPos2 + extendSize;

						flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
						if(flag_tmp) overlap_type = ONE_OVERLAP_FOUR_BP;
					}
					if(flag_tmp) overlap_type_vec.push_back(overlap_type);
				}

				if(item1->startPos!=item1->endPos and item1->endPos!=0){
					// two_overlap_three_bp
					startPos1 = item1->endPos - extendSize;
					if(startPos1<1) startPos1 = 1;
					endPos1 = item1->endPos + extendSize;

					if(item2->startPos2!=0){
						startPos2 = item2->startPos2 - extendSize;
						if(startPos2<1) startPos2 = 1;
						endPos2 = item2->startPos2 + extendSize;
						flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
					}

					if(flag_tmp) overlap_type = TWO_OVERLAP_THREE_BP;
					else if(item2->startPos2!=item2->endPos2 and item2->endPos2!=0) { // two_overlap_four_bp
						startPos2 = item2->endPos2 - extendSize;
						if(startPos2<1) startPos2 = 1;
						endPos2 = item2->endPos2 + extendSize;

						flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
						if(flag_tmp) overlap_type = TWO_OVERLAP_FOUR_BP;
					}
					if(flag_tmp) overlap_type_vec.push_back(overlap_type);
				}
			}
		}

		if(item1->chrname2.size()>0){
			if(IsSameChrname(item1->chrname2, item2->chrname)){
				if(item1->startPos2!=0){
					// 3: three_overlap_one_bp
					startPos1 = item1->startPos2 - extendSize;
					if(startPos1<1) startPos1 = 1;
					endPos1 = item1->startPos2 + extendSize;

					if(item2->startPos!=0){
						startPos2 = item2->startPos - extendSize;
						if(startPos2<1) startPos2 = 1;
						endPos2 = item2->startPos + extendSize;
						flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
					}

					if(flag_tmp) overlap_type = THREE_OVERLAP_ONE_BP;
					else if(item2->startPos!=item2->endPos and item2->endPos!=0){
						// three_overlap_two_bp
						startPos2 = item2->endPos - extendSize;
						if(startPos2<1) startPos2 = 1;
						endPos2 = item2->endPos + extendSize;

						flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
						if(flag_tmp) overlap_type = THREE_OVERLAP_TWO_BP;
					}
					if(flag_tmp) overlap_type_vec.push_back(overlap_type);
				}

				if(item1->startPos2!=item1->endPos2 and item1->endPos2!=0){
					// 4: four_overlap_one_bp
					startPos1 = item1->endPos2 - extendSize;
					if(startPos1<1) startPos1 = 1;
					endPos1 = item1->endPos2 + extendSize;

					if(item2->startPos!=0){
						startPos2 = item2->startPos - extendSize;
						if(startPos2<1) startPos2 = 1;
						endPos2 = item2->startPos + extendSize;
						flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
					}

					if(flag_tmp) overlap_type = FOUR_OVERLAP_ONE_BP;
					else if(item2->startPos!=item2->endPos and item2->endPos!=0){
						// four_overlap_two_bp
						startPos2 = item2->endPos - extendSize;
						if(startPos2<1) startPos2 = 1;
						endPos2 = item2->endPos + extendSize;

						flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
						if(flag_tmp) overlap_type = FOUR_OVERLAP_TWO_BP;
					}
					if(flag_tmp) overlap_type_vec.push_back(overlap_type);
				}
			}else if(IsSameChrname(item1->chrname2, item2->chrname2)){
				if(item1->startPos2!=0){
					// three_overlap_three_bp
					startPos1 = item1->startPos2 - extendSize;
					if(startPos1<1) startPos1 = 1;
					endPos1 = item1->startPos2 + extendSize;

					if(item2->startPos2!=0){
						startPos2 = item2->startPos2 - extendSize;
						if(startPos2<1) startPos2 = 1;
						endPos2 = item2->startPos2 + extendSize;
						flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
					}

					if(flag_tmp) overlap_type = THREE_OVERLAP_THREE_BP;
					else if(item2->startPos2!=item2->endPos2 and item2->endPos2!=0) {
						// three_overlap_four_bp
						startPos2 = item2->endPos2 - extendSize;
						if(startPos2<1) startPos2 = 1;
						endPos2 = item2->endPos2 + extendSize;

						flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
						if(flag_tmp) overlap_type = THREE_OVERLAP_FOUR_BP;
					}
					if(flag_tmp) overlap_type_vec.push_back(overlap_type);
				}

				if(item1->startPos2!=item1->endPos2 and item1->endPos2!=0){
					// four_overlap_three_bp
					startPos1 = item1->endPos2 - extendSize;
					if(startPos1<1) startPos1 = 1;
					endPos1 = item1->endPos2 + extendSize;

					if(item2->startPos2!=0){
						startPos2 = item2->startPos2 - extendSize;
						if(startPos2<1) startPos2 = 1;
						endPos2 = item2->startPos2 + extendSize;
						flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
					}

					if(flag_tmp) overlap_type = FOUR_OVERLAP_THREE_BP;
					else if(item2->startPos2!=item2->endPos2 and item2->endPos2!=0) { // four_overlap_four_bp
						startPos2 = item2->endPos2 - extendSize;
						if(startPos2<1) startPos2 = 1;
						endPos2 = item2->endPos2 + extendSize;

						flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
						if(flag_tmp) overlap_type = FOUR_OVERLAP_FOUR_BP;
					}
					if(flag_tmp) overlap_type_vec.push_back(overlap_type);
				}
			}
		}
	}

	return overlap_type_vec;
}

vector<bam1_t*> load_region_data(const string& chrname, int start_pos, int end_pos, int min_mapq = 20) {
    vector<bam1_t*> records;
    static bool log_level_configured = false;
    if(!log_level_configured){
    	hts_set_log_level(HTS_LOG_ERROR);
    	log_level_configured = true;
    }
    // Open the BAM file
    samFile* bam_file = sam_open(BamFilePath.c_str(), "r");
    if (!bam_file) {
        cerr << "Error opening BAM file: " << BamFilePath << endl;
        return records;
    }

    // Load the BAM header
    bam_hdr_t* header = sam_hdr_read(bam_file);
    if (!header) {
        cerr << "Error reading BAM header" << endl;
        sam_close(bam_file);
        return records;
    }

    // Load the index
    hts_idx_t* idx = sam_index_load(bam_file, BamFilePath.c_str());
    if (!idx) {
        cerr << "Error loading BAM index" << endl;
        bam_hdr_destroy(header);
        sam_close(bam_file);
        return records;
    }

    // Construct the region string
    ostringstream region_stream;
    region_stream << chrname << ":" << start_pos << "-" << end_pos;
    string region = region_stream.str();

    // Create a region iterator
    hts_itr_t* iter = sam_itr_querys(idx, header, region.c_str());
    if (!iter) {
        cerr << "Error parsing region: " << region << endl;
        hts_idx_destroy(idx);
        bam_hdr_destroy(header);
        sam_close(bam_file);
        return records;
    }

    // Traverse the record
    bam1_t* bam_record = bam_init1();
    int record_count = 0;
    const int max_records = 50;

    while (sam_itr_next(bam_file, iter, bam_record) >= 0 && record_count < max_records) {
        if (bam_record->core.qual < min_mapq) continue;
        records.push_back(bam_dup1(bam_record));
        record_count++;
        if (record_count >= max_records) break;
    }

    // Clear resources
    bam_destroy1(bam_record);
    hts_itr_destroy(iter);
    hts_idx_destroy(idx);
    bam_hdr_destroy(header);
    sam_close(bam_file);

    return records;
}

void IsmergeJudge(size_t user_pos, SV_item* user, SV_item* bench, vector<SV_item*> subset1, faidx_t *fai){
	vector<SV_item*> merge_users_vec, merge_item_vec;
	vector<vector<SV_item*>> clusters;
	SV_item *Maymerge_user, *Merge_users, *Merge_item, *Lp_Merge_item, *Definite_Merge_item;
	size_t i = user_pos + 1;
	int32_t min_startpos, max_endpos, refseq_len_tmp;
	int64_t pos_update;
	bool flag = false;
	size_t startpos1, startpos2, endpos1, endpos2;
	string seq1="", seq2="", seq_extact="", reg_str, aln_seq1, aln_seq2;
	double consistency;
	char *seq;
	vector<bam1_t*> regionReads;
	int64_t sv_length;

	merge_users_vec.push_back(user);
	pos_update = user->endPos;

	if(i<subset1.size()){
		for(i = user_pos + 1; i<subset1.size(); i++){
			Maymerge_user = subset1.at(i);
			if(user->sv_type == VAR_INS and bench->sv_type == VAR_INS){
				if(abs(Maymerge_user->startPos - user->startPos) > 500) break;
				else{
					if(Maymerge_user->sv_type == user->sv_type)
						merge_users_vec.push_back(Maymerge_user);
				}
			}else if(user->sv_type == VAR_DEL and bench->sv_type == VAR_DEL){
				//distance
				if(abs(Maymerge_user->startPos - pos_update) > 400) break;
				else{
					if(Maymerge_user->sv_type == user->sv_type){
						merge_users_vec.push_back(Maymerge_user);
						pos_update = Maymerge_user->endPos;
					}
				}
			}else break;
		}
		//BAM
		if(BamFileSign and merge_users_vec.size()>1){
			//Area range calculation
			int32_t region_start = std::max(static_cast<int32_t>(0), static_cast<int32_t>(merge_users_vec.front()->startPos - 200));
			int32_t region_end = merge_users_vec.back()->endPos + 200;
			regionReads = load_region_data(user->chrname, region_start, region_end);
//			cout<<"111111111111:"<<regionReads.size()<<endl;
			 // Build SV location index (variation for fast query coverage)
			map<pair<int64_t, int64_t>, SV_item*> sv_pos_map;
			for (SV_item* sv : merge_users_vec) {
				if(sv->sv_type == VAR_INS)
					sv_pos_map[{sv->startPos, sv->endPos + sv->sv_len}] = sv;
				else
					sv_pos_map[{sv->startPos, sv->endPos}] = sv;
			}

			// Iterate through all reads and record the covered SV set
			unordered_map<string, vector<SV_item*>> read_sv_groups;
			for (bam1_t* read : regionReads) {
				// Gets the alignment position and direction of the read
//				size_t read_start = read->core.pos + 1;
//				size_t read_end = bam_endpos(read);

				// Analyze the CIGAR and calculate the coverage area (more accurate way)
				uint32_t* cigar = bam_get_cigar(read);
				int n_cigar = read->core.n_cigar;
				int64_t ref_pos = read->core.pos;
				int64_t read_len = 0;

				// Traverse the CIGAR operation to calculate the reference genome coverage interval
				vector<pair<int64_t, int64_t>> ref_blocks;
				for (int i = 0; i < n_cigar; ++i) {
				    int op = bam_cigar_op(cigar[i]);
				    int len = bam_cigar_oplen(cigar[i]);
				    switch(op) {
				        case BAM_CDEL:
				        	if(len >= 20){
				        		ref_blocks.push_back({ref_pos, ref_pos + len});
				        	}
				            ref_pos += len;
				            break;
				        case BAM_CINS:
				        	if(len >= 20){
				        		ref_blocks.push_back({ref_pos, ref_pos + len});
				        	}
				            read_len += len;
				            break;
				        case BAM_CMATCH:
				        case BAM_CEQUAL:
				        case BAM_CDIFF:
				        case BAM_CREF_SKIP:
				        case BAM_CSOFT_CLIP:
				        case BAM_CHARD_CLIP:
				            // Update reference coordinates or read lengths without recording intervals
				            if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
				                ref_pos += len;
				                read_len += len;
				            } else if (op == BAM_CREF_SKIP) {
				                ref_pos += len;
				            } else {
				                read_len += len;
				            }
				            break;
				    }
				}

				// Check whether each SV is overwritten by the current read
				vector<SV_item*> covered_svs;
				for (auto& sv_pair : sv_pos_map) {
				    int64_t sv_start = sv_pair.first.first;
				    int64_t sv_end = sv_pair.first.second;
				    SV_item* sv = sv_pair.second;
				    if(sv->sv_type == VAR_INS)
				    	sv_length = sv_end - sv_start;
				    else
				    	sv_length = sv_end - sv_start +1;
				    double ratio = 0;
				    // Calculate the total length covered
				    int64_t covered_length = 0;
				    for (auto& block : ref_blocks) {
				        //Intersection of computation interval
				        int64_t overlap_start = max(sv_start, block.first);
				        int64_t overlap_end = min(sv_end, block.second);
				        if (overlap_start - 50 < overlap_end + 50) {
				        	if(abs(sv_start-block.first) > 250) continue;
				        	if(overlap_start < overlap_end)
				        		covered_length += (overlap_end - overlap_start);
				        	else
				        		covered_length += (block.second - block.first);
				            if(sv_length >= (block.second - block.first)){
				            	ratio = (double)(block.second - block.first)/sv_length;
				            }else{
				            	ratio = (double)sv_length/(block.second - block.first);
				            }
				            if((double)covered_length/(block.second - block.first) >= 0 and ratio >= 0.9){	//  0.7/0.95
				            	covered_svs.push_back(sv);
				            }
				        }
				    }
				}
				// Records reads covering multiple SVS
				if (covered_svs.size() >= 2) {
				    string qname = bam_get_qname(read);
				    read_sv_groups[qname] = covered_svs;
				}
			}

			// Merge SVS that share the same read (build clusters)
			unordered_map<SV_item*, int> sv_to_cluster;
			clusters.clear();
			int cluster_id = 0;

			for (auto& group : read_sv_groups) {
				vector<SV_item*>& svs = group.second;
				vector<int> existing_clusters;

				// Checks whether the current SV already belongs to a cluster
				for (SV_item* sv : svs) {
					if (sv_to_cluster.find(sv) != sv_to_cluster.end()) {
						existing_clusters.push_back(sv_to_cluster[sv]);
					}
				}

				//Merge or create a new cluster
				if (existing_clusters.empty()) {
					for (SV_item* sv : svs) {
						sv_to_cluster[sv] = cluster_id;
					}
					clusters.push_back(svs);
					cluster_id++;
				} else {
					int target_cluster = existing_clusters[0];
					// duplicate removal
					for (int cid : existing_clusters) {
						if (cid != target_cluster) {
							for (SV_item* sv : clusters[cid]) {
								sv_to_cluster[sv] = target_cluster;
							}
							clusters[target_cluster].insert(
								clusters[target_cluster].end(),
								clusters[cid].begin(),
								clusters[cid].end()
							);
							clusters[cid].clear();
						}
					}
					// Add the current SV to the target cluster
					for (SV_item* sv : svs) {
						if (sv_to_cluster.find(sv) == sv_to_cluster.end()) {
							sv_to_cluster[sv] = target_cluster;
							clusters[target_cluster].push_back(sv);
						}
					}
				}
			}

			// Clean empty clustering
			vector<vector<SV_item*>> final_clusters;
			for (auto& cluster : clusters) {
				if (!cluster.empty()) {
					// Sort by startPos ascending order
					sort(cluster.begin(), cluster.end(), [](SV_item* a, SV_item* b) {
						return a->startPos < b->startPos;
					});
					final_clusters.push_back(cluster);
				}
			}
			clusters = final_clusters;

//			cout << "Found " << clusters.size() << " clusters supported by reads" << endl;

		}

		for (bam1_t* read : regionReads) {
		    bam_destroy1(read);
		}
		regionReads.clear();

		if(BamFileSign){
			for(size_t m =0; m < clusters.size(); m++){
				if(clusters[m].size()>1){
					//Merge verification: Whether it is consistent after the merger
					auto merge_result = findClosestSubarray(clusters[0], bench->sv_len);
					merge_item_vec = merge_result.second;
					if(merge_item_vec.size() > 1){
						//Verify heterozygous variation overlap
						for(size_t j=1; j < merge_item_vec.size(); j++){
							if(bench->sv_type == VAR_INS){
								startpos1 = merge_item_vec.at(j-1)->startPos;
								startpos2 = merge_item_vec.at(j)->startPos;
								endpos1 = merge_item_vec.at(j-1)->startPos + merge_item_vec.at(j-1)->ref_seq.length()-1;
								endpos2 = merge_item_vec.at(j)->endPos + merge_item_vec.at(j)->ref_seq.length()-1;
							}else if(bench->sv_type == VAR_DEL){
								startpos1 = merge_item_vec.at(j-1)->startPos;
								startpos2 = merge_item_vec.at(j)->startPos;
								endpos1 = merge_item_vec.at(j-1)->endPos;
								endpos2 = merge_item_vec.at(j)->endPos;
							}
							if(endpos1 >= startpos2 && endpos2 >= startpos1) return ;
						}
						if(bench->sv_type == VAR_INS){
							min_startpos = merge_item_vec.at(0)->startPos <= bench->startPos ? merge_item_vec.at(0)->startPos : bench->startPos;
							max_endpos = (merge_item_vec.at(0)->startPos + merge_item_vec.at(0)-> ref_seq.length()-1) >= (bench->startPos + bench->ref_seq.length()-1) ? (merge_item_vec.at(0)->startPos + merge_item_vec.at(0)-> ref_seq.length()-1) :(bench->startPos + bench->ref_seq.length()-1);
						}else if(bench->sv_type == VAR_DEL){
							min_startpos = merge_item_vec.at(0)->startPos <= bench->startPos ? merge_item_vec.at(0)->startPos : bench->startPos;
							max_endpos = merge_item_vec.at(0)->endPos >= bench->endPos ? merge_item_vec.at(0)->endPos : bench->endPos;
						}
						for(size_t j=0; j < merge_item_vec.size(); j++){
							Merge_users = merge_item_vec.at(j);
							if(Merge_users->sv_type == VAR_INS and bench->sv_type == VAR_INS){
								min_startpos = Merge_users->startPos <= min_startpos ? Merge_users->startPos : min_startpos;
								max_endpos = (Merge_users->startPos + Merge_users->ref_seq.length()-1) >= (size_t)max_endpos ? (Merge_users->startPos + Merge_users->ref_seq.length()-1) : max_endpos;
							}else if(user->sv_type == VAR_DEL and bench->sv_type == VAR_DEL){
								min_startpos = Merge_users->startPos <= min_startpos ? Merge_users->startPos : min_startpos;
								max_endpos = Merge_users->endPos >= max_endpos ? Merge_users->endPos : max_endpos;
							}
						}
						//Extract seq1 sequence
						for(size_t j=0; j < merge_item_vec.size(); j++){
							Merge_item = merge_item_vec.at(j);
							if(Merge_item->sv_type == VAR_INS){
								if(j==0){
									if(min_startpos == Merge_item->startPos){
										seq1 += Merge_item->alt_seq;
										/*reg_str = Merge_item->chrname + ":" + to_string(Merge_item->startPos+1) + "-" + to_string(merge_item_vec[j+1]->startPos-1);
										pthread_mutex_lock(&mutex_mem);
										seq = fai_fetch(fai, reg_str.c.str(), &refseq_len_tmp);
										pthread_mutex_unlock(&mutex_mem);
										seq_extact = seq;
										free(seq);
										seq1 += seq_extact;*/
									}else{
										reg_str = Merge_item->chrname + ":" + to_string(min_startpos) + "-" + to_string(Merge_item->startPos);
										pthread_mutex_lock(&mutex_mem);
										seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
										pthread_mutex_unlock(&mutex_mem);
										seq_extact = seq;
										free(seq);
										seq1 += seq_extact;
										seq1 += Merge_item->alt_seq;
									}
								}else if(j==merge_item_vec.size()-1){
									if((size_t)max_endpos == Merge_item->startPos+Merge_item->alt_seq.length()-1){
										reg_str = Merge_item->chrname + ":" + to_string(merge_item_vec[j-1]->startPos) + "-" + to_string(max_endpos);
										pthread_mutex_lock(&mutex_mem);
										seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
										pthread_mutex_unlock(&mutex_mem);
										seq_extact = seq;
										free(seq);
										seq1 += seq_extact;
										seq1 += Merge_item->alt_seq;
									}else{
										reg_str = Merge_item->chrname + ":" + to_string(merge_item_vec[j-1]->startPos) + "-" + to_string(Merge_item->startPos);
										pthread_mutex_lock(&mutex_mem);
										seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
										pthread_mutex_unlock(&mutex_mem);
										seq_extact = seq;
										free(seq);
										seq1 += seq_extact;
										seq1 += Merge_item->alt_seq;
										reg_str = Merge_item->chrname + ":" + to_string(Merge_item->startPos) + "-" + to_string(max_endpos);
										pthread_mutex_lock(&mutex_mem);
										seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
										pthread_mutex_unlock(&mutex_mem);
										seq_extact = seq;
										free(seq);
										seq1 += seq_extact;
									}
								}else{
									if(j+1 < merge_item_vec.size()){
										reg_str = Merge_item->chrname + ":" + to_string(merge_item_vec[j-1]->startPos) + "-" + to_string(Merge_item->startPos);
										pthread_mutex_lock(&mutex_mem);
										seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
										pthread_mutex_unlock(&mutex_mem);
										seq_extact = seq;
										free(seq);
										seq1 += seq_extact;
										seq1 += Merge_item->alt_seq;
									}
								}
							}else if(Merge_item->sv_type == VAR_DEL){
								if(j==0){
									if(min_startpos == Merge_item->startPos){
										seq1 += Merge_item->ref_seq;
										/*reg_str = Merge_item->chrname + ":" + to_string(Merge_item->startPos+1) + "-" + to_string(merge_item_vec[j+1]->startPos-1);
										pthread_mutex_lock(&mutex_mem);
										seq = fai_fetch(fai, reg_str.c.str(), &refseq_len_tmp);
										pthread_mutex_unlock(&mutex_mem);
										seq_extact = seq;
										free(seq);
										seq1 += seq_extact;*/
									}else{
										reg_str = Merge_item->chrname + ":" + to_string(min_startpos) + "-" + to_string(Merge_item->startPos);
										pthread_mutex_lock(&mutex_mem);
										seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
										pthread_mutex_unlock(&mutex_mem);
										seq_extact = seq;
										free(seq);
										seq1 += seq_extact;
										seq1 += Merge_item->ref_seq;
									}
								}else if(j==merge_item_vec.size()-1){
									if(max_endpos == Merge_item->endPos){
										reg_str = Merge_item->chrname + ":" + to_string(merge_item_vec[j-1]->endPos) + "-" + to_string(Merge_item->startPos);
										pthread_mutex_lock(&mutex_mem);
										seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
										pthread_mutex_unlock(&mutex_mem);
										seq_extact = seq;
										free(seq);
										seq1 += seq_extact;
										seq1 += Merge_item->ref_seq;
									}else{
										reg_str = Merge_item->chrname + ":" + to_string(merge_item_vec[j-1]->endPos) + "-" + to_string(Merge_item->startPos);
										pthread_mutex_lock(&mutex_mem);
										seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
										pthread_mutex_unlock(&mutex_mem);
										seq_extact = seq;
										free(seq);
										seq1 += seq_extact;
										seq1 += Merge_item->ref_seq;
										reg_str = Merge_item->chrname + ":" + to_string(Merge_item->endPos) + "-" + to_string(max_endpos);
										pthread_mutex_lock(&mutex_mem);
										seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
										pthread_mutex_unlock(&mutex_mem);
										seq_extact = seq;
										free(seq);
										seq1 += seq_extact;
									}
								}else{
									if(j+1 < merge_item_vec.size()){
										reg_str = Merge_item->chrname + ":" + to_string(merge_item_vec[j-1]->endPos) + "-" + to_string(Merge_item->startPos);
										pthread_mutex_lock(&mutex_mem);
										seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
										pthread_mutex_unlock(&mutex_mem);
										seq_extact = seq;
										free(seq);
										seq1 += seq_extact;
										seq1 += Merge_item->ref_seq;
									}
								}
							}
						}
						//Extract seq2 sequence
						if(bench->sv_type == VAR_INS){
							if(min_startpos == bench->startPos){
								seq2 += bench->alt_seq;
								reg_str = bench->chrname + ":" + to_string(min_startpos) + "-" + to_string(max_endpos);
								pthread_mutex_lock(&mutex_mem);
								seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
								pthread_mutex_unlock(&mutex_mem);
								seq_extact = seq;
								free(seq);
								seq2 += seq_extact;
							}else if(min_startpos < bench->startPos and max_endpos > bench->startPos){
								reg_str = bench->chrname + ":" + to_string(min_startpos) + "-" + to_string(bench->startPos);
								pthread_mutex_lock(&mutex_mem);
								seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
								pthread_mutex_unlock(&mutex_mem);
								seq_extact = seq;
								free(seq);
								seq2 += seq_extact;
								seq2 += bench->alt_seq;
								reg_str = bench->chrname + ":" + to_string(bench->startPos) + "-" + to_string(max_endpos);
								pthread_mutex_lock(&mutex_mem);
								seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
								pthread_mutex_unlock(&mutex_mem);
								seq_extact = seq;
								free(seq);
								seq2 += max_endpos;
							}else{
								reg_str = bench->chrname + ":" + to_string(min_startpos) + "-" + to_string(max_endpos);
								pthread_mutex_lock(&mutex_mem);
								seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
								pthread_mutex_unlock(&mutex_mem);
								seq_extact = seq;
								free(seq);
								seq2 += seq_extact;
								seq2 += bench->alt_seq;
							}
						}else if(bench->sv_type == VAR_DEL){
							if(min_startpos == bench->startPos){
								seq2 += bench->ref_seq;
								reg_str = bench->chrname + ":" + to_string(bench->endPos) + "-" + to_string(max_endpos);
								pthread_mutex_lock(&mutex_mem);
								seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
								pthread_mutex_unlock(&mutex_mem);
								seq_extact = seq;
								free(seq);
								seq2 += seq_extact;
							}else if(min_startpos < bench->startPos and max_endpos > bench->startPos){
								reg_str = bench->chrname + ":" + to_string(min_startpos) + "-" + to_string(bench->startPos);
								pthread_mutex_lock(&mutex_mem);
								seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
								pthread_mutex_unlock(&mutex_mem);
								seq_extact = seq;
								free(seq);
								seq2 += seq_extact;
								seq2 += bench->ref_seq;
								reg_str = bench->chrname + ":" + to_string(bench->endPos) + "-" + to_string(max_endpos);
								pthread_mutex_lock(&mutex_mem);
								seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
								pthread_mutex_unlock(&mutex_mem);
								seq_extact = seq;
								free(seq);
								seq2 += max_endpos;
							}else{
								reg_str = bench->chrname + ":" + to_string(min_startpos) + "-" + to_string(bench->startPos);
								pthread_mutex_lock(&mutex_mem);
								seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
								pthread_mutex_unlock(&mutex_mem);
								seq_extact = seq;
								free(seq);
								seq2 += seq_extact;
								seq2 += bench->ref_seq;
							}
						}
						//Sequence identity calculation
						upperSeq(seq1);
						upperSeq(seq2);
						needleman_wunsch(seq1, seq2, MATCH_SCORE, MISMATCH_SCORE, GAP_PENALTY, aln_seq1, aln_seq2);
						consistency = calculate_consistency(aln_seq1, aln_seq2);
						if(consistency >= SEQ_CONSISTENCY){
							for(size_t j=0; j < merge_item_vec.size(); j++){
								//LP flag
								Lp_Merge_item = merge_item_vec.at(j);
								if(Lp_Merge_item->LpFlag){
									flag = true;
									break;
								}
							}
							if(!flag) {LpNum += 1;
		//					cout << "merge Num: " << merge_item_vec.size() <<"	pos:" << merge_item_vec.at(0)->startPos <<"	LpNum = " << LpNum << endl;
							}
							for(size_t j=0; j < merge_item_vec.size(); j++){
								Definite_Merge_item = merge_item_vec.at(j);
								Definite_Merge_item->mergeMark = true;
								Definite_Merge_item->LpFlag = true;
							}
							bench->overlapped = true;
							bench->mergeMark = true;
		//					bench->seqcons = consistency;
						}
					}
				}
			}
		}else{
			//Merge verification: Whether it is consistent after the merger
			auto merge_result = findClosestSubarray(merge_users_vec, bench->sv_len);
			merge_item_vec = merge_result.second;
			if(merge_item_vec.size() > 1){
				//Verify heterozygous variation overlap
				for(size_t j=1; j < merge_item_vec.size(); j++){
					if(bench->sv_type == VAR_INS){
						startpos1 = merge_item_vec.at(j-1)->startPos;
						startpos2 = merge_item_vec.at(j)->startPos;
						endpos1 = merge_item_vec.at(j-1)->startPos + merge_item_vec.at(j-1)->ref_seq.length()-1;
						endpos2 = merge_item_vec.at(j)->endPos + merge_item_vec.at(j)->ref_seq.length()-1;
					}else if(bench->sv_type == VAR_DEL){
						startpos1 = merge_item_vec.at(j-1)->startPos;
						startpos2 = merge_item_vec.at(j)->startPos;
						endpos1 = merge_item_vec.at(j-1)->endPos;
						endpos2 = merge_item_vec.at(j)->endPos;
					}
					if(endpos1 >= startpos2 && endpos2 >= startpos1) return ;
				}
				if(bench->sv_type == VAR_INS){
					min_startpos = merge_item_vec.at(0)->startPos <= bench->startPos ? merge_item_vec.at(0)->startPos : bench->startPos;
					max_endpos = (merge_item_vec.at(0)->startPos + merge_item_vec.at(0)-> ref_seq.length()-1) >= (bench->startPos + bench->ref_seq.length()-1) ? (merge_item_vec.at(0)->startPos + merge_item_vec.at(0)-> ref_seq.length()-1) :(bench->startPos + bench->ref_seq.length()-1);
				}else if(bench->sv_type == VAR_DEL){
					min_startpos = merge_item_vec.at(0)->startPos <= bench->startPos ? merge_item_vec.at(0)->startPos : bench->startPos;
					max_endpos = merge_item_vec.at(0)->endPos >= bench->endPos ? merge_item_vec.at(0)->endPos : bench->endPos;
				}
				for(size_t j=0; j < merge_item_vec.size(); j++){
					Merge_users = merge_item_vec.at(j);
					if(Merge_users->sv_type == VAR_INS and bench->sv_type == VAR_INS){
						min_startpos = Merge_users->startPos <= min_startpos ? Merge_users->startPos : min_startpos;
						max_endpos = (Merge_users->startPos + Merge_users->ref_seq.length()-1) >= (size_t)max_endpos ? (Merge_users->startPos + Merge_users->ref_seq.length()-1) : max_endpos;
					}else if(user->sv_type == VAR_DEL and bench->sv_type == VAR_DEL){
						min_startpos = Merge_users->startPos <= min_startpos ? Merge_users->startPos : min_startpos;
						max_endpos = Merge_users->endPos >= max_endpos ? Merge_users->endPos : max_endpos;
					}
				}
				//Extract seq1 sequence
				for(size_t j=0; j < merge_item_vec.size(); j++){
					Merge_item = merge_item_vec.at(j);
					if(Merge_item->sv_type == VAR_INS){
						if(j==0){
							if(min_startpos == Merge_item->startPos){
								seq1 += Merge_item->alt_seq;
								/*reg_str = Merge_item->chrname + ":" + to_string(Merge_item->startPos+1) + "-" + to_string(merge_item_vec[j+1]->startPos-1);
								pthread_mutex_lock(&mutex_mem);
								seq = fai_fetch(fai, reg_str.c.str(), &refseq_len_tmp);
								pthread_mutex_unlock(&mutex_mem);
								seq_extact = seq;
								free(seq);
								seq1 += seq_extact;*/
							}else{
								reg_str = Merge_item->chrname + ":" + to_string(min_startpos) + "-" + to_string(Merge_item->startPos);
								pthread_mutex_lock(&mutex_mem);
								seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
								pthread_mutex_unlock(&mutex_mem);
								seq_extact = seq;
								free(seq);
								seq1 += seq_extact;
								seq1 += Merge_item->alt_seq;
							}
						}else if(j==merge_item_vec.size()-1){
							if((size_t)max_endpos == Merge_item->startPos+Merge_item->alt_seq.length()-1){
								reg_str = Merge_item->chrname + ":" + to_string(merge_item_vec[j-1]->startPos) + "-" + to_string(max_endpos);
								pthread_mutex_lock(&mutex_mem);
								seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
								pthread_mutex_unlock(&mutex_mem);
								seq_extact = seq;
								free(seq);
								seq1 += seq_extact;
								seq1 += Merge_item->alt_seq;
							}else{
								reg_str = Merge_item->chrname + ":" + to_string(merge_item_vec[j-1]->startPos) + "-" + to_string(Merge_item->startPos);
								pthread_mutex_lock(&mutex_mem);
								seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
								pthread_mutex_unlock(&mutex_mem);
								seq_extact = seq;
								free(seq);
								seq1 += seq_extact;
								seq1 += Merge_item->alt_seq;
								reg_str = Merge_item->chrname + ":" + to_string(Merge_item->startPos) + "-" + to_string(max_endpos);
								pthread_mutex_lock(&mutex_mem);
								seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
								pthread_mutex_unlock(&mutex_mem);
								seq_extact = seq;
								free(seq);
								seq1 += seq_extact;
							}
						}else{
							if(j+1 < merge_item_vec.size()){
								reg_str = Merge_item->chrname + ":" + to_string(merge_item_vec[j-1]->startPos) + "-" + to_string(Merge_item->startPos);
								pthread_mutex_lock(&mutex_mem);
								seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
								pthread_mutex_unlock(&mutex_mem);
								seq_extact = seq;
								free(seq);
								seq1 += seq_extact;
								seq1 += Merge_item->alt_seq;
							}
						}
					}else if(Merge_item->sv_type == VAR_DEL){
						if(j==0){
							if(min_startpos == Merge_item->startPos){
								seq1 += Merge_item->ref_seq;
								/*reg_str = Merge_item->chrname + ":" + to_string(Merge_item->startPos+1) + "-" + to_string(merge_item_vec[j+1]->startPos-1);
								pthread_mutex_lock(&mutex_mem);
								seq = fai_fetch(fai, reg_str.c.str(), &refseq_len_tmp);
								pthread_mutex_unlock(&mutex_mem);
								seq_extact = seq;
								free(seq);
								seq1 += seq_extact;*/
							}else{
								reg_str = Merge_item->chrname + ":" + to_string(min_startpos) + "-" + to_string(Merge_item->startPos);
								pthread_mutex_lock(&mutex_mem);
								seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
								pthread_mutex_unlock(&mutex_mem);
								seq_extact = seq;
								free(seq);
								seq1 += seq_extact;
								seq1 += Merge_item->ref_seq;
							}
						}else if(j==merge_item_vec.size()-1){
							if(max_endpos == Merge_item->endPos){
								reg_str = Merge_item->chrname + ":" + to_string(merge_item_vec[j-1]->endPos) + "-" + to_string(Merge_item->startPos);
								pthread_mutex_lock(&mutex_mem);
								seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
								pthread_mutex_unlock(&mutex_mem);
								seq_extact = seq;
								free(seq);
								seq1 += seq_extact;
								seq1 += Merge_item->ref_seq;
							}else{
								reg_str = Merge_item->chrname + ":" + to_string(merge_item_vec[j-1]->endPos) + "-" + to_string(Merge_item->startPos);
								pthread_mutex_lock(&mutex_mem);
								seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
								pthread_mutex_unlock(&mutex_mem);
								seq_extact = seq;
								free(seq);
								seq1 += seq_extact;
								seq1 += Merge_item->ref_seq;
								reg_str = Merge_item->chrname + ":" + to_string(Merge_item->endPos) + "-" + to_string(max_endpos);
								pthread_mutex_lock(&mutex_mem);
								seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
								pthread_mutex_unlock(&mutex_mem);
								seq_extact = seq;
								free(seq);
								seq1 += seq_extact;
							}
						}else{
							if(j+1 < merge_item_vec.size()){
								reg_str = Merge_item->chrname + ":" + to_string(merge_item_vec[j-1]->endPos) + "-" + to_string(Merge_item->startPos);
								pthread_mutex_lock(&mutex_mem);
								seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
								pthread_mutex_unlock(&mutex_mem);
								seq_extact = seq;
								free(seq);
								seq1 += seq_extact;
								seq1 += Merge_item->ref_seq;
							}
						}
					}
				}
				//Extract seq2 sequence
				if(bench->sv_type == VAR_INS){
					if(min_startpos == bench->startPos){
						seq2 += bench->alt_seq;
						reg_str = bench->chrname + ":" + to_string(min_startpos) + "-" + to_string(max_endpos);
						pthread_mutex_lock(&mutex_mem);
						seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
						pthread_mutex_unlock(&mutex_mem);
						seq_extact = seq;
						free(seq);
						seq2 += seq_extact;
					}else if(min_startpos < bench->startPos and max_endpos > bench->startPos){
						reg_str = bench->chrname + ":" + to_string(min_startpos) + "-" + to_string(bench->startPos);
						pthread_mutex_lock(&mutex_mem);
						seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
						pthread_mutex_unlock(&mutex_mem);
						seq_extact = seq;
						free(seq);
						seq2 += seq_extact;
						seq2 += bench->alt_seq;
						reg_str = bench->chrname + ":" + to_string(bench->startPos) + "-" + to_string(max_endpos);
						pthread_mutex_lock(&mutex_mem);
						seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
						pthread_mutex_unlock(&mutex_mem);
						seq_extact = seq;
						free(seq);
						seq2 += max_endpos;
					}else{
						reg_str = bench->chrname + ":" + to_string(min_startpos) + "-" + to_string(max_endpos);
						pthread_mutex_lock(&mutex_mem);
						seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
						pthread_mutex_unlock(&mutex_mem);
						seq_extact = seq;
						free(seq);
						seq2 += seq_extact;
						seq2 += bench->alt_seq;
					}
				}else if(bench->sv_type == VAR_DEL){
					if(min_startpos == bench->startPos){
						seq2 += bench->ref_seq;
						reg_str = bench->chrname + ":" + to_string(bench->endPos) + "-" + to_string(max_endpos);
						pthread_mutex_lock(&mutex_mem);
						seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
						pthread_mutex_unlock(&mutex_mem);
						seq_extact = seq;
						free(seq);
						seq2 += seq_extact;
					}else if(min_startpos < bench->startPos and max_endpos > bench->startPos){
						reg_str = bench->chrname + ":" + to_string(min_startpos) + "-" + to_string(bench->startPos);
						pthread_mutex_lock(&mutex_mem);
						seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
						pthread_mutex_unlock(&mutex_mem);
						seq_extact = seq;
						free(seq);
						seq2 += seq_extact;
						seq2 += bench->ref_seq;
						reg_str = bench->chrname + ":" + to_string(bench->endPos) + "-" + to_string(max_endpos);
						pthread_mutex_lock(&mutex_mem);
						seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
						pthread_mutex_unlock(&mutex_mem);
						seq_extact = seq;
						free(seq);
						seq2 += max_endpos;
					}else{
						reg_str = bench->chrname + ":" + to_string(min_startpos) + "-" + to_string(bench->startPos);
						pthread_mutex_lock(&mutex_mem);
						seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
						pthread_mutex_unlock(&mutex_mem);
						seq_extact = seq;
						free(seq);
						seq2 += seq_extact;
						seq2 += bench->ref_seq;
					}
				}
				//Sequence identity calculation
				upperSeq(seq1);
				upperSeq(seq2);
				needleman_wunsch(seq1, seq2, MATCH_SCORE, MISMATCH_SCORE, GAP_PENALTY, aln_seq1, aln_seq2);
				consistency = calculate_consistency(aln_seq1, aln_seq2);
				if(consistency >= SEQ_CONSISTENCY){
					for(size_t j=0; j < merge_item_vec.size(); j++){
						//LP flag
						Lp_Merge_item = merge_item_vec.at(j);
						if(Lp_Merge_item->LpFlag){
							flag = true;
							break;
						}
					}
					if(!flag) {LpNum += 1;
//					cout << "merge Num: " << merge_item_vec.size() <<"	pos:" << merge_item_vec.at(0)->startPos <<"	LpNum = " << LpNum << endl;
					}
					for(size_t j=0; j < merge_item_vec.size(); j++){
						Definite_Merge_item = merge_item_vec.at(j);
						Definite_Merge_item->mergeMark = true;
						Definite_Merge_item->LpFlag = true;
					}
					bench->overlapped = true;
					bench->mergeMark = true;
//					bench->seqcons = consistency;
				}
			}
		}
	}
	return ;
}
// Finds the contiguous SV_item subarray that is closest to the target value and satisfies the condition
pair<int, vector<SV_item*>> findClosestSubarray(vector<SV_item*> items, int32_t target) {
    int closestSVlenSum = 0; // Used to record the nearest SVlen Sum
    int minDiff = numeric_limits<int>::max();
    vector<SV_item*> bestSubarray;
    double merge_len_ratio;

    // Go through all possible starting positions
    for (size_t i = 0; i < items.size(); ++i) {
        int currentSVlenSum = 0;
        vector<SV_item*> currentSubarray;
        // Go through all possible end positions
        for (size_t j = i; j < items.size(); ++j) {
        	currentSVlenSum += items[j]->sv_len;
            currentSubarray.push_back(items[j]);
            merge_len_ratio = currentSVlenSum >= target ? (double)target / currentSVlenSum : (double)currentSVlenSum / target;
            // Determines whether the SV size condition is met
            if (merge_len_ratio >= 0.9) {
                int currentDiff = abs(currentSVlenSum - target);
                //Update the nearest SVlen and subarray
                if (currentDiff < minDiff) {
                    minDiff = currentDiff;
                    closestSVlenSum = currentSVlenSum;
                    bestSubarray = currentSubarray;
                }else if(currentDiff == minDiff){
                	break;
                }
            }
        }
    }
    return { closestSVlenSum, bestSubarray }; // Returns the nearest sum and subarray
}

int32_t minDistance(const string &seq1, const string &seq2) {
	int64_t m = seq1.size() + 1, n = seq2.size() + 1;
	int64_t i, j, arrSize;
	int32_t *dp;

	// Initialize the score matrix
	arrSize = m * n;
	dp = (int32_t*) calloc (arrSize, sizeof(int32_t));
	if(dp==NULL){
		cerr << "line=" << __LINE__ << ", m=" << m << ", n=" << n << ", cannot allocate memory, error!" << endl;
		exit(1);
	}

	// initialize the first row and first colum
	for (i = 0; i < m; ++i) dp[i*n] = i;
	for (j = 0; j < n; ++j) dp[j] = j;

	// compute distance by dynamic program
	for (i = 1; i < m; ++i) {
		for (j = 1; j < n; ++j) {
			if (seq1[i - 1] == seq2[j - 1]) {
				dp[i*n+j] = dp[(i-1)*n+j-1];
			} else {
				dp[i*n+j] = 1 + min(dp[(i-1)*n+j], dp[i*n+j-1], dp[(i-1)*n+j-1]);
			}
		}
	}
	int32_t result = dp[(m-1)*n+n-1];
	free(dp);

//	// initialize the matrix with size: m*n
//	vector<vector<int>> dp(m, vector<int>(n, 0));
//
//	// initialize the first row and first colum
//	for (int i = 0; i <= m; ++i)
//		dp[i][0] = i;
//	for (int j = 0; j <= n; ++j)
//		dp[0][j] = j;
//
//	// compute distance by dynamic program
//	for (int i = 1; i <= m; ++i) {
//		for (int j = 1; j <= n; ++j) {
//			if (seq1[i - 1] == seq2[j - 1]) {
//				dp[i][j] = dp[i - 1][j - 1];
//			} else {
//				dp[i][j] = 1 + min({dp[i - 1][j], dp[i][j - 1], dp[i - 1][j - 1]});
//			}
//		}
//	}
	return result;
}

double computeVarseqConsistency(SV_item *item1, SV_item *item2, faidx_t *fai, double SeqIdentity){
	string seq1, seq2, aln_seq1, aln_seq2,  AlignSeq,  AlignSeq1, aln_seq3, aln_seq4;
	double consistency,consistency_tmp1, consistency_tmp2;
	vector<Minimizer> minimizers, minimizers1;
	seq1 = seq2 = aln_seq1 = aln_seq2 = aln_seq3 = aln_seq4 = "";

	if(item1->sv_type == VAR_INS and item2->sv_type == VAR_INS){
		if(item1->alt_seq.size() < MAX_SEQ_LEN and item2->alt_seq.size() < MAX_SEQ_LEN) {
			seq1 = item1->alt_seq;	seq2 = item2->alt_seq;
			upperSeq(seq1);
			upperSeq(seq2);
//			cout << "seq1: "<< seq1 << endl;
//			cout << "seq2: "<< seq2 << endl;
			needleman_wunsch(seq1, seq2, MATCH_SCORE, MISMATCH_SCORE, GAP_PENALTY, aln_seq1, aln_seq2);
			consistency_tmp1 = calculate_consistency(aln_seq1, aln_seq2);
			extractRefSeq(item1, item2, seq1, seq2, fai);
			needleman_wunsch(seq1, seq2, MATCH_SCORE, MISMATCH_SCORE, GAP_PENALTY, aln_seq3, aln_seq4);
			consistency_tmp2 = calculate_consistency(aln_seq3, aln_seq4);
//			if(consistency_tmp1 > consistency_tmp2)
//				cout << "chr: " << item1->chrname << " startPos: " << item1->startPos << " SVtype: " << item1->sv_type << "SVLen: " << item1->sv_len << " item2_startpos: " << item2->startPos << " consistency_tmp1:" <<consistency_tmp1<<" consistency_tmp2:" << consistency_tmp2 << endl;
			consistency = (consistency_tmp1 > consistency_tmp2) ? consistency_tmp1 : consistency_tmp2;
//			cout << "seq1: "<< seq1 << endl;
//			cout << "seq2: "<< seq2 << endl;
			if(consistency >= SeqIdentity)	return consistency;
			else{
				int distance = minDistance(seq1, seq2);
				consistency =  1 - ((double)distance / (seq1.length() + seq2.length()));
				if(consistency >= SeqIdentity) return consistency;
				else seq1 = seq2 = aln_seq1 = aln_seq2 = "";
			}
		}else{
			seq1 = item1->alt_seq;	seq2 = item2->alt_seq;
			upperSeq(seq1);
			upperSeq(seq2);
			consistency = MinimizerMethodOp(seq1, seq2, aln_seq1, aln_seq2, AlignSeq, AlignSeq1);
			if(consistency >= SeqIdentity)	return consistency;
			else	seq1 = seq2 = aln_seq1 = aln_seq2 = AlignSeq = AlignSeq1 = "";
		}
	}
	// extract the nucleotide sequences and add to both sides of variant sequences
	extractRefSeq(item1, item2, seq1, seq2, fai);
//	cout<<"exactSeq1:"<<seq1<<endl;
//	cout<<"exactSeq2:"<<seq2<<endl;
	//Determine which one to execute based on the sequence size
	if(seq1.size() >= MAX_SEQ_LEN or seq2.size() >= MAX_SEQ_LEN){
		upperSeq(seq1);
		upperSeq(seq2);
		consistency = MinimizerMethodOp(seq1, seq2, aln_seq1, aln_seq2, AlignSeq, AlignSeq1);
		/*//Find minimizer
		minimizers = findMinimizers(seq1, WINDOWSIZE, KMERSIZE);
		minimizers1 = findMinimizers(seq2, WINDOWSIZE, KMERSIZE);
		//Find a minimizer with a consistent sequence
		retainCommonHashValuesMinimizers(minimizers, minimizers1);
		//Look for minimizers in similar positions
		findSimilarityPosMinimizers(minimizers, minimizers1);
		//Filter minimizers with dissimilar positions.
		FilterDissimilaPosMinimizers(minimizers, minimizers1);
		//Segmentation sequence alignment
		LongSequenceSplitAlignment(minimizers, minimizers1,  KMERSIZE,  AlignSeq,  AlignSeq1,  seq1,  seq2);
		//Check whether minimizer exists after filtering
		if (minimizers.empty() or minimizers1.empty()) {
		        needleman_wunsch(seq1, seq2, MATCH_SCORE, MISMATCH_SCORE,  GAP_PENALTY, aln_seq1, aln_seq2);
		        consistency = calculate_consistency(aln_seq1, aln_seq2);
		}else
			consistency = calculate_consistency(AlignSeq, AlignSeq1);*/
	}else{
	// compute the optimal match result by NW alignment
		upperSeq(seq1);
		upperSeq(seq2);
		needleman_wunsch(seq1, seq2, MATCH_SCORE, MISMATCH_SCORE, GAP_PENALTY, aln_seq1, aln_seq2);
	// compute the consistency according to alignment result
		consistency = calculate_consistency(aln_seq1, aln_seq2);
	}
//	consistency = 0.9;
	/*if(item1->sv_type == VAR_INS and consistency < percentSeqIdentity){
			cout<<"111111111111: "<<item1->chrname<<" :"<<item1->startPos<<endl;
		}*/
	if(consistency < SeqIdentity and seq1.size() <= EDLIBLEN){
		int distance = minDistance(seq1, seq2);
		consistency =  1 - ((double)distance / (seq1.length() + seq2.length()));
	}

	// compute the edit distance
//	int32_t distance = minDistance(seq1, seq2);
//	int32_t seq1_len = seq1.size(); int seq2_len = seq2.size();
//	double discrepancy_score = (double)distance/(seq1_len+seq2_len);
//	consistency = 1-discrepancy_score;

	return consistency;
}

void SeqConsNumStat(double consistency){
	SeqconsSum += consistency;
	SeqconsNum += 1;
//	cout << "num: " << SeqconsNum << "value: " << consistency <<" SUM:"<< SeqconsSum << endl;
}

// extract the reference sequences
void extractRefSeq(SV_item* item1, SV_item* item2, string &seq_new1, string &seq_new2, faidx_t *fai) {
	string seq_left, seq_right, reg_str;
	int startpos1, startpos2, endpos1, endpos2;
	int32_t refseq_len_tmp;
	char *seq;
	/*int s,e;
	string q;*/

	// If any of the alternate sequences is "-", no manipulation needed
//	if (item1->alt_seq.compare("-") == 0 or item2->alt_seq.compare("-") == 0) {
//		return;
//	}
	/*if(item1->sv_type == VAR_INS || item2->sv_type == VAR_INS)
		cout<<"INS:11111111111111"<<endl;*/
	if(item1->sv_type == VAR_DEL || item2->sv_type == VAR_DEL){
		seq_new1 = item1->ref_seq;
		seq_new2 = item2->ref_seq;
	}else{
		seq_new1 = item1->alt_seq;
		seq_new2 = item2->alt_seq;
	}

	startpos1 = item1->startPos;
	startpos2 = item2->startPos;
	if((item1->sv_type == VAR_INS and item2->sv_type == VAR_INS) or (item1->sv_type == VAR_DEL || item2->sv_type == VAR_DEL)){
		endpos1 = item1->startPos + item1->ref_seq.length() - 1;
		endpos2 = item2->startPos + item2->ref_seq.length() - 1;
	}else{
		endpos1 = item1->endPos;
		endpos2 = item2->endPos;
	}
	/*if(startpos1<=startpos2)
		s = startpos1;
	else
		s = startpos2;
	if(endpos1<=endpos2)
		e=endpos2;
	else
		e=endpos1;
	reg_str = item1->chrname + ":" + to_string(s-20) + "-" + to_string(e+20);
	pthread_mutex_lock(&mutex_mem);
	seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
	pthread_mutex_unlock(&mutex_mem);
	q = seq;
	cout<<"The extended sequence is："<<q<<endl;
	free(seq);*/
//	if (item1->sv_type == VAR_INS and item1->startPos == item1->endPos) {
//		endpos1 = item1->endPos + item1->sv_len;
//	}else{
//		endpos1 = item1->endPos;
//	}
//	if (item2->sv_type == VAR_INS and item2->startPos == item2->endPos) {
//		endpos2 = item2->endPos + item2->sv_len;
//	} else {
//		endpos2 = item2->endPos;
//	}

	if (startpos1 < startpos2) {
//		if (startpos1 == startpos2)
//			reg_str = item1->chrname + ":" + to_string(startpos1) + "-" + to_string(startpos2);
//		else
		reg_str = item1->chrname + ":" + to_string(startpos1) + "-" + to_string(startpos2 - 1);
		pthread_mutex_lock(&mutex_mem);
		seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
		pthread_mutex_unlock(&mutex_mem);
		seq_left = seq;
		free(seq);

		if (item1->sv_type == VAR_INV and item2->sv_type == VAR_INV)
			reverseComplement(seq_left);
		seq_new2 = seq_left + seq_new2;
	} else if(startpos1 > startpos2) {
		reg_str = item1->chrname + ":" + to_string(startpos2) + "-" + to_string(startpos1 - 1);
		pthread_mutex_lock(&mutex_mem);
		seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
		pthread_mutex_unlock(&mutex_mem);
		seq_left = seq;
		free(seq);

		if (item1->sv_type == VAR_INV and item2->sv_type == VAR_INV)
			reverseComplement(seq_left);

		seq_new1 = seq_left + seq_new1;
	}

	if (endpos1 < endpos2) {
//		if (endpos1 == endpos2)
//			reg_str = item1->chrname + ":" + to_string(endpos1) + "-" + to_string(endpos2);
//		else
		reg_str = item1->chrname + ":" + to_string(endpos1 + 1) + "-" + to_string(endpos2);
		pthread_mutex_lock(&mutex_mem);
		seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
		pthread_mutex_unlock(&mutex_mem);
		seq_right = seq;
		free(seq);

		if (item1->sv_type == VAR_INV and item2->sv_type == VAR_INV)
			reverseComplement(seq_right);

		seq_new1 = seq_new1 + seq_right;
	} else if(endpos1 > endpos2){
		reg_str = item1->chrname + ":" + to_string(endpos2 + 1) + "-" + to_string(endpos1);
		pthread_mutex_lock(&mutex_mem);
		seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
		pthread_mutex_unlock(&mutex_mem);
		seq_right = seq;
		free(seq);

		if (item1->sv_type == VAR_INV and item2->sv_type == VAR_INV)
			reverseComplement(seq_right);

		seq_new2 = seq_new2 + seq_right;
	}
}

// compute the maximum value of three integers
int32_t max(int32_t a, int32_t b, int32_t c) {
	return max(a, max(b, c));
}

// compute the minimum value of three integers
int32_t min(int32_t a, int32_t b, int32_t c) {
	return min(a, min(b, c));
}

double MinimizerMethodOp(string& seq1, string& seq2, string& aln_seq1, string& aln_seq2, string& AlignSeq, string& AlignSeq1){
	vector<Minimizer> minimizers,  minimizers1;
	double result;
	//Find minimizer
	minimizers = findMinimizers(seq1, WINDOWSIZE, KMERSIZE);
	minimizers1 = findMinimizers(seq2, WINDOWSIZE, KMERSIZE);
	//Find a minimizer with a consistent sequence
	retainCommonHashValuesMinimizers(minimizers, minimizers1);
	//Look for minimizers in similar positions
	findSimilarityPosMinimizers(minimizers, minimizers1);
	//Filter minimizers with dissimilar positions.
	FilterDissimilaPosMinimizers(minimizers, minimizers1);
	//Segmentation sequence alignment
//	LongSequenceSplitAlignment(minimizers, minimizers1,  KMERSIZE,  AlignSeq,  AlignSeq1,  seq1,  seq2);
	//Check whether minimizer exists after filtering
	if (minimizers.empty() or minimizers1.empty()) {
		needleman_wunsch(seq1, seq2, MATCH_SCORE, MISMATCH_SCORE,  GAP_PENALTY, aln_seq1, aln_seq2);
		result = calculate_consistency(aln_seq1, aln_seq2);
	}else{
		//Segmentation sequence alignment
		LongSequenceSplitAlignment(minimizers, minimizers1,  KMERSIZE,  AlignSeq,  AlignSeq1,  seq1,  seq2);
		result = calculate_consistency(AlignSeq, AlignSeq1);
	}
	return result;
}

size_t customHashFunction(const string& kmer, size_t kmerSize) {
    size_t hashValue = 0;
    size_t base = 4;  //

    for (size_t i = 0; i < kmer.length(); ++i) {
        char baseChar = kmer[i];

        //
        size_t baseValue;
        switch (baseChar) {
        case 'A':
            baseValue = 0;
            break;
        case 'C':
            baseValue = 1;
            break;
        case 'G':
            baseValue = 2;
            break;
        case 'T':
            baseValue = 3;
            break;
        default:
            //The abnormal character is assigned to the k-mer maximum value
//            cerr << "Error: Invalid character in kmer." << endl;
            hashValue = INT_MAX;
            return hashValue;  //
        }

        //
        hashValue += baseValue * pow(base, kmerSize-i-1);
    }
    return hashValue;
}

vector<Minimizer> findMinimizers(const string& sequence, size_t windowSize, size_t kmerSize) {
    vector<Minimizer> minimizers;

    for (size_t i = 0; i <= sequence.length() - windowSize- kmerSize+1; ++i) {
        string window = sequence.substr(i, windowSize+ kmerSize-1);
        //size_t Maximum value of the data type
        size_t minHash = numeric_limits<size_t>::max();
        string minKmer = "";
        size_t pos;

        for (size_t j = 0; j <= window.size() - kmerSize; ++j) {
            string kmer = window.substr(j, kmerSize);
            size_t hashValue = customHashFunction(kmer, kmerSize);

            if (hashValue < minHash) {
                minHash = hashValue;
                minKmer = kmer;
                pos = i + j;
            }

        }
        Minimizer minimizer(pos, minHash, minKmer, false);
        if (!minimizers.empty()) {
            const Minimizer& lastMinimizer = minimizers.back();
            if(minHash!=lastMinimizer.hashValue or pos != lastMinimizer.position)
                minimizers.push_back(minimizer);
        }
        else {
            minimizers.push_back(minimizer);
        }
    }
    return minimizers;
}

void retainCommonHashValuesMinimizers(vector<Minimizer>& containerA, vector<Minimizer>& containerB) {
    //
    unordered_multiset<size_t> setB;
    for (const auto& element : containerB) {
        setB.insert(element.hashValue);
    }

    containerA.erase(remove_if(containerA.begin(), containerA.end(),
        [&setB](const Minimizer& element) {
            return setB.find(element.hashValue) == setB.end();
        }),
        containerA.end());

    unordered_multiset<size_t> setA;
    for (const auto& element : containerA) {
        setA.insert(element.hashValue);
    }
    //
    containerB.erase(remove_if(containerB.begin(), containerB.end(),
        [&setA](const Minimizer& element) {
            return setA.find(element.hashValue) == setA.end();
        }),
        containerB.end());
}

void findSimilarityPosMinimizers(vector<Minimizer>& minimizers, vector<Minimizer>& minimizers1) {

    int record_pos = -1;
    if (minimizers.empty() or minimizers1.empty())  return;
    if (minimizers.size() <= minimizers1.size()) {
        for (int i = 0; i < (int)minimizers.size(); i++) {
            int start = max(0,i - EXTEND_NUM);
            int end = min(i + EXTEND_NUM, (int)minimizers1.size()-1);

            if (record_pos != -1 and i != 0  ) {
                if (minimizers[record_pos].matchflag and abs((int)minimizers[i].position - (int)minimizers[record_pos].position) <= MAXLEN) continue;
            }
            //another
            int MinPos = INT_MAX;
            int j_subscript = -1;

            for (int j = start; j <= end; j++) {

                if (minimizers1[j].matchflag == true)   continue;
                if ((minimizers[i].hashValue == minimizers1[j].hashValue) and (!minimizers[i].matchflag and !minimizers1[j].matchflag)) {
                    if (abs((int)minimizers1[i].position - (int)minimizers[j].position) < MinPos)
                    {
                        MinPos = abs((int)minimizers1[i].position - (int)minimizers[j].position);
                        j_subscript = j;
                    }
                }
            }
            //another
            if (j_subscript != -1 and MinPos < MAXLEN ){
                minimizers[i].matchflag = minimizers1[j_subscript].matchflag = true;
                 record_pos = i;
             }
        }
    } else {
        for (int i = 0; i < (int)minimizers1.size(); i++) {
            int start = max(0, i - EXTEND_NUM);
            int end = min(i + EXTEND_NUM, (int)minimizers.size()-1);

            if (record_pos != -1 and i != 0 ) {
                if (minimizers1[record_pos].matchflag and abs((int)minimizers1[i].position - (int)minimizers1[record_pos].position) <= 15) continue;
            }
            //another
            int MinPos = INT_MAX;
            int j_subscript = -1;
            for (int j = start; j <= end; j++) {
                if (minimizers[j].matchflag == true)   continue;
                    //another
                if ((minimizers1[i].hashValue == minimizers[j].hashValue) and (!minimizers1[i].matchflag and !minimizers[j].matchflag)) {
                    if (abs((int)minimizers1[i].position - (int)minimizers[j].position) < MinPos)
                    {
                        MinPos = abs((int)minimizers1[i].position - (int)minimizers[j].position);
                        j_subscript = j;
                    }
                }
            }
            //another
            if (j_subscript != -1 and MinPos < MAXLEN) {
                minimizers1[i].matchflag = minimizers[j_subscript].matchflag = true;
                record_pos = i;
            }
        }
    }
}

void FilterDissimilaPosMinimizers(vector<Minimizer>& minimizers, vector<Minimizer>& minimizers1) {
    //
    minimizers.erase(remove_if(minimizers.begin(), minimizers.end(),
        [](const Minimizer& m) { return !m.matchflag; }),
        minimizers.end());

    minimizers1.erase(remove_if(minimizers1.begin(), minimizers1.end(),
        [](const Minimizer& m) { return !m.matchflag; }),
        minimizers1.end());
}

void needleman_wunsch(const string &seq1, const string &seq2, int32_t match_score, int32_t mismatch_score, int32_t gap_penalty, string &seq1_new, string &seq2_new){
	int64_t rowsNum = seq1.size() + 1, colsNum = seq2.size() + 1;
	int64_t arrSize = rowsNum * colsNum;
	int64_t mem_cost = (arrSize * sizeof(int) + (seq1.size() + seq2.size())*2) >> 10; // divided by 1024
	bool aln_flag = false;

	while(1){
		pthread_mutex_lock(&mutex_mem);
		if(mem_seqAln+mem_cost<=mem_total*mem_use_block_factor+extend_total*extend_use_block_factor){ // prepare for alignment computation
			mem_seqAln += mem_cost;
			aln_flag = true;
		}
		pthread_mutex_unlock(&mutex_mem);

		if(aln_flag==false){ // block the alignment computation
			//cout << "\t" << __func__ << ": wait " << mem_wait_seconds << " seconds, rowsNum=" << rowsNum << ", colsNum=" << colsNum << endl;
			sleep(mem_wait_seconds);
		}else break;
	}
	if(aln_flag){
		needleman_wunschOp(seq1, seq2, match_score, mismatch_score, gap_penalty, seq1_new, seq2_new);
		// update memory consumption
		pthread_mutex_lock(&mutex_mem);
		mem_seqAln -= mem_cost;
		if(mem_seqAln<0){
			cerr << "line=" << __LINE__ << ", mem_seqAln=" << mem_seqAln << ", error." << endl;
			exit(1);
		}
		pthread_mutex_unlock(&mutex_mem);
	}
}

// compute the optimal path by NW (Needleman Wunsch) alignment
void needleman_wunschOp(const string &seq1, const string &seq2, int32_t match_score, int32_t mismatch_score, int32_t gap_penalty, string &seq1_new, string &seq2_new){
	int64_t rowsNum = seq1.size() + 1, colsNum = seq2.size() + 1, arrSize;
	int32_t i, j, *score_matrix;
	int32_t score, diag_score, up_score;

	// Initialize the score matrix
	arrSize = rowsNum * colsNum;
	score_matrix = (int32_t*) calloc (arrSize, sizeof(int32_t));
	if(score_matrix==NULL){
		cerr << "line=" << __LINE__ << ", rowsNum=" << rowsNum << ", colsNum=" << colsNum << ", cannot allocate memory, error!" << endl;
		exit(1);
	}

	// Initialize the first row and first column
	for (i=1; i<rowsNum; i++) score_matrix[i*colsNum] = score_matrix[(i-1)*colsNum] + gap_penalty;
	for (j=1; j<colsNum; j++) score_matrix[j] = score_matrix[j-1] + gap_penalty;

	// Fill in the rest of the score matrix
	int32_t match, delete_gap, insert_gap;
	for (i=1; i<rowsNum; i++) {
		for (j=1; j<colsNum; j++) {
			match = score_matrix[(i-1)*colsNum+j-1] + (seq1[i-1]==seq2[j-1] ? match_score : mismatch_score);
			delete_gap = score_matrix[(i-1)*colsNum+j] + gap_penalty;
			insert_gap = score_matrix[i*colsNum+j-1] + gap_penalty;
			score_matrix[i*colsNum+j] = max(match, delete_gap, insert_gap);
		}
	}

	// Traceback to find the alignment
	i = rowsNum - 1, j = colsNum - 1;
	while (i>0 && j>0) {
		score = score_matrix[i*colsNum+j];
		diag_score = score_matrix[(i-1)*colsNum+j-1];
		//left_score = score_matrix[i*colsNum+j-1];
		up_score = score_matrix[(i-1)*colsNum+j];

		if (score == diag_score + (seq1[i-1] == seq2[j-1] ? match_score : mismatch_score)) {
			seq1_new = seq1[i-1] + seq1_new;
			seq2_new = seq2[j-1] + seq2_new;
			i--;
			j--;
		}else if (score == up_score + gap_penalty) {
			seq1_new = seq1[i-1] + seq1_new;
			seq2_new = "-" + seq2_new;
			i--;
		}else {
			seq1_new = "-" + seq1_new;
			seq2_new = seq2[j-1] + seq2_new;
			j--;
		}
	}

	while (i > 0) {
		seq1_new = seq1[i-1] + seq1_new;
		seq2_new = "-" + seq2_new;
		i--;
	}

	while (j > 0) {
		seq1_new = "-" + seq1_new;
		seq2_new = seq2[j-1] + seq2_new;
		j--;
	}

	free(score_matrix);


//	int64_t rows = seq1.length() + 1, cols = seq2.length() + 1, arrSize;
//
//	// Initialize the score matrix
//	vector<vector<int>> score_matrix(rows, vector<int>(cols, 0));
//
//	// Initialize the first row and first column
//	for (int i = 1; i < rows; i++) {
//		score_matrix[i][0] = score_matrix[i - 1][0] + gap_penalty;
//	}
//
//	for (int j = 1; j < cols; j++) {
//		score_matrix[0][j] = score_matrix[0][j - 1] + gap_penalty;
//	}
//
//	// Fill in the rest of the score matrix
//	int match, delete_gap, insert_gap;
//	for (int i = 1; i < rows; i++) {
//		for (int j = 1; j < cols; j++) {
//			match = score_matrix[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? match_score : mismatch_score);
//			delete_gap = score_matrix[i - 1][j] + gap_penalty;
//			insert_gap = score_matrix[i][j - 1] + gap_penalty;
//			score_matrix[i][j] = max(match, delete_gap, insert_gap);
//		}
//	}
//
//	// Traceback to find the alignment
//	int64_t i = rows - 1, j = cols - 1;
//	int score, diag_score, left_score, up_score;
//	while (i > 0 && j > 0) {
//		score = score_matrix[i][j];
//		diag_score = score_matrix[i - 1][j - 1];
//		left_score = score_matrix[i][j - 1];
//		up_score = score_matrix[i - 1][j];
//
//		if (score == diag_score + (seq1[i - 1] == seq2[j - 1] ? match_score : mismatch_score)) {
//			seq1_new = seq1[i - 1] + seq1_new;
//			seq2_new = seq2[j - 1] + seq2_new;
//			i--;
//			j--;
//		}
//		else if (score == up_score + gap_penalty) {
//			seq1_new = seq1[i - 1] + seq1_new;
//			seq2_new = "-" + seq2_new;
//			i--;
//		}
//		else {
//			seq1_new = "-" + seq1_new;
//			seq2_new = seq2[j - 1] + seq2_new;
//			j--;
//		}
//	}
//
//	while (i > 0) {
//		seq1_new = seq1[i - 1] + seq1_new;
//		seq2_new = "-" + seq2_new;
//		i--;
//	}
//
//	while (j > 0) {
//		seq1_new = "-" + seq1_new;
//		seq2_new = seq2[j - 1] + seq2_new;
//		j--;
//	}

//	cout << "Aligned Sequence 1: " << seq1_new << endl;
//	cout << "Aligned Sequence 2: " << seq2_new << endl;
//	cout << "Aligned Sequence 1(size):" << seq1_new.size() << endl;
//	cout << "Aligned Sequence 2(size):" << seq2_new.size() << endl;
}

void LongSequenceSplitAlignment(vector<Minimizer>& minimizers, vector<Minimizer>& minimizers1, int k, string& AlignSeq, string& AlignSeq1, string sequence, string sequence1) {

    int Exact_len, Exact_len1;
    string seq = "", seq1 = "", temporarily_seq = "", temporarily_seq1 ="";
    if (minimizers.size() != minimizers1.size() or (minimizers.empty() or minimizers1.empty())) return;
    int Length = minimizers.size(), sequece_len = sequence.size(), seqence_len1 = sequence1.size();
    for (int i = 0; i <= Length; i++) {
        if (i == 0) {
            seq = sequence.substr(0 , minimizers[i].position);
            seq1 = sequence.substr(0, minimizers1[i].position);
        }
        else if (i == Length) {
            Exact_len = sequece_len - (minimizers[i - 1].position + k);
            Exact_len1 = seqence_len1 - (minimizers1[i - 1].position + k);
            if(Exact_len < 0) Exact_len = 0;
            if(Exact_len1 < 0) Exact_len1 = 0;
            seq = sequence.substr(minimizers[i - 1].position + k, Exact_len);
            seq1 = sequence1.substr(minimizers1[i - 1].position + k, Exact_len1);
        }
        else {
            Exact_len = minimizers[i].position - (minimizers[i - 1].position + k);
            Exact_len1 = minimizers1[i].position - (minimizers1[i - 1].position + k);
            if(Exact_len < 0) Exact_len = 0;
            if(Exact_len1 < 0) Exact_len1 = 0;
            seq = sequence.substr(minimizers[i-1].position + k, Exact_len);
            seq1 = sequence1.substr(minimizers1[i - 1].position + k, Exact_len1);
        }
        //
        needleman_wunschOp(seq, seq1, MATCH_SCORE, MISMATCH_SCORE, GAP_PENALTY, temporarily_seq, temporarily_seq1);

        if (i != Length) {
            temporarily_seq = temporarily_seq + minimizers[i].sequence;
            temporarily_seq1 = temporarily_seq1 + minimizers1[i].sequence;
            AlignSeq += temporarily_seq;
            AlignSeq1 += temporarily_seq1;

        }
        else if (i == Length) {
            AlignSeq += temporarily_seq;
            AlignSeq1 += temporarily_seq1;
        }

        seq = "", seq = "", temporarily_seq = "", temporarily_seq1 = "";
    }
}

// compute consistency according to alignment result
double calculate_consistency(const string& seq1, const string& seq2) {
	// number of matching characters
	int64_t matching_chars = 0;
	size_t i, seq_length = max(seq1.size(), seq2.size());

	for (i = 0; i < seq_length; ++i) {
		if (i < seq1.size() && i < seq2.size() && seq1[i] == seq2[i]) {
			matching_chars ++;
			//cout << "matching_chars:" << matching_chars << endl;
		}
	}

	// relief factor
	double relief = 0;

	// compute the number of consecutive gap
	int32_t consecutive_gap_count = 0;
	//int max_consecutive_gap = 0;
	for (i = 0; i <= seq_length; ++i) {
		if (i < seq1.length() and i < seq2.length() and (seq1[i] == '-' || seq2[i] == '-')) {
			consecutive_gap_count++;
			//cout << consecutive_gap_count << endl;
		}else {
			if (consecutive_gap_count >= 3) {
				 relief += (consecutive_gap_count * RELIEF_FACTOR);
				 //cout << "relief:" << relief << endl;
			}
			// reset if there is no consecutive gap
			consecutive_gap_count = 0;
		}
	}

	// compute consistency and return
	return (double)(matching_chars + relief) / seq_length;
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

void destroyResultData(vector<vector<SV_item*>> &result){
	vector<vector<SV_item*>>::iterator it;
	for(it=result.begin(); it!=result.end(); it++)
		destroyData(*it);
	vector<vector<SV_item*>>().swap(result);
}

// SV number stat operation for TRA
void SVNumStatTraOp(string &user_file, string &benchmark_file, string &dirname){
	vector<SV_item*> user_data, benchmark_data, long_sv_data;
	string description_str_long, file_prefix;

	// non-TRA variants
	user_data = loadDataTra(user_file);
	benchmark_data = loadDataTra(benchmark_file);

	cout << "TRA: user data size: " << user_data.size() << endl;
	cout << "TRA: benchmark data size: " << benchmark_data.size() << endl;

	outStatScreenFile << "TRA: user data size: " << user_data.size() << endl;
	outStatScreenFile << "TRA: benchmark data size: " << benchmark_data.size() << endl;

	file_prefix = dirname + "num_stat_TRA_BND";
	computeNumStatTra(user_data, benchmark_data, file_prefix);

	destroyData(user_data);
	destroyData(benchmark_data);
}

void computeNumStatTra(vector<SV_item*> &user_data, vector<SV_item*> &benchmark_data, string &file_prefix){
	vector<Breakpoint_t*> bp_vec_user, bp_vec_benchmark;
	string filename_bp_user, filename_bp_benchmark;

	filename_bp_user = file_prefix + "_user";
	filename_bp_benchmark = file_prefix + "_benchmark";

	// compute intersection
	computeOverlapTra(user_data, benchmark_data);

	// construct breakpoint vectors
	bp_vec_user = constructBPVec(user_data);
	bp_vec_benchmark = constructBPVec(benchmark_data);

	// output breakpoints to file
	outputBP2File(filename_bp_user, bp_vec_user);
	outputBP2File(filename_bp_benchmark, bp_vec_benchmark);

	// compute statistics
	computeBPNumStatTra(bp_vec_user, bp_vec_benchmark, file_prefix);

	// destroy data
	destroyBPData(bp_vec_user);
	destroyBPData(bp_vec_benchmark);
}

void computeOverlapTra(vector<SV_item*> &user_data, vector<SV_item*> &benchmark_data){
	SV_item *item1, *item2;
	size_t i, j, k, overlap_type;
	vector<size_t> overlap_type_vec;
	vector<Breakpoint_t*> bp_vec_user, bp_vec_benchmark;

	for(i=0; i<user_data.size(); i++) {
		user_data.at(i)->overlapped = false;
		for(j=0; j<4; j++) user_data.at(i)->traOverlappedArr[j] = false;
	}
	for(i=0; i<benchmark_data.size(); i++) {
		benchmark_data.at(i)->overlapped = false;
		for(j=0; j<4; j++) benchmark_data.at(i)->traOverlappedArr[j] = false;
	}

	for(i=0; i<user_data.size(); i++){
		item1 = user_data.at(i);
		for(j=0; j<benchmark_data.size(); j++){
			item2 = benchmark_data.at(j);
			overlap_type_vec = computeOverlapType(item1, item2);
			for(k=0; k<overlap_type_vec.size(); k++){
				overlap_type = overlap_type_vec.at(k);
				switch(overlap_type){
					case NO_OVERLAP: break;
					case WHOLE_REG_OVERLAP: item1->overlapped = true; item2->overlapped = true; break;

					case ONE_OVERLAP_ONE_BP: item1->traOverlappedArr[0] = item2->traOverlappedArr[0] = true; break;
					case ONE_OVERLAP_TWO_BP: item1->traOverlappedArr[0] = item2->traOverlappedArr[1] = true; break;
					case ONE_OVERLAP_THREE_BP: item1->traOverlappedArr[0] = item2->traOverlappedArr[2] = true; break;
					case ONE_OVERLAP_FOUR_BP: item1->traOverlappedArr[0] = item2->traOverlappedArr[3] = true; break;

					case TWO_OVERLAP_ONE_BP: item1->traOverlappedArr[1] = item2->traOverlappedArr[0] = true; break;
					case TWO_OVERLAP_TWO_BP: item1->traOverlappedArr[1] = item2->traOverlappedArr[1] = true; break;
					case TWO_OVERLAP_THREE_BP: item1->traOverlappedArr[1] = item2->traOverlappedArr[2] = true; break;
					case TWO_OVERLAP_FOUR_BP: item1->traOverlappedArr[1] = item2->traOverlappedArr[3] = true; break;

					case THREE_OVERLAP_ONE_BP: item1->traOverlappedArr[2] = item2->traOverlappedArr[0] = true; break;
					case THREE_OVERLAP_TWO_BP: item1->traOverlappedArr[2] = item2->traOverlappedArr[1] = true; break;
					case THREE_OVERLAP_THREE_BP: item1->traOverlappedArr[2] = item2->traOverlappedArr[2] = true; break;
					case THREE_OVERLAP_FOUR_BP: item1->traOverlappedArr[2] = item2->traOverlappedArr[3] = true; break;

					case FOUR_OVERLAP_ONE_BP: item1->traOverlappedArr[3] = item2->traOverlappedArr[0] = true; break;
					case FOUR_OVERLAP_TWO_BP: item1->traOverlappedArr[3] = item2->traOverlappedArr[1] = true; break;
					case FOUR_OVERLAP_THREE_BP: item1->traOverlappedArr[3] = item2->traOverlappedArr[2] = true; break;
					case FOUR_OVERLAP_FOUR_BP: item1->traOverlappedArr[3] = item2->traOverlappedArr[3] = true; break;

					default: cerr << "line=" << __LINE__ << ", overlap_type=" << overlap_type << ", error!" << endl; exit(1);
				}
			}
		}
	}
}

// cpmpute breakpoint statistics for TRA
void computeBPNumStatTra(vector<Breakpoint_t*> &bp_vec_user, vector<Breakpoint_t*> &bp_vec_benchmark, string &file_prefix){
	size_t i, total_bp_num_user, overlapped_bp_num_user, total_bp_num_benchmark, overlapped_bp_num_benchmark;
	vector<Breakpoint_t*> bp_vec_intersect_user, bp_vec_intersect_benchmark, bp_vec_private_user, bp_vec_private_benchmark;
	Breakpoint_t *bp_item, *bp_item_tmp;

	int32_t TP_user, TP_benchmark, FP, FN;
	float recall, precision_user, F1_score_user, precision_benchmark, F1_score_benchmark;
	string filename_intersect_user, filename_intersect_benchmark, filename_private_user, filename_private_benchmark;

	filename_intersect_user = file_prefix + "_intersect_user";
	filename_intersect_benchmark = file_prefix + "_intersect_benchmark";
	filename_private_user = file_prefix + "_private_user";
	filename_private_benchmark = file_prefix + "_private_benchmark";

	total_bp_num_user = bp_vec_user.size();
	total_bp_num_benchmark = bp_vec_benchmark.size();

	// construct intersect vector and private vectors
	overlapped_bp_num_user = 0;
	for(i=0; i<bp_vec_user.size(); i++){
		bp_item = bp_vec_user.at(i);
		bp_item_tmp = allocateBPItem(bp_item->chrname, bp_item->bp_loc);
		if(bp_item->overlappedFlag){
			bp_vec_intersect_user.push_back(bp_item_tmp);
			overlapped_bp_num_user ++;
		}else // user-private
			bp_vec_private_user.push_back(bp_item_tmp);
	}
	overlapped_bp_num_benchmark = 0;
	for(i=0; i<bp_vec_benchmark.size(); i++){
		bp_item = bp_vec_benchmark.at(i);
		bp_item_tmp = allocateBPItem(bp_item->chrname, bp_item->bp_loc);
		if(bp_item->overlappedFlag) {
			bp_vec_intersect_benchmark.push_back(bp_item_tmp);
			overlapped_bp_num_benchmark ++;
		}else // benchmark-private
			bp_vec_private_benchmark.push_back(bp_item_tmp);
	}

	// output breakpoints to file
	outputBP2File(filename_intersect_user, bp_vec_intersect_user);
	outputBP2File(filename_intersect_benchmark, bp_vec_intersect_benchmark);
	outputBP2File(filename_private_user, bp_vec_private_user);
	outputBP2File(filename_private_benchmark, bp_vec_private_benchmark);

	// compute TRA breakpoint statistics
	TP_user = overlapped_bp_num_user;
	TP_benchmark = overlapped_bp_num_benchmark;
	FP = total_bp_num_user - overlapped_bp_num_user;
	FN = total_bp_num_benchmark - overlapped_bp_num_benchmark;

	if(total_bp_num_benchmark>0) recall = (float)TP_benchmark / total_bp_num_benchmark;
	else recall = 0;
	if(total_bp_num_user>0) precision_user = (float)TP_user / total_bp_num_user;
	else precision_user = 0;
	if(total_bp_num_user>0) precision_benchmark = (float)TP_benchmark / total_bp_num_user;
	else precision_benchmark = 0;

	if(recall+precision_user>0) F1_score_user = 2.0 * (recall * precision_user) / (recall + precision_user);
	else F1_score_user = 0;
	if(recall+precision_benchmark>0) F1_score_benchmark = 2.0 * (recall * precision_benchmark) / (recall + precision_benchmark);
	else F1_score_benchmark = 0;

	cout << "Total user breakpoint data size: " << total_bp_num_user << endl;
	cout << "Total benchmark breakpoint data size: " << total_bp_num_benchmark << endl;
	cout << "User intersection breakpoint data size:" << overlapped_bp_num_user << endl;
	cout << "Benchmark intersection breakpoint data size:" << overlapped_bp_num_benchmark << endl;
	cout << "User private breakpoint data size:" << bp_vec_private_user.size() << endl;
	cout << "Benchmark private breakpoint data size:" << bp_vec_private_benchmark.size() << endl;
	cout << "TP_user=" << TP_user << ", TP_benchmark=" << TP_benchmark << ", FP=" << FP << ", FN=" << FN << endl;
	cout << "Recall=" << recall << endl;
	cout << "precision_user=" << precision_user << ", precision_benchmark=" << precision_benchmark << endl;
	cout << "F1 score_user=" << F1_score_user << ", F1 score_benchmark=" << F1_score_benchmark << endl;

	outStatScreenFile << "Total user breakpoint data size: " << total_bp_num_user << endl;
	outStatScreenFile << "Total benchmark breakpoint data size: " << total_bp_num_benchmark << endl;
	outStatScreenFile << "User intersection breakpoint data size:" << overlapped_bp_num_user << endl;
	outStatScreenFile << "Benchmark intersection breakpoint data size:" << overlapped_bp_num_benchmark << endl;
	outStatScreenFile << "User private breakpoint data size:" << bp_vec_private_user.size() << endl;
	outStatScreenFile << "Benchmark private breakpoint data size:" << bp_vec_private_benchmark.size() << endl;
	outStatScreenFile << "TP_user=" << TP_user << ", TP_benchmark=" << TP_benchmark << ", FP=" << FP << ", FN=" << FN << endl;
	outStatScreenFile << "Recall=" << recall << endl;
	outStatScreenFile << "precision_user=" << precision_user << ", precision_benchmark=" << precision_benchmark << endl;
	outStatScreenFile <<"F1 score_user=" << F1_score_user << ", F1 score_benchmark=" << F1_score_benchmark << endl;

	destroyBPData(bp_vec_intersect_user);
	destroyBPData(bp_vec_intersect_benchmark);
	destroyBPData(bp_vec_private_user);
	destroyBPData(bp_vec_private_benchmark);
}

Breakpoint_t *allocateBPItem(string &chrname, size_t bp_loc){
	Breakpoint_t *bp_item = new Breakpoint_t();
	bp_item->chrname = chrname;
	bp_item->bp_loc = bp_loc;
	bp_item->overlappedFlag = false;
	return bp_item;
}

// destory breakpoint vector
void destroyBPData(vector<Breakpoint_t*> &bp_vec){
	for(size_t i=0; i<bp_vec.size(); i++) delete bp_vec.at(i);
	vector<Breakpoint_t*>().swap(bp_vec);
}

// construct breakpoint vector according to overlapped SV items
vector<Breakpoint_t*> constructBPVec(vector<SV_item*> &sv_data){
	vector<Breakpoint_t*> bp_vec;
	size_t i;
	SV_item *item;
	Breakpoint_t *bp_item;

	for(i=0; i<sv_data.size(); i++){
		item = sv_data.at(i);
		if(item->sv_type==VAR_TRA or item->sv_type==VAR_BND){
			if(item->chrname.size()>0){
//				if(item->startPos!=item->endPos){ // two breakpoints
					if(item->startPos!=0){
						bp_item = allocateBPItem(item->chrname, item->startPos);
						if(item->traOverlappedArr[0]) bp_item->overlappedFlag = true;
						bp_vec.push_back(bp_item);
					}
					if(item->endPos!=0){
						bp_item = allocateBPItem(item->chrname, item->endPos);
						if(item->traOverlappedArr[1]) bp_item->overlappedFlag = true;
						bp_vec.push_back(bp_item);
					}
//				}else{ // one breakpoint
//					bp_item = allocateBPItem(item->chrname, item->startPos);
//					if(item->traOverlappedArr[0]) bp_item->overlappedFlag = true;
//					bp_vec.push_back(bp_item);
//				}
			}
			if(item->chrname2.size()>0){
//				if(item->startPos2!=item->endPos2){ // two breakpoints
					if(item->startPos2!=0){
						bp_item = allocateBPItem(item->chrname2, item->startPos2);
						if(item->traOverlappedArr[2]) bp_item->overlappedFlag = true;
						bp_vec.push_back(bp_item);
					}
					if(item->endPos2!=0){
						bp_item = allocateBPItem(item->chrname2, item->endPos2);
						if(item->traOverlappedArr[3]) bp_item->overlappedFlag = true;
						bp_vec.push_back(bp_item);
					}
//				}else{ // one breakpoint
//					bp_item = allocateBPItem(item->chrname2, item->startPos2);
//					if(item->traOverlappedArr[2]) bp_item->overlappedFlag = true;
//					bp_vec.push_back(bp_item);
//				}
			}
		}
	}

	return bp_vec;
}

// output breakpoints to file
void outputBP2File(string &filename, vector<Breakpoint_t*> &bp_vec){
	ofstream outfile;
	string line;
	Breakpoint_t *bp_item;

	outfile.open(filename);
	if(!outfile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << filename << endl;
		exit(1);
	}

	line = "#chr\trefPos";
	outfile << line << endl;

	for(size_t i=0; i<bp_vec.size(); i++){
		bp_item = bp_vec.at(i);
		line = bp_item->chrname + "\t" + to_string(bp_item->bp_loc);
		outfile << line << endl;
	}

	outfile.close();
}
