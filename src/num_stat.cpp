#include "dataLoader.h"
#include "num_stat.h"
#include "extvab.h"


void SVNumStat(string &user_file, string &benchmark_file, int32_t max_valid_reg_thres, string &outputPathname){
	numStatDirname = outputPathname + numStatDirname;
	mkdir(numStatDirname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	// non-TRA
	if(max_valid_reg_thres>0){
		cout << ">>>>>>>>> Before filtering long SV regions: <<<<<<<<<" << endl;
		outStatScreenFile << ">>>>>>>>> Before filtering long SV regions: <<<<<<<<<" << endl;
	}
	SVNumStatOp(user_file, benchmark_file, 0, numStatDirname);

	if(max_valid_reg_thres>0){
		cout << "\n>>>>>>>>> After filtering long SV regions: <<<<<<<<<" << endl;
		outStatScreenFile << "\n>>>>>>>>> After filtering long SV regions: <<<<<<<<<" << endl;
		SVNumStatOp(user_file, benchmark_file, max_valid_reg_thres, numStatDirname);
	}

	// TRA
	cout << "\n>>>>>>>>> TRA breakpoint statistics: <<<<<<<<<" << endl;
	outStatScreenFile << "\n>>>>>>>>> TRA breakpoint statistics: <<<<<<<<<" << endl;
	SVNumStatTraOp(user_file, benchmark_file, numStatDirname);
}

// SV number stat operation
void SVNumStatOp(string &user_file, string &benchmark_file, int32_t max_valid_reg_thres, string &dirname){
	vector<SV_item*> user_data, benchmark_data, long_sv_data;
	string description_str_long, file_prefix, longFilename_tmp;

	longFilename_tmp = dirname + longSVFilename;

	// non-TRA variants
	user_data = loadDataWithoutTra(user_file);
	benchmark_data = loadDataWithoutTra(benchmark_file);

	if(max_valid_reg_thres>0){
		long_sv_data = getLongSVReg(user_data, max_valid_reg_thres);
		output2File(longFilename_tmp, long_sv_data, outStatScreenFile);
		cout << "non-TRA: user data size: " << user_data.size() << endl;
		cout << "non-TRA: benchmark data size: " << benchmark_data.size() << endl;
		cout << "non-TRA: long_sv_data.size: " << long_sv_data.size() << endl;
		outStatScreenFile << "non-TRA: user data size: " << user_data.size() << endl;
		outStatScreenFile << "non-TRA: benchmark data size: " << benchmark_data.size() << endl;
		outStatScreenFile << "non-TRA: long_sv_data.size: " << long_sv_data.size() << endl;
		description_str_long = "Size statistics for long SV regions: ";
		computeLenStat(long_sv_data, description_str_long);
	}else{
		cout << "non-TRA: user data size: " << user_data.size() << endl;
		cout << "non-TRA: benchmark data size: " << benchmark_data.size() << endl;
		outStatScreenFile << "non-TRA: user data size: " << user_data.size() << endl;
		outStatScreenFile << "non-TRA: benchmark data size: " << benchmark_data.size() << endl;
	}

	file_prefix = dirname + "num_stat";
	computeNumStat(user_data, benchmark_data, file_prefix);

	destroyData(user_data);
	destroyData(benchmark_data);
	if(max_valid_reg_thres>0) destroyData(long_sv_data);
}

void computeNumStat(vector<SV_item*> &user_data, vector<SV_item*> &benchmark_data, string &file_prefix){
	vector< vector<SV_item*> > result;
	int32_t TP, FP, FN, positive_num_called;
	float recall, precision, F1_score, sv_num_per_reg, percent;
	string filename_intersect_user, filename_intersect_benchmark, filename_private_user, filename_private_benchmark, out_str;

	filename_intersect_user = file_prefix + "_intersect_user";
	filename_intersect_benchmark = file_prefix + "_intersect_benchmark";
	filename_private_user = file_prefix + "_private_user";
	filename_private_benchmark = file_prefix + "_private_benchmark";

	// compute intersection
	result = intersect(user_data, benchmark_data);

	if(user_data.size()>0) percent = (double)result.at(0).size() / user_data.size();
	else percent = 0;
	out_str = "user intersection data size: " + to_string(result.at(0).size()) + ", percent: " + to_string(percent);
	cout << out_str << endl;
	outStatScreenFile << out_str << endl;

	if(benchmark_data.size()>0) percent = (double)result.at(1).size() / benchmark_data.size();
	else percent = 0;
	out_str = "benchmark intersection data size: " + to_string(result.at(1).size()) + ", percent: " + to_string(percent);
	cout << out_str << endl;
	outStatScreenFile << out_str << endl;

	if(user_data.size()>0) percent = (double)result.at(1).size()/user_data.size();
	else percent = 0;
	out_str = "user private data size: " + to_string(result.at(1).size()) + ", percent: " + to_string(percent);
	cout << out_str << endl;
	outStatScreenFile << out_str << endl;

	if(benchmark_data.size()>0) percent = (double)result.at(2).size()/benchmark_data.size();
	else percent = 0;
	out_str = "benchmark private data size: " + to_string(result.at(2).size()) + ", percent: " + to_string(percent);
	cout << out_str << endl;
	outStatScreenFile << out_str << endl;

	// output to file
	output2File(filename_intersect_user, result.at(0), outStatScreenFile);
	output2File(filename_intersect_benchmark, result.at(1), outStatScreenFile);
	output2File(filename_private_user, result.at(2), outStatScreenFile);
	output2File(filename_private_benchmark, result.at(3), outStatScreenFile);

	//TP = user_data.size() - result.at(1).size();
	TP = result.at(1).size();
	FP = result.at(2).size();
	FN = result.at(3).size();
	//positive_num_called = user_data.size() - result.at(2).size();
	positive_num_called = result.at(0).size();
	if(benchmark_data.size()>0) recall = (float)TP / benchmark_data.size();
	else recall = 0;
	//precision = (float)TP / user_data.size();
	if(TP+result.at(2).size()>0) precision = (float)TP / (TP + result.at(2).size());
	else precision = 0;

	if(recall+precision>0) F1_score = 2.0 * (recall * precision) / (recall + precision);
	else F1_score = 0;

	if(positive_num_called>0) sv_num_per_reg = (float)TP / positive_num_called;
	else sv_num_per_reg = 0;

	cout << "TP=" << TP << ", FP=" << FP << ", FN=" << FN << endl;
	cout << "Identified SV regions: " << positive_num_called << endl;
	cout << "Recall=" << recall << ", precision=" << precision << ", F1 score=" << F1_score << ", sv_num_per_reg=" << sv_num_per_reg << endl;

	outStatScreenFile << "TP=" << TP << ", FP=" << FP << ", FN=" << FN << endl;
	outStatScreenFile << "Identified SV regions: " << positive_num_called << endl;
	outStatScreenFile << "Recall=" << recall << ", precision=" << precision << ", F1 score=" << F1_score << ", sv_num_per_reg=" << sv_num_per_reg << endl;

	destroyResultData(result);
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

vector<vector<SV_item*>> intersect(vector<SV_item*> &user_data, vector<SV_item*> &benchmark_data){
	vector< vector<SV_item*> > result;
	vector<SV_item*> intersect_vec_user, intersect_vec_benchmark, private_vec_user, private_vec_benchmark;
	SV_item *item1, *item2, *item;
	size_t i, j;
	vector<size_t> overlap_type_vec;

	for(i=0; i<user_data.size(); i++) user_data.at(i)->overlapped = false;
	for(i=0; i<benchmark_data.size(); i++) benchmark_data.at(i)->overlapped = false;

	for(i=0; i<user_data.size(); i++){
		item1 = user_data.at(i);
		for(j=0; j<benchmark_data.size(); j++){
			item2 = benchmark_data.at(j);
			overlap_type_vec = computeOverlapType(item1, item2);
			if(overlap_type_vec.size()>0 and overlap_type_vec.at(0)!=NO_OVERLAP){
				item1->overlapped = true;
				item2->overlapped = true;
			}
		}
	}

	for(i=0; i<benchmark_data.size(); i++){
		item2 = benchmark_data.at(i);
		item = itemdup(item2);
		if(item2->overlapped)
			intersect_vec_benchmark.push_back(item);
		else
			private_vec_benchmark.push_back(item);
	}

	for(i=0; i<user_data.size(); i++){
		item1 = user_data.at(i);
		item = itemdup(item1);
		if(!item1->overlapped)
			private_vec_user.push_back(item);
		else
			intersect_vec_user.push_back(item);
	}

	result.push_back(intersect_vec_user);
	result.push_back(intersect_vec_benchmark);
	result.push_back(private_vec_user);
	result.push_back(private_vec_benchmark);

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
	for(size_t i=0; i<4; i++) item_new->traOverlappedArr[i] = item->traOverlappedArr[i];
	return item_new;
}

vector<size_t> computeOverlapType(SV_item* item1, SV_item* item2){
	bool flag = false, flag_tmp = false;
	int32_t startPos1, endPos1, startPos2, endPos2;
	vector<size_t> overlap_type_vec;
	size_t overlap_type;

	if(item1->sv_type!=VAR_TRA and item1->sv_type!=VAR_BND){ // INS, DEL, DUP, INV
		if(item1->chrname.size()>0 and item1->chrname.compare(item2->chrname)==0){
			startPos1 = item1->startPos - extendSize;
			if(startPos1<1) startPos1 = 1;
			endPos1 = item1->endPos + extendSize;
			startPos2 = item2->startPos - extendSize;
			if(startPos2<1) startPos2 = 1;
			endPos2 = item2->endPos + extendSize;

			flag = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
		}

		if(flag==false and (item2->sv_type==VAR_TRA or item2->sv_type==VAR_BND)){
			if(item1->chrname.size()>0 and item1->chrname.compare(item2->chrname2)==0){
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
			if(item1->chrname.compare(item2->chrname)==0){
				// 1: one_overlap_one_bp
				startPos1 = item1->startPos - extendSize;
				if(startPos1<1) startPos1 = 1;
				endPos1 = item1->startPos + extendSize;
				startPos2 = item2->startPos - extendSize;
				if(startPos2<1) startPos2 = 1;
				endPos2 = item2->startPos + extendSize;

				flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
				if(flag_tmp) overlap_type = ONE_OVERLAP_ONE_BP;
				else if(item2->startPos!=item2->endPos){
					// one_overlap_two_bp
					startPos2 = item2->endPos - extendSize;
					if(startPos2<1) startPos2 = 1;
					endPos2 = item2->endPos + extendSize;

					flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
					if(flag_tmp) overlap_type = ONE_OVERLAP_TWO_BP;
				}
				if(flag_tmp) overlap_type_vec.push_back(overlap_type);

				if(item1->startPos!=item1->endPos){
					// 2: two_overlap_one_bp
					startPos1 = item1->endPos - extendSize;
					if(startPos1<1) startPos1 = 1;
					endPos1 = item1->endPos + extendSize;
					startPos2 = item2->startPos - extendSize;
					if(startPos2<1) startPos2 = 1;
					endPos2 = item2->startPos + extendSize;

					flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
					if(flag_tmp) overlap_type = TWO_OVERLAP_ONE_BP;
					else if(item2->startPos!=item2->endPos){
						// two_overlap_two_bp
						startPos2 = item2->endPos - extendSize;
						if(startPos2<1) startPos2 = 1;
						endPos2 = item2->endPos + extendSize;

						flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
						if(flag_tmp) overlap_type = TWO_OVERLAP_TWO_BP;
					}
					if(flag_tmp) overlap_type_vec.push_back(overlap_type);
				}
			}else if(item1->chrname.compare(item2->chrname2)==0){
				// one_overlap_three_bp
				startPos1 = item1->startPos - extendSize;
				if(startPos1<1) startPos1 = 1;
				endPos1 = item1->startPos + extendSize;
				startPos2 = item2->startPos2 - extendSize;
				if(startPos2<1) startPos2 = 1;
				endPos2 = item2->startPos2 + extendSize;

				flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
				if(flag_tmp) overlap_type = ONE_OVERLAP_THREE_BP;
				else if(item2->startPos2!=item2->endPos2){
					// one_overlap_four_bp
					startPos2 = item2->endPos2 - extendSize;
					if(startPos2<1) startPos2 = 1;
					endPos2 = item2->endPos2 + extendSize;

					flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
					if(flag_tmp) overlap_type = ONE_OVERLAP_FOUR_BP;
				}
				if(flag_tmp) overlap_type_vec.push_back(overlap_type);

				if(item1->startPos!=item1->endPos){
					// two_overlap_three_bp
					startPos1 = item1->endPos - extendSize;
					if(startPos1<1) startPos1 = 1;
					endPos1 = item1->endPos + extendSize;
					startPos2 = item2->startPos2 - extendSize;
					if(startPos2<1) startPos2 = 1;
					endPos2 = item2->startPos2 + extendSize;

					flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
					if(flag_tmp) overlap_type = TWO_OVERLAP_THREE_BP;
					else if(item2->startPos2!=item2->endPos2) { // two_overlap_four_bp
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
			if(item1->chrname2.compare(item2->chrname)==0){
				// 3: three_overlap_one_bp
				startPos1 = item1->startPos2 - extendSize;
				if(startPos1<1) startPos1 = 1;
				endPos1 = item1->startPos2 + extendSize;
				startPos2 = item2->startPos - extendSize;
				if(startPos2<1) startPos2 = 1;
				endPos2 = item2->startPos + extendSize;

				flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
				if(flag_tmp) overlap_type = THREE_OVERLAP_ONE_BP;
				else if(item2->startPos!=item2->endPos){
					// three_overlap_two_bp
					startPos2 = item2->endPos - extendSize;
					if(startPos2<1) startPos2 = 1;
					endPos2 = item2->endPos + extendSize;

					flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
					if(flag_tmp) overlap_type = THREE_OVERLAP_TWO_BP;
				}
				if(flag_tmp) overlap_type_vec.push_back(overlap_type);

				if(item1->startPos2!=item1->endPos2){
					// 4: four_overlap_one_bp
					startPos1 = item1->endPos2 - extendSize;
					if(startPos1<1) startPos1 = 1;
					endPos1 = item1->endPos2 + extendSize;
					startPos2 = item2->startPos - extendSize;
					if(startPos2<1) startPos2 = 1;
					endPos2 = item2->startPos + extendSize;

					flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
					if(flag_tmp) overlap_type = FOUR_OVERLAP_ONE_BP;
					else if(item2->startPos!=item2->endPos){
						// four_overlap_two_bp
						startPos2 = item2->endPos - extendSize;
						if(startPos2<1) startPos2 = 1;
						endPos2 = item2->endPos + extendSize;

						flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
						if(flag_tmp) overlap_type = FOUR_OVERLAP_TWO_BP;
					}
					if(flag_tmp) overlap_type_vec.push_back(overlap_type);
				}
			}else if(item1->chrname2.compare(item2->chrname2)==0){
				// three_overlap_three_bp
				startPos1 = item1->startPos2 - extendSize;
				if(startPos1<1) startPos1 = 1;
				endPos1 = item1->startPos2 + extendSize;
				startPos2 = item2->startPos2 - extendSize;
				if(startPos2<1) startPos2 = 1;
				endPos2 = item2->startPos2 + extendSize;

				flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
				if(flag_tmp) overlap_type = THREE_OVERLAP_THREE_BP;
				else if(item2->startPos2!=item2->endPos2) {
					// three_overlap_four_bp
					startPos2 = item2->endPos2 - extendSize;
					if(startPos2<1) startPos2 = 1;
					endPos2 = item2->endPos2 + extendSize;

					flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
					if(flag_tmp) overlap_type = THREE_OVERLAP_FOUR_BP;
				}
				if(flag_tmp) overlap_type_vec.push_back(overlap_type);

				if(item1->startPos2!=item1->endPos2){
					// four_overlap_three_bp
					startPos1 = item1->endPos2 - extendSize;
					if(startPos1<1) startPos1 = 1;
					endPos1 = item1->endPos2 + extendSize;
					startPos2 = item2->startPos2 - extendSize;
					if(startPos2<1) startPos2 = 1;
					endPos2 = item2->startPos2 + extendSize;

					flag_tmp = isOverlappedPos(startPos1, endPos1, startPos2, endPos2);
					if(flag_tmp) overlap_type = FOUR_OVERLAP_THREE_BP;
					else if(item2->startPos2!=item2->endPos2) { // four_overlap_four_bp
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

	// compute intersection
	computeOverlapTra(user_data, benchmark_data);

	// construct breakpoint vectors
	bp_vec_user = constructBPVec(user_data);
	bp_vec_benchmark = constructBPVec(benchmark_data);

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
	vector<Breakpoint_t*> bp_vec_intersect, bp_vec_private_user, bp_vec_private_benchmark;
	Breakpoint_t *bp_item, *bp_item_tmp;

	int32_t TP, FP, FN, positive_num_called;
	float recall, precision, F1_score;
	string filename_intersect, filename_private_user, filename_private_benchmark;

	filename_intersect = file_prefix + "_BP_intersect";
	filename_private_user = file_prefix + "_BP_private_user";
	filename_private_benchmark = file_prefix + "_BP_private_benchmark";

	total_bp_num_user = bp_vec_user.size();
	total_bp_num_benchmark = bp_vec_benchmark.size();

	// construct intersect vector and private vectors
	overlapped_bp_num_user = 0;
	for(i=0; i<bp_vec_user.size(); i++){
		bp_item = bp_vec_user.at(i);
		if(bp_item->overlappedFlag) overlapped_bp_num_user ++;
		else{ // user-private
			bp_item_tmp = allocateBPItem(bp_item->chrname, bp_item->bp_loc);
			bp_vec_private_user.push_back(bp_item_tmp);
		}
	}
	overlapped_bp_num_benchmark = 0;
	for(i=0; i<bp_vec_benchmark.size(); i++){
		bp_item = bp_vec_benchmark.at(i);
		if(bp_item->overlappedFlag) {
			bp_item_tmp = allocateBPItem(bp_item->chrname, bp_item->bp_loc);
			bp_vec_intersect.push_back(bp_item_tmp);
			overlapped_bp_num_benchmark ++;
		}else{ // user-private
			bp_item_tmp = allocateBPItem(bp_item->chrname, bp_item->bp_loc);
			bp_vec_private_benchmark.push_back(bp_item_tmp);
		}
	}

	// output breakpoints to file
	outputBP2File(filename_intersect, bp_vec_intersect);
	outputBP2File(filename_private_user, bp_vec_private_user);
	outputBP2File(filename_private_benchmark, bp_vec_private_benchmark);

	// compute TRA breakpoint statistics
	TP = overlapped_bp_num_benchmark;
	FP = total_bp_num_user - overlapped_bp_num_user;
	FN = total_bp_num_benchmark - overlapped_bp_num_benchmark;
	positive_num_called = overlapped_bp_num_user;
	if(total_bp_num_benchmark>0) recall = (float)TP / total_bp_num_benchmark;
	else recall = 0;
	if(total_bp_num_user>0) precision = (float)TP / total_bp_num_user;
	else precision = 0;

	if(recall+precision>0) F1_score = 2.0 * (recall * precision) / (recall + precision);
	else F1_score = 0;

	cout << "user breakpoint data size: " << total_bp_num_user << endl;
	cout << "benchmark breakpoint data size: " << total_bp_num_benchmark << endl;
	cout << "TP=" << TP << ", FP=" << FP << ", FN=" << FN << endl;
	cout << "Identified SV regions: " << positive_num_called << endl;
	cout << "Recall=" << recall << ", precision=" << precision << ", F1 score=" << F1_score << endl;

	outStatScreenFile << "user breakpoint data size: " << total_bp_num_user << endl;
	outStatScreenFile << "benchmark breakpoint data size: " << total_bp_num_benchmark << endl;
	outStatScreenFile << "TP=" << TP << ", FP=" << FP << ", FN=" << FN << endl;
	outStatScreenFile << "Identified SV regions: " << positive_num_called << endl;
	outStatScreenFile << "Recall=" << recall << ", precision=" << precision << ", F1 score=" << F1_score << endl;

	destroyBPData(bp_vec_intersect);
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
				if(item->startPos!=item->endPos){ // two breakpoints
					bp_item = allocateBPItem(item->chrname, item->startPos);
					if(item->traOverlappedArr[0]) bp_item->overlappedFlag = true;
					bp_vec.push_back(bp_item);

					bp_item = allocateBPItem(item->chrname, item->endPos);
					if(item->traOverlappedArr[1]) bp_item->overlappedFlag = true;
					bp_vec.push_back(bp_item);
				}else{ // one breakpoint
					bp_item = allocateBPItem(item->chrname, item->startPos);
					if(item->traOverlappedArr[0]) bp_item->overlappedFlag = true;
					bp_vec.push_back(bp_item);
				}
			}
			if(item->chrname2.size()>0){
				if(item->startPos2!=item->endPos2){ // two breakpoints
					bp_item = allocateBPItem(item->chrname2, item->startPos2);
					if(item->traOverlappedArr[2]) bp_item->overlappedFlag = true;
					bp_vec.push_back(bp_item);

					bp_item = allocateBPItem(item->chrname2, item->endPos2);
					if(item->traOverlappedArr[3]) bp_item->overlappedFlag = true;
					bp_vec.push_back(bp_item);
				}else{ // one breakpoint
					bp_item = allocateBPItem(item->chrname2, item->startPos2);
					if(item->traOverlappedArr[2]) bp_item->overlappedFlag = true;
					bp_vec.push_back(bp_item);
				}
			}
		}
	}

	return bp_vec;
}

// output breakpoints to file
void outputBP2File(string &filename, vector<Breakpoint_t*> &bp_vec){
	ofstream outfile;
	string data_out;
	Breakpoint_t *bp_item;

	outfile.open(filename);
	if(!outfile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << filename << endl;
		exit(1);
	}

	for(size_t i=0; i<bp_vec.size(); i++){
		bp_item = bp_vec.at(i);
		data_out = bp_item->chrname + "\t" + to_string(bp_item->bp_loc);
		outfile << data_out << endl;
	}

	outfile.close();
}
