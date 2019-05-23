#include "dataLoader.h"
#include "num_stat.h"
#include "extvab.h"


void SVNumStat(string &sv_file1, string &sv_file2, int32_t max_valid_reg_thres, string &outputPathname){
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
		description_str_long = "Size statistics for long SV regions: ";
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

void destroyResultData(vector<vector<SV_item*>> &result){
	vector<vector<SV_item*>>::iterator it;
	for(it=result.begin(); it!=result.end(); it++)
		destroyData(*it);
	vector<vector<SV_item*>>().swap(result);
}
