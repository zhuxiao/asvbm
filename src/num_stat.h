#ifndef SRC_NUM_STAT_H_
#define SRC_NUM_STAT_H_

#include <iostream>
#include <cstring>
#include <vector>
#include <limits.h>

#include "structure.h"

using namespace std;

#define ADJUST_SIZE		10


void SVNumStat(string &sv_file1, string &sv_file2, int32_t max_valid_reg_thres, string &outputPathname);
void SVNumStatOp(string &sv_file1, string &sv_file2, int32_t max_valid_reg_thres, string &dirname);
void computeNumStat(vector<SV_item*> &sv_data1, vector<SV_item*> &sv_data2, string &file_prefix);

void computeLenStat(vector<SV_item*> &data, string &description_str);

vector<vector<SV_item*>> intersect(vector<SV_item*> &data1, vector<SV_item*> &data2);
SV_item* itemdup(SV_item* item);
bool haveOverlap(SV_item* item1, SV_item* item2);
bool isOverlappedPos(size_t startPos1, size_t endPos1, size_t startPos2, size_t endPos2);

void destroyResultData(vector<vector<SV_item*>> &result);

#endif /* SRC_NUM_STAT_H_ */
