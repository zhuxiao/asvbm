#ifndef SRC_TYPE_NUM_STAT_H_
#define SRC_TYPE_NUM_STAT_H_

#include <iostream>
#include <cstring>
#include <vector>
#include <sys/stat.h>

#include "structure.h"
#include "dataLoader.h"

using namespace std;

void SVTypeNumStat(string &sv_file1, string &sv_file2, int32_t max_valid_reg_thres);
void SVTypeNumStatOp(string &sv_file1, string &sv_file2, vector<size_t> &size_div_vec, int32_t max_valid_reg_thres);
vector<vector<SV_item*>> sizeDivideSV(vector<SV_item*> &sv_data, vector<size_t> &size_div_vec);
void computeTypeNumStat(vector<vector<SV_item*>> &divided_vec1, vector<vector<SV_item*>> &divided_vec2);

void destroySizeDividedData(vector< vector<SV_item*> > &divided_vec);

extern SV_item* itemdup(SV_item* item);
extern void computeNumStat(vector<SV_item*> &sv_data1, vector<SV_item*> &sv_data2, string &file_prefix);
extern void destroyResultData(vector<vector<SV_item*>> &result);

#endif /* SRC_TYPE_NUM_STAT_H_ */
