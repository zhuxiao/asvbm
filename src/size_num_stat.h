#ifndef SRC_SIZE_NUM_STAT_H_
#define SRC_SIZE_NUM_STAT_H_

#include <iostream>
#include <cstring>
#include <vector>
#include <sys/stat.h>

#include "structure.h"
#include "dataLoader.h"

using namespace std;

void SVSizeNumStat(string &user_file, string &benchmark_file, int32_t max_valid_reg_thres);
void SVSizeNumStatOp(string &user_file, string &benchmark_file, vector<size_t> &size_div_vec, int32_t max_valid_reg_thres, string &dirname);
vector<vector<SV_item*>> sizeDivideSV(vector<SV_item*> &sv_data, vector<size_t> &size_div_vec);
void computeSizeNumStat(vector<vector<SV_item*>> &divided_vec1, vector<vector<SV_item*>> &divided_vec2, string &dirname);

void destroySizeDividedData(vector< vector<SV_item*> > &divided_vec);

extern SV_item* itemdup(SV_item* item);
extern void computeNumStat(vector<SV_item*> &sv_data1, vector<SV_item*> &sv_data2, string &file_prefix);
extern void computeNumStatTra(vector<SV_item*> &sv_data1, vector<SV_item*> &sv_data2, string &file_prefix);
extern void destroyResultData(vector<vector<SV_item*>> &result);

#endif /* SRC_SIZE_NUM_STAT_H_ */
