#ifndef SRC_SIZE_NUM_STAT_H_
#define SRC_SIZE_NUM_STAT_H_

#include <iostream>
#include <cstring>
#include <vector>
#include <sys/stat.h>

#include "structure.h"
#include "dataLoader.h"

using namespace std;

void FolderInit();
void SVSizeNumStat(string &user_file, string &benchmark_file, string &ref_file, int32_t max_valid_reg_thres, int32_t min_valid_reg_thres, vector<string> &sv_files1);
void SVSizeNumStatOp(string &user_file, string &benchmark_file, string &ref_file, vector<size_t> &size_div_vec, int32_t max_valid_reg_thres, int32_t min_valid_reg_thres,string &dirname);
vector<vector<SV_item*>> sizeDivideSV(vector<SV_item*> &sv_data, vector<size_t> &size_div_vec);
void computeSizeNumStat(vector<vector<SV_item*>> &divided_vec1, vector<vector<SV_item*>> &divided_vec2, string &dirname, string &ref_file);
void computeSizeNumStatFromFile(vector<SV_item*> user_data, vector<SV_item*> benchmark_data, vector<SV_item*> TP_bench, vector<SV_item*> TP_user, vector<SV_item*> FP, vector<SV_item*> FN, string &dirname);

void destroySizeDividedData(vector< vector<SV_item*> > &divided_vec);
void InformationStorage(size_t num, vector<float> &Data);

extern SV_item* itemdup(SV_item* item);
extern void computeNumStatFromFile(vector<SV_item*> TPbench_data, vector<SV_item*> TPuser_data, vector<SV_item*> FP_data, vector<SV_item*> FN_data, string &file_prefix, int32_t endpos, int Markers);
extern void computeNumStat(vector<SV_item*> &sv_data1, vector<SV_item*> &sv_data2, string &file_prefix, faidx_t *fai, int Markers);
extern void computeNumStatTra(vector<SV_item*> &sv_data1, vector<SV_item*> &sv_data2, string &file_prefix);
extern void destroyResultData(vector<vector<SV_item*>> &result);

#endif /* SRC_SIZE_NUM_STAT_H_ */
