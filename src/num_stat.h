#ifndef SRC_NUM_STAT_H_
#define SRC_NUM_STAT_H_

#include <iostream>
#include <cstring>
#include <vector>
#include <limits.h>
#include <sys/stat.h>
#include <set>
#include <algorithm>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>

#include "constants.h"
#include "structure.h"
#include "util.h"

using namespace std;


void SVNumStat(string &user_file, string &benchmark_file, int32_t max_valid_reg_thres, string &outputPathname);
void SVNumStatOp(string &user_file, string &benchmark_file, int32_t max_valid_reg_thres, string &dirname);
void computeNumStat(vector<SV_item*> &sv_data1, vector<SV_item*> &sv_data2, string &file_prefix);

void computeLenStat(vector<SV_item*> &data, string &description_str);

vector<vector<SV_item*>> intersect(vector<SV_item*> &data1, vector<SV_item*> &data2);
vector<vector<SV_item*>> constructSubsetByChr(vector<SV_item*> &user_data, vector<SV_item*> &benchmark_data);
set<string> getChrnames(vector<SV_item*> &dataset);
set<string> getChrUnion(set<string> &chrname_set1, set<string> &chrname_set2);
vector<string> sortChrnames(set<string> &chrname_set);
vector<vector<SV_item*>> constructSubsetByChrOp(vector<SV_item*> &user_data, vector<SV_item*> &benchmark_data, vector<string> &chrname_vec);
void sortSubsets(vector<vector<SV_item*>> &subsets);
void* sortSubsetOp(void *arg);
bool sortFunSameChr(const SV_item *item1, const SV_item *item2);
void checkOrder(vector<vector<SV_item*>> &subsets);
vector<vector<SV_item*>> intersectOp(vector<vector<SV_item*>> &subsets);
void* intersectSubset(void *arg);
vector<SV_item*> getItemsByChr(string &chrname, vector<SV_item*> &dataset);
SV_item* itemdup(SV_item* item);
bool IsSameChrname(string &chrname1, string &chrname2);
vector<size_t> computeOverlapType(SV_item* item1, SV_item* item2);
bool isOverlappedPos(size_t startPos1, size_t endPos1, size_t startPos2, size_t endPos2);

void destroyResultData(vector<vector<SV_item*>> &result);

void SVNumStatTraOp(string &user_file, string &benchmark_file, string &dirname);
void computeNumStatTra(vector<SV_item*> &sv_data1, vector<SV_item*> &sv_data2, string &file_prefix);
void computeOverlapTra(vector<SV_item*> &data1, vector<SV_item*> &data2);
void computeBPNumStatTra(vector<Breakpoint_t*> &bp_vec_user, vector<Breakpoint_t*> &bp_vec_benchmark, string &file_prefix);
Breakpoint_t *allocateBPItem(string &chrname, size_t bp_loc);
void destroyBPData(vector<Breakpoint_t*> &bp_vec);
vector<Breakpoint_t*> constructBPVec(vector<SV_item*> &sv_data);
void outputBP2File(string &filename, vector<Breakpoint_t*> &bp_vec);

#endif /* SRC_NUM_STAT_H_ */
