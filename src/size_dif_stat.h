#ifndef SRC_SIZE_DIF_STAT_H_
#define SRC_SIZE_DIF_STAT_H_

#include <iostream>
#include <cstring>
#include <vector>
#include <limits.h>
#include <sys/stat.h>

#include "constants.h"
#include "structure.h"
#include "dataLoader.h"

using namespace std;

#define MIN_DIF_SIZE		20
#define MAX_DIF_SIZE		1000

#define SUB_DIF_REG_SIZE	50


void SVSizeDifStat(string &user_file, string &benchmark_file, int32_t max_valid_reg_thres);
void SVSizeDifStatOp(string &user_file, string &benchmark_file, int32_t max_valid_reg_thres, string &dirname);
vector<SV_pair*> computeOverlapSVPair(vector<SV_item*> &data1, vector<SV_item*> &data2);
void outputPairDataToFile(string &filename, vector<SV_pair*> &sv_pair_vec);
vector< vector<int32_t> > computeDifStatVec(vector<SV_pair*> &sv_pair_vec);
void outputDifStatToFile(string &svSizeDifStatFilename, vector< vector<int32_t> > &dif_stat_vec);
vector<size_t> computeRatioStatVec(vector<SV_pair*> &sv_pair_vec, vector<double> &ratio_div_vec);
void outputRatioStatToFile(string &svSizeRatioStatFilename_tmp, vector<size_t> &ratio_stat_vec, vector<double> &ratio_div_vec);

void destroyPairData(vector<SV_pair*> &sv_vec);

extern vector<size_t> computeOverlapType(SV_item* item1, SV_item* item2);

#endif /* SRC_SIZE_DIF_STAT_H_ */
