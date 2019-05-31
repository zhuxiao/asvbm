#ifndef SRC_REF_REG_SIZE_STAT_H_
#define SRC_REF_REG_SIZE_STAT_H_

#include <iostream>
#include <cstring>
#include <vector>
#include <sys/stat.h>

#include "structure.h"
#include "dataLoader.h"

using namespace std;

#define SV_SIZE_ARR_SIZE	2000


void refRegSizeStat(string &user_file, string &benchmark_file, int32_t max_valid_reg_thres);
void refRegSizeStatOp(string &refRegSizeFinename, string &sv_file, int32_t max_valid_reg_thres, string &dirname);

#endif /* SRC_REF_REG_SIZE_STAT_H_ */
