#ifndef SV_STAT_H_
#define SV_STAT_H_

#include <iostream>
#include <cstring>
#include <fstream>
#include <vector>
#include <getopt.h>

#include "structure.h"

#include "ref_reg_size_stat.h"
#include "num_stat.h"
#include "size_dif_stat.h"
#include "type_num_stat.h"

using namespace std;


void showUsage();
void showUsageConvert();
void showUsageStat();
int parseConvert(int argc, char **argv);
int parseStat(int argc, char **argv);
void SVStat(string &user_file, string &benchmark_file);

#endif /* SV_STAT_H_ */
