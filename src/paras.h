#ifndef SRC_PARAS_H_
#define SRC_PARAS_H_

#include <iostream>
#include <cstring>
#include <fstream>
#include <vector>
#include <getopt.h>

#include "structure.h"

#include "convert.h"

#include "ref_reg_size_stat.h"
#include "num_stat.h"
#include "size_dif_stat.h"
#include "type_num_stat.h"
#include "util.h"

using namespace std;

void showUsage();
void showUsageConvert();
void showUsageStat();
int parseConvert(int argc, char **argv);
int parseStat(int argc, char **argv);

void convert(string &infilename, string &outfilename, bool remove_dup_falg, string &sv_format);
void SVStat(string &user_file, string &benchmark_file);

void printConvertParas(string &infilename, string &outfilename, bool remove_dup_falg, string &sv_format);
void printStatParas(string &user_file, string &benchmark_file);

#endif /* SRC_PARAS_H_ */
