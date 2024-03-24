#ifndef SRC_PARAS_H_
#define SRC_PARAS_H_

#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>
#include <getopt.h>
#include <unistd.h>

#include "structure.h"
#include "convert.h"
#include "ref_reg_size_stat.h"
#include "num_stat.h"
#include "size_dif_stat.h"
#include "size_num_stat.h"
#include "util.h"
#include "meminfo.h"
#include "gnuplotcall.h"

using namespace std;

// program variables
#define PROG_NAME		"SV_STAT"
#define PROG_DESC		"A tool for Structural Variant Statistics Benchmarking"
#define PROG_VERSION	"1.0.1"



void show_version();
void showUsage();
void showUsageConvert();
void showUsageStat();
int parseConvert(int argc, char **argv);
int parseStat(int argc, char **argv);

string getProgramVersion(const string &cmd_str);
void convert(string &infilename, string &outfilename, string &reffilename, string &mate_filename, string &snv_filename, string &sv_format);
//void SVStat(string &ref_file, string &user_file, string &benchmark_file);
void SVStat(string &ref_file, string &user_file, string &benchmark_file, vector<string> &sv_files1, vector<string> &tool_names);
void SVStatOp(string &ref_file, string &sv_file1, string &sv_file2, vector<string> &sv_files1, vector<string> tool_names);

void printConvertParas(string &infilename, string &outfilename, string &reffilename, string &mate_filename, string &snv_filename, string &sv_format);
void printStatParas(string &user_file, string &benchmark_file, string &ref_file);

#endif /* SRC_PARAS_H_ */
