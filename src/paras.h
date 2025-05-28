#ifndef SRC_PARAS_H_
#define SRC_PARAS_H_

#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>
#include <getopt.h>
#include <unistd.h>
#include <random>
#include <cstdlib>
#include <htslib/sam.h>

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
#define PROG_NAME		"ASVBM"
#define PROG_DESC		"A tool for multi-Allele-aware Structural Variants Statistics Benchmarking for Multiple callsets"
#define PROG_VERSION	"1.3.2"

#define CMD_STAT_STR					"stat"
#define CMD_CREATE_STR					"create"


void show_version();
void show();
void showUsage();
void showUsageConvert();
void showUsageCreate();
void showUsageStat();
int parseConvert(int argc, char **argv);
int parseStat(int argc, char **argv);
int parseCreate(int argc, char **argv);

bool check_bam_average_read_length();
vector<string> generateSVregion();
string getProgramVersion(const string &cmd_str);
void convert(string &infilename, string &outfilename, string &reffilename, string &mate_filename, string &snv_filename, string &sv_format);
//void SVStat(string &ref_file, string &user_file, string &benchmark_file);
void SVStat(string &ref_file, string &user_file, string &benchmark_file, vector<string> &sv_files1, vector<string> &tool_names);
void SVStatOp(string &ref_file, string &sv_file1, string &sv_file2, vector<string> &sv_files1, vector<string> tool_names);

void printConvertParas(string &infilename, string &outfilename, string &reffilename, string &mate_filename, string &snv_filename, string &sv_format);
void printStatParas(string &user_file, string &benchmark_file, string &ref_file);

#endif /* SRC_PARAS_H_ */
