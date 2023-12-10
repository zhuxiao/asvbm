#ifndef SRC_GLOBAL_H_
#define SRC_GLOBAL_H_

#include <iostream>
#include <cstring>
#include <vector>
#include <fstream>

using namespace std;

// global variables
string outputPathname = "output";
string outputInsideToolDirname = "";
string outputBasicMetricschart = "figures";
string outputDiffRangeBasicMetricschart = "different_ranges";
string SVsizeratiofile = "sv_size_ratio";
string refRegSizeStatDirname = "1_ref_reg_size_stat/";
string numStatDirname = "2_num_stat/";
string sizeDifStatDirname = "3_size_dif_stat/";
string sizeNumStatDirname = "4_size_num_stat/";

string program_cmd_str = "";
string convertScreenFilename = "convert_screen";
string statScreenFilename = "stat_screen";
ofstream outConvertScreenFile;
ofstream outStatScreenFile;

int32_t num_threads;
int32_t extendSize;		// overlap extend size

int32_t minSizeLargeSV;
int32_t extendSizeLargeSV;

string typeMatchLevel;

float svlenRatio;
vector<float> data;
vector<int> data1;
vector< vector<float> > MeticsValues;
vector< vector<int> > MeticsValues1;
vector<float> data_4;
vector<int> data1_4;
vector< vector<float> > MeticsValues_4;
vector< vector<int> > MeticsValues1_4;

vector< vector<float> > MeticsValues4_0;
vector< vector<float> > MeticsValues4_1;
vector< vector<float> > MeticsValues4_2;
vector< vector<float> > MeticsValues4_3;
vector< vector<float> > MeticsValues4_4;
vector< vector<float> > MeticsValues4_5;
vector< vector<float> > MeticsValues4_6;

vector<vector<int>> SizeRatioV = {{},{},{},{},{},{},{},{},{}};
vector<vector<double>> SizeRatio_V = {{},{},{},{},{},{},{},{},{}};

vector<string> SVSizeRatioFile;

int32_t maxValidRegThres;
string mateItemFilename = "duplicated_mate_items.bed";
string snvFilename = "snv_items.bed";
string longSVFilename = "long_sv_reg.bed";
string svSizeDifRatioFilename = "sv_size_dif_ratio";
string svSizeDifStatFilename = "sv_size_dif_stat";
string svSizeRatioStatFilename = "sv_size_ratio_stat";

vector<size_t> size_div_vec = {100, 250, 500, 1000, 5000, 10000};
vector<double> ratio_div_vec = {0.5, 0.7, 1.2, 2, 5, 10, 50, 100};

pthread_mutex_t mtx_overlap = PTHREAD_MUTEX_INITIALIZER;


#endif /* SRC_GLOBAL_H_ */
