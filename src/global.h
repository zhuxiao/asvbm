#ifndef SRC_GLOBAL_H_
#define SRC_GLOBAL_H_

#include <iostream>
#include <cstring>
#include <vector>
#include <fstream>
#include <unordered_map>

using namespace std;

// global variables
string outputPathname = "output";
string outputInsideToolDirname = "";
string suboutputDirname = "convert";
string suboutputDirnamePath;
string outputBasicMetricschart = "figures";
string outputCommonFN = "SharedFN_RefinedBenchmark";
string outputBasicMetricschartPath;
string outputDiffRangeBasicMetricschart = "different_ranges";
string outputUpSetRchart = "UpSetR";
string SVsizeratiofile = "sv_size_ratio";
string SVdistributionDirname = "sv_distribution";
string SVdistributionPath;
string refRegSizeStatDirname = "1_ref_reg_size_stat/";
string numStatDirname = "2_num_stat/";
string sizeDifStatDirname = "3_size_dif_stat/";
string sizeNumStatDirname = "4_size_num_stat/";

string acquiesce_count = "method";
string method_name = "method";
string program_cmd_str = "";
string convertScreenFilename = "convert_screen";
string statScreenFilename = "stat_screen";
string Rscriptfilename = "UpSetR.R";
string RscriptfilePNGname = "upset_plot.png";
string RscriptfilePNGnamePath;
string RscriptfilenamePath;
ofstream outConvertScreenFile;
ofstream outStatScreenFile;

int32_t num_threads;
int32_t extendSize;		// overlap extend size

int32_t minSizeLargeSV;
int32_t extendSizeLargeSV;

vector<string> chromosomeSet;

string typeMatchLevel;

vector<string> folderPng1;
vector<string> folderPng3;
vector<string> folderPng4;
vector<string> folderPng5;
vector<string> alltoolnames;
vector<string> allmetric;
vector< vector<string>> allmetrics = {{"Tool", "#SVs_bench", "#SVs_user", "#SVs_filtered_user", "#TP_bench", "#TP_user", "#FP", "#FN", "Identity", "Recall", "Precision", "F1 score"}};
vector<string> centerDistance;
vector< vector<string>> centerDistances = {{"Tool", "-200--151", "-150--101", "-100--51", "-50--1", "0-50", "51-100", "101-150", "151-200"}};
vector<string> sizeratio;
vector< vector<string>> sizeratios = {{"Tool", "0.0-0.5", "0.5-0.7", "0.7-1.2", "1.2-2.0", "2.0-5.0", "5.0-10.0", "10.0-50.0", "50.0-100.0", ">100.0"}};
vector<string> region = {"1-100bp", "101-250bp", "251-500bp", "501-1000bp", "1001-2500bp", "2501-5000bp", "5001-10000bp", ">10000bp"};
vector<string> regionmetric;
vector<vector<vector<string>>> allregionmetrics;
vector< vector<string>> regionmetrics = {{"Region", "#TP_bench", "#TP_user", "#FP", "#FN", "Identity", "Recall", "Precision", "F1 score"}};
vector<string> regSizeFiles;
unordered_map<std::string, int> benchmarklineMap;
vector<string> benchmarkannotationLines;

double Seqcons;
double SeqconsSum;
int32_t SeqconsNum;
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
vector< vector<float> > MeticsValues4_7;

vector<vector<int>> SizeRatioV = {{},{},{},{},{},{},{},{},{}};
vector<vector<double>> SizeRatio_V = {{},{},{},{},{},{},{},{},{}};

vector<string> SVSizeRatioFile;
vector<string> SVcallernames;
vector<string> TPfilesPath;
vector<string> FNfilesPath;

int32_t maxValidRegThres;
string mateItemFilename = "duplicated_mate_items.bed";
string snvFilename = "snv_items.bed";
string longSVFilename = "long_sv_reg.bed";
string svSizeDifRatioFilename = "sv_size_dif_ratio";
string svSizeDifStatFilename = "sv_size_dif_stat";
string svSizeRatioStatFilename = "sv_size_ratio_stat";
string htmlFilename = "asvbm_reports.html";

vector<size_t> size_div_vec = {100, 250, 500, 1000, 2500, 5000, 10000};
vector<double> ratio_div_vec = {0.5, 0.7, 1.2, 2, 5, 10, 50, 100};

pthread_mutex_t mtx_overlap = PTHREAD_MUTEX_INITIALIZER;

#endif /* SRC_GLOBAL_H_ */
