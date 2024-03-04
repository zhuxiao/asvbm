#ifndef SRC_EXTVAB_H_
#define SRC_EXTVAB_H_

#include <iostream>
#include <cstring>
#include <vector>

using namespace std;

// global variables

extern string outputPathname;
extern string outputInsideToolDirname;
extern string suboutputDirname;
extern string suboutputDirnamePath;
extern string outputBasicMetricschart;
extern string outputBasicMetricschartPath;
extern string outputDiffRangeBasicMetricschart;
extern string SVsizeratiofile;
extern string refRegSizeStatDirname;
extern string numStatDirname;
extern string sizeDifStatDirname;
extern string sizeNumStatDirname;

extern string acquiesce_count;
extern string method_name;
extern string program_cmd_str;
extern string convertScreenFilename;
extern string statScreenFilename;
extern ofstream outConvertScreenFile;
extern ofstream outStatScreenFile;

extern int32_t num_threads;
extern int32_t extendSize;

extern int32_t minSizeLargeSV;
extern int32_t extendSizeLargeSV;

extern vector<string> folderPng1;
extern vector<string> folderPng3;
extern vector<string> folderPng4;
extern vector<string> folderPng5;
extern vector<string> alltoolnames;
extern vector<string> allmetric;
extern vector< vector<string>> allmetrics;
extern vector<string> centerDistance;
extern vector< vector<string>> centerDistances;
extern vector<string> sizeratio;
extern vector< vector<string>> sizeratios;
extern vector<string> region;
extern vector<string> regionmetric;
extern vector<vector<vector<string>>> allregionmetrics;
extern vector< vector<string>> regionmetrics;

extern double Seqcons;
extern double SeqconsSum;
extern int32_t SeqconsNum;
extern float svlenRatio;
extern vector<string> chromosomeSet;
extern string typeMatchLevel;
extern vector<float> data;
extern vector<int> data1;
extern vector< vector<float> > MeticsValues;
extern vector< vector<int> > MeticsValues1;
extern vector<float> data_4;
extern vector<int> data1_4;
extern vector< vector<float> > MeticsValues_4;
extern vector< vector<int> > MeticsValues1_4;

extern vector< vector<float> > MeticsValues4_0;
extern vector< vector<float> > MeticsValues4_1;
extern vector< vector<float> > MeticsValues4_2;
extern vector< vector<float> > MeticsValues4_3;
extern vector< vector<float> > MeticsValues4_4;
extern vector< vector<float> > MeticsValues4_5;
extern vector< vector<float> > MeticsValues4_6;

extern vector<vector<int>> SizeRatioV;
extern vector<vector<double>> SizeRatio_V;

extern vector<string> SVSizeRatioFile;

extern int32_t maxValidRegThres;
extern string snvFilename;
extern string mateItemFilename;
extern string longSVFilename;
extern string svSizeDifRatioFilename;
extern string svSizeDifStatFilename;
extern string svSizeRatioStatFilename;
extern string htmlFilename;

extern vector<size_t> size_div_vec;
extern vector<double> ratio_div_vec;

extern pthread_mutex_t mtx_overlap;

#endif /* SRC_EXTVAB_H_ */
