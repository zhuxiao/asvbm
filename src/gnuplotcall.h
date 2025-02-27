#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdio>
#include "extvab.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <map>
#include <unordered_map>
#include <algorithm>
#include "htmlGenerate.h"
#include "convert.h"
#include "dataLoader.h"
#include "shared_fp.h"
#include "benchmark_integration.h"
using namespace std;

void replaceUnderscoreWithDot(string& str);
string getContentAfterSlash(const string& inputString);
void ResultPresentation(vector<string> &sv_files1, string &outputPathname, vector<string> &tool_names, string &outputBasicMetricschart, vector< vector<float> > MeticsValues, vector< vector<int> > MeticsValues1, string &ref_file);
void multipledataset(vector< vector<float> > MeticsValues, vector<string> &sv_files1, vector<string> &tool_names, string &outputBasicMetricschart);
void multipledataset(vector< vector<int> > MeticsValues, vector<string> &sv_files1, vector<string> &tool_names, string &outputBasicMetricschart);
void AddfileInformation(string &FileNamePath, string &Info);
void Histogram_drawing(vector< vector<float> > MeticsValues, string &outputPathname, string &outputBasicMetricschart);
void Histogram_drawing(vector< vector<int> > MeticsValues1, string &outputPathname, string &outputBasicMetricschart);
void SvNumberDistributionGraph(int max_valid_reg_thres, string &refRegSizeFinename, string &refRegSizeFinename_tmp);
void CenterdistanceAndAreasizeratio(string &sizeDifStatDirname);
void SVsizeAndNumstatistics(string &sizeNumStatDirname, vector< vector<float> > MeticsValues_4);
void SVsizeAndNumstatistics(string &sizeNumStatDirname, vector< vector<int> > MeticsValues1_4);
void ComparisonofMetricsindifferentrangesOp(string &outputDiffRangeBasicMetricschart, vector<string> &tool_names, vector<string> &sv_files1);
void ComparisonofMetricsindifferentranges(string &DiffRangeStatDirname, vector< vector<float> > meticsvalues, vector<string> &tool_names, string &FilenName, string &outputfilename, vector<string> &sv_files1);

void GenerateFileOp(string &DiffRangeStatDirname, vector<string> &tool_names, vector<string> &filenames, vector<string> &sv_files1);
void GenerateFile(string &DiffRangeStatDirname, vector<string> &tool_names, string &FileName, vector< vector<float> > meticsvalues, vector<string> &filenames, vector<string> &sv_files1);


void GenerateMultiBarCharts(string &outputBasicMetricschart, vector<string>& fileNames, vector<string> &tool_names, vector<string> &sv_files1);

void GenerateSVsizeRatioFileGraph(vector<string> svfilesV, vector<string> toolnameV, string &FileSavingPath);

void RscriptGeneration(string &filePath);
void TPfilemanipulation(string &filePath);
void RunRscriptoperation(string &RscriptPath, string &TPfilesPath);
void GenerateSharedFNfile(string& outputFile);
void findCommonFN(string& outputFile, vector<string>& inputFiles);

void SvNumberDistributionGraph(const vector<string> &fileNamesPaths, string& svdistributionDirPath);
//extern void findSharedFP(std::string& outputFile, std::vector<std::string>& sharedFP);
bool numericSort(const pair<string, int>& lhs, const pair<string, int>& rhs);
