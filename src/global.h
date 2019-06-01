#ifndef SRC_GLOBAL_H_
#define SRC_GLOBAL_H_

#include <iostream>
#include <cstring>
#include <vector>
#include <fstream>

using namespace std;

// global variables
string refRegSizeStatDirname = "1_ref_reg_size_stat/";
string numStatDirname = "2_num_stat/";
string sizeDifStatDirname = "3_size_dif_stat/";
string typeNumStatDirname = "4_type_num_stat/";

string statScreenFilename = "stat_screen";
ofstream outStatScreenFile;

int32_t extendSize;

int32_t maxValidRegThres;
string outputPathname = "output";
string longSVFilename = "long_sv_reg.bed";
string svSizeDifRatioFilename = "sv_size_dif_ratio";
string svSizeDifStatFilename = "sv_size_dif_stat";
string svSizeRatioStatFilename = "sv_size_ratio_stat";

vector<size_t> size_div_vec = {100, 250, 500, 1000, 5000, 10000};
vector<double> ratio_div_vec = {0.5, 0.7, 1.2, 2, 5, 10, 50, 100};


#endif /* SRC_GLOBAL_H_ */
