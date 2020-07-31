#ifndef SRC_EXTVAB_H_
#define SRC_EXTVAB_H_

#include <iostream>
#include <cstring>
#include <vector>

using namespace std;

// global variables
extern string refRegSizeStatDirname;
extern string numStatDirname;
extern string sizeDifStatDirname;
extern string sizeNumStatDirname;

extern string convertScreenFilename;
extern string statScreenFilename;
extern ofstream outConvertScreenFile;
extern ofstream outStatScreenFile;

extern int32_t extendSize;

extern int32_t maxValidRegThres;
extern string outputPathname;
extern string redundantItemFilename;
extern string longSVFilename;
extern string svSizeDifRatioFilename;
extern string svSizeDifStatFilename;
extern string svSizeRatioStatFilename;

extern vector<size_t> size_div_vec;
extern vector<double> ratio_div_vec;


#endif /* SRC_EXTVAB_H_ */
