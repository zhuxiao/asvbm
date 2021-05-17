#ifndef SRC_CONVERT_H_
#define SRC_CONVERT_H_

#include <iostream>
#include <cstring>
#include <vector>
#include <fstream>

#include "constants.h"
#include "structure.h"
#include "dataLoader.h"
#include "extvab.h"

using namespace std;


void convertBed(const string &infilename, const string &outfilename, string &redundant_filename, string &snv_filename);
void convertVcf(const string &infilename, const string &outfilename, string &redundant_filename, string &snv_filename);
void convertCsv(const string &infilename, const string &outfilename, string &redundant_filename, string &snv_filename);
void convertNm(const string &infilename, const string &outfilename, string &redundant_filename, string &snv_filename);
vector<string> getSVType(vector<string> &str_vec);
vector<int32_t> getSVLen(vector<string> &str_vec, string &sv_type);
bool isSeq(string &seq);
bool isComma(string &seq);

SV_item *allocateSVItem(string &chrname, size_t startPos, size_t endPos, string &chrname2, size_t startPos2, size_t endPos2, string &sv_type_str, int32_t sv_len);
void removeRedundantSVItems(string &redundant_filename, vector<SV_item*> &sv_item_vec);
void removeSNVItems(string &snv_filename, vector<SV_item*> &sv_item_vec);
void releaseSVItemVec(vector<SV_item*> &sv_item_vec);


#endif /* SRC_CONVERT_H_ */
