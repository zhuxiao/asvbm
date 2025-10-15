#ifndef SRC_CONVERT_H_
#define SRC_CONVERT_H_

#include <iostream>
#include <cstring>
#include <vector>
#include <fstream>
#include <algorithm>

#include "constants.h"
#include "structure.h"
#include "dataLoader.h"
#include "extvab.h"

using namespace std;

#define Max_SeqLen    50000  //The maximum length of the sequence


string Pathquerybackslash(string filename);
string PathqueryDot(string filename);
bool isExistChromosomeSet(string &chrname);
bool isDecoyChr(string &chrname);
void processFile(const string &infilename, bool &isVCF, bool &isBED);
void convertBed(const string &infilename, const string &outfilename, const string &reffilename, string &mate_filename, string &snv_filename, string &label);
void convertVcf(const string &infilename, const string &outfilename, const string &reffilename, string &mate_filename, string &snv_filename, string &label);
void convertCsv(const string &infilename, const string &outfilename, const string &reffilename, string &mate_filename, string &snv_filename, string &label);
void convertNm(const string &infilename, const string &outfilename, const string &reffilename, string &mate_filename, string &snv_filename, string &label);
vector<string> getSVType(vector<string> &str_vec);
vector<int32_t> getSVLen(vector<string> &str_vec, string &sv_type);
vector<int32_t> getSVLen(vector<string> &str_vec, vector<string> &sv_type_vec);
bool isSeq(string &seq);
bool isComma(string &seq);

SV_item *allocateSVItem(string &chrname, size_t startPos, size_t endPos, string &chrname2, size_t startPos2, size_t endPos2, string &sv_type_str, int32_t sv_len, string &ref_seq, string &alt_seq, string &line, string &tab);
void removemateSVItems(string &mate_filename, vector<SV_item*> &sv_item_vec);
void removeSNVItems(string &snv_filename, vector<SV_item*> &sv_item_vec);
void releaseSVItemVec(vector<SV_item*> &sv_item_vec);


#endif /* SRC_CONVERT_H_ */
