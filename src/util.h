#ifndef SRC_UTIL_H_
#define SRC_UTIL_H_

#include <iostream>
#include <cstring>
#include <vector>

#include "structure.h"

using namespace std;

vector<string> split(const string& s, const string& delim);
string double2Str(double num, int32_t fixedPrecision);
bool isBase(const char ch);
void printSV(vector<SV_item*> &data);
string getProgramCmdStr(int argc, char *argv[]);

#endif /* SRC_UTIL_H_ */
