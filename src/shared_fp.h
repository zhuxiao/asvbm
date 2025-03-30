/*
 * find_shared_fp.h
 *
 *  Created on: Oct 23, 2024
 *      Author: node
 */

#ifndef SRC_SHARED_FP_H_
#define SRC_SHARED_FP_H_
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <set>
#include <unordered_set>
#include <tuple>
#include <thread>
#include <future>
#include "extvab.h"
#include "gnuplotcall.h"
#include "num_stat.h"
#include "vcf_header.h"




void addEntryMap(size_t index, SV_item* entry);
vector<string> split(const string& str, char delimiter);
string appendCustomInfo(const string& originalInfo, const string& customInfo);
int calculateOverlap(SV_item* item1, SV_item* item2);
void releaseSVItems(vector<SV_item*>& items);
bool hasMatchingPairs(const unordered_map<int, set<pair<int, int>>>& map, int index1, int index2);
void match(vector<vector<SV_item*>>& sets, size_t l, vector<vector<int>>& subsetStartIndex, faidx_t *fai, size_t totalEntries, vector<int>& entriesStartIndex, vector<vector<int>>& connectedComponents);
void findSharedFP(string &outputFile, string &ref_file, string &benchfilename, vector<string> &sv_files1);
#endif /* SRC_SHARED_FP_H_ */
