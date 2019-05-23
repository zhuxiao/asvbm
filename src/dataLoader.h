#ifndef SRC_DATALOADER_H_
#define SRC_DATALOADER_H_

#include <iostream>
#include <cstring>
#include <fstream>
#include <vector>

#include "structure.h"

using namespace std;

vector<SV_item*> loadData(const string &filename);
vector<SV_item*> loadDataWithoutTra(const string &filename);
vector<SV_item*> loadDataTra(const string &filename);
SV_item *constructSVItem(string &line);

void destroyData(vector<SV_item*> &sv_vec);

vector<SV_item*> getLongSVReg(vector<SV_item*> &dataset, int32_t thres);
void output2File(const string &filename, vector<SV_item*> &data);

#endif /* SRC_DATALOADER_H_ */
