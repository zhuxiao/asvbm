/*
 * find_shared_fp.cpp
 *
 *  Created on: Oct 23, 2024
 *      Author: node
 */
#include "shared_fp.h"
int lcs(const string& a, const string& b){
	int m, n, i, j;
	m = a.length();
	n = b.length();
	vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));
	for(i=1; i<=m; i++){
		for(j=1; j<=n; j++){
			if(a[i - 1] == b[j - 1]){
				dp[i][j] = dp[i - 1][j - 1] + 1;
			}else{
				dp[i][j] = max(dp[i - 1][j], dp[i][j - 1]);
			}
		}
	}
	return dp[m][n];
}
double lcsRatio(const string& a, const string& b){
	int lcs_length, min_length;
	lcs_length = lcs(a, b);
	min_length = min(a.length(), b.length());
	return static_cast<double>(lcs_length) / min_length;
}
//void sort(vector<SV_item*>& shared_fp_vec){
//	map<string, int> chrname_to_index;
//	int i;
//	for(i=1; i<=22; i++){
//		chrname_to_index[to_string(i)] = i;
//	}
//	chrname_to_index["X"] = 23;
//	chrname_to_index["Y"] = 24;
//	sort(shared_fp_vec.begin(), shared_fp_vec.end(),[&chrname_to_index](const SV_item* a, const SV_item* b){
//		int indexA = chrname_to_index.at(a->chrname);
//		int indexB = chrname_to_index.at(b->chrname);
//		if(indexA != indexB){
//			return indexA < indexB;
//		}
//		return a->startPos < b->startPos;
//	});
//}
vector<string> split(const string& str, char delimiter){
	vector<string> tokens;
	string token;
	for(char ch : str){
		if(ch != delimiter)
			token += ch;
		else{
			if(!token.empty())
				tokens.push_back(token);
			token.clear();
		}
	}
	if(!token.empty())
		tokens.push_back(token);
	return tokens;
}
string appendCustomInfo(const string& originalInfo, const string& customInfo){
	if(originalInfo.empty()){
		return customInfo;
	}else{
		return originalInfo + ";" + customInfo;
	}
}
void releaseSVItems(vector<SV_item*>& items){
	for(SV_item* item : items){
		delete item;
	}
	items.clear();
}
int calculateOverlap(SV_item* item1, SV_item* item2){
	if(item1->sv_type==VAR_INS){
		item1->endPos=item1->startPos+item1->sv_len;
		item2->endPos=item2->startPos+item2->sv_len;
	}
	return (max(item1->startPos, item2->startPos) < min(item1->endPos, item2->endPos)) ? (min(item1->endPos, item2->endPos) - max(item1->startPos, item2->startPos)) : 0;
}
bool hasMatchingPairs(const unordered_map<int, set<pair<int, int>>>& map, int index1, int index2){
	if(map.find(index1) != map.end() && map.find(index2) != map.end() && !map.at(index1).empty() && !map.at(index2).empty()){
		for(const auto& pair1 : map.at(index1)){
			for(const auto& pair2 : map.at(index2)){
				if((pair1.first == pair2.first || pair1.first == pair2.second || pair1.second == pair2.first || pair1.second == pair2.second)){
					return true;
				}
			}
		}
	}
	return false;
}
//entryMap:<index,entry>
void addEntryMap(size_t index, SV_item *entry){
	if(index >= 0 && entry != nullptr){
		entryMap.push_back(entry);
	}else{
		cerr<< "Index out of bounds or null pointer in addEntryMap function: "<< index << endl;
	}
}
void filter(vector<vector<SV_item*>>& sets, size_t l, vector<int>& subset2StartIndex, vector<SV_item*>& shared_fp_vec){

}
void match(vector<vector<SV_item*>>& sets, size_t l, vector<vector<int>>& subsetStartIndex, faidx_t *fai, size_t totalEntries, vector<int>& entriesStartIndex, vector<vector<int>>& connectedComponents){
	size_t i, j, m, n;
	vector<size_t> overlap_type_vec;
	vector<SV_item*> subset1, subset2;
	int64_t start_idx, end_idx, new_start_idx, start_search_pos, end_search_pos;
	int64_t start_idx_tmp0, end_idx_tmp0, sv_len1, sv_len2;
	float svlenRatio_tmp1;
	bool flag, foundIndex;
	double consistency;
	size_t index1, index2;
	pair<int, int> connectedPair;
	unordered_map<int, set<pair<int, int>>> map;
	for(m=0; m<sharedFPFilenames.size(); m++){
		for(n=m+1; n<sharedFPFilenames.size(); n++){
			subset1=sets[m];
			subset2=sets[n];
			SV_item *item1, *item2;
			start_idx = end_idx = -1;
			for(i=0; i<subset1.size(); i++){
				item1 = subset1.at(i);
				start_idx_tmp0 = start_idx;
				end_idx_tmp0 = end_idx;
				start_search_pos = item1->startPos - tolerance;
				end_search_pos = item1->endPos + tolerance;
				if(start_search_pos < 1)
					start_search_pos = 1;
				if(start_idx_tmp0 == -1)
					start_idx = 0;
				else{
					new_start_idx = -1;
					for(j=start_idx_tmp0; j<subset2.size(); j++){
						item2 = subset2.at(j);
						overlap_type_vec = computeOverlapType(item1, item2);
						if(overlap_type_vec.size() > 0 && overlap_type_vec.at(0) != NO_OVERLAP){
							new_start_idx = j;
							break;
						}else if(item2->endPos >= start_search_pos){
							new_start_idx = j;
							break;
						}
					}
					if(new_start_idx != -1)
						start_idx = new_start_idx;
					else
						start_idx = end_idx_tmp0;
				}
				if(start_idx != -1){
					end_idx = -1;
					for(j=start_idx; j<subset2.size(); j++){
						item2 = subset2.at(j);
						flag = false;
						if(item2->startPos < start_search_pos)
							continue;
						foundIndex = false;
						index1 = entriesStartIndex[m] + subsetStartIndex[l][m] + i;
						index2 = entriesStartIndex[n] + subsetStartIndex[l][n] + j;
						if(hasMatchingPairs(map, index1, index2)){
							foundIndex = true;
							flag=true;
						}

						if(!foundIndex){
							overlap_type_vec = computeOverlapType(item1, item2);
							if(overlap_type_vec.size() > 0 && overlap_type_vec.at(0) != NO_OVERLAP){
								sv_len1 = item1->sv_len;
								sv_len2 = item2->sv_len;
								if(sv_len1 < 0)
									sv_len1 = -sv_len1;
								if(sv_len2 < 0)
									sv_len2 = -sv_len2;
								if(sv_len1 <= sv_len2)
									svlenRatio_tmp1 = (float) sv_len1 / sv_len2;
								else
									svlenRatio_tmp1 = (float) sv_len2 / sv_len1;
								if(svlenRatio_tmp1 >= 0.7){
									if(item1->sv_type == item2->sv_type or (item1->sv_type == VAR_INS and item2->sv_type== VAR_DUP) or (item1->sv_type == VAR_DUP and item2->sv_type== VAR_INS)){
	//    									if(i!=subset1.size()-1 and j!=subset2.size()-1){
	//    										if((subset1.at(i+1)->startPos - item1->startPos <= ALLELIC_DISTANCE and subset1.at(i+1)->sv_type == item1->sv_type) or (i>0 and item1->startPos - subset1.at(i-1)->startPos <= ALLELIC_DISTANCE and subset1.at(i-1)->sv_type == item1->sv_type))	alleleflag = true;
	//    										if((subset2.at(j+1)->startPos - item2->startPos <= ALLELIC_DISTANCE and subset2.at(j+1)->sv_type == item2->sv_type) or (j>0 and item2->startPos - subset2.at(j-1)->startPos <= ALLELIC_DISTANCE and subset2.at(j-1)->sv_type == item2->sv_type))	alleleflag = true; //or item2->startPos - subset2.at(j-1)->startPos <=10
	//    									}else{
	//    										if(i>0 and item1->startPos - subset1.at(i-1)->startPos <= ALLELIC_DISTANCE and subset1.at(i-1)->sv_type == item1->sv_type)	alleleflag = true;
	//    										if(j>0 and item2->startPos - subset2.at(j-1)->startPos <= ALLELIC_DISTANCE and subset2.at(j-1)->sv_type == item2->sv_type)	alleleflag = true;
	//    									}
	//    									if(alleleflag){
	//    										}
										if(item1->alt_seq.compare("-") == 0 or item2->alt_seq.compare("-") == 0){
											flag = true;
										}else{
											consistency = computeVarseqConsistency(item1, item2, fai, percentAlleleSeqSim);
											if(consistency >= percentAlleleSeqSim){
												//SeqConsNumStat(consistency);
												flag = true;
											}
										}
									}
								}
							}
						}
						if(flag){
							item1->overlapped = true;
							item2->overlapped = true;
							end_idx = j;
//							index1 = entriesStartIndex[m] + subsetStartIndex[l][m] + i;
//							index2 = entriesStartIndex[n] + subsetStartIndex[l][n] + j;
							if(index1 < totalEntries && index2 < totalEntries){
								if(entryMap[index1] != nullptr && entryMap[index2] != nullptr){
									if(index1 < entryMap.size() && index2 < entryMap.size()){
//											connectedComponents[index1].push_back(entryMap[index2]);
//											connectedComponents[index2].push_back(entryMap[index1]);
										connectedComponents[index1].push_back(index2);
										connectedComponents[index2].push_back(index1);
										if(map.find(index1) == map.end()){
											map[index1] = set<pair<int, int>>();
										}
										if(map.find(index2) == map.end()){
											map[index2] = set<pair<int, int>>();
										}
										connectedPair = make_pair(index1,index2);
										map[index1].insert(connectedPair);
										map[index2].insert(connectedPair);
									}else{
										cerr<< "Index out of bounds: index1 "<< index1 << " index2 "<< index2 << endl;
									}
								}else{
									cerr << "Null pointer " << i<< endl;
								}
								// unite(index1, index2);
							}else{
								cerr << "Index out of bounds: index1 "<< index1 << " index2 " << index2<< " total size: " << totalEntries<< endl;
							}
						}
						if(item2->startPos > end_search_pos){
							end_idx = j - 1;
							if(end_idx < 0)
								end_idx = 0;
							if(end_idx < start_idx)
								end_idx = start_idx;
							break;
						}
//						start_idx_tmp = start_idx;
//						end_idx_tmp = end_idx;
					}
				}
			}
		}
	}
}

void findSharedFP(string &outputFile, string &ref_file, string &benchfilename, vector<string> &sv_files1){
	vector<string> headers = mergeHeaders(sv_files1);
	faidx_t *fai;
	fai = fai_load(ref_file.c_str());
	string newFilename = outputFile + "shared_fp.vcf";
	string newFilename1 = outputFile + "IntegrateBenchmark.vcf";
//	string newFilename2 = outputFile + "refinedshared_fp.vcf";
	//string newFilename3 = outputFile + "benchmarkAddTP_user.bed";
	vector<int> entriesStartIndex(sharedFPFilenames.size(), 0);
	//unordered_map<int, vector<SV_item*>> connectedComponents;
	//vector<vector<SV_item*>> connectedComponents;
	vector<vector<int>> connectedComponents;
	unordered_map<string, int> benchmarkMap_nonFNS;
	pair<string, int> sharedFPPair;
	size_t i, j, m, n, l, k, o;
	size_t totalEntries = 0;
	int entryIndex = 0;
	size_t selectToolsNum;
//	size_t input;
//	bool validInput;

//	vector<SV_item*> sv_item_bench_vec, sv_item_vec1;
//	vector<SV_item*> tp_vec;
//	pair<SV_item*, int> tp_vecPair;
//	unordered_map<SV_item*, int> tp_vecMap;
//	sv_item_bench_vec = loadData(benchfilename);
//	for(i=0; i<benchmarkAddTP_user.size(); i++){
//		sv_item_vec1 = loadData(benchmarkAddTP_user[i]);
//		for(SV_item* item : sv_item_vec1){
//			tp_vec.push_back(item);
//		}
//	}
//	sort(tp_vec);

	selectToolsNum = ceil(sharedFPFilenames.size() * 0.6);
	cout << "Do you want to set the number of selected tools? (1 for Yes, 0 for No): ";
//	while(true){
//		cin >> input;
//		if(cin.fail() || input < 0 || input > 1){
//			cin.clear();
//			cin.ignore(numeric_limits<streamsize>::max(), '\n');
//			cout << "Invalid choice. Please enter 1 for Yes or 0 for No." << endl;
//		}else{
//			break;
//		}
//	}
//	if(input == 1){
//		validInput = false;
//		while(!validInput){
//			cout << "Set your number (between 2 and " << sharedFPFilenames.size() << "): ";
//			cin >> input;
//			if(cin.fail() || input < 2 || input > sharedFPFilenames.size()){
//				cin.clear();
//				cin.ignore(numeric_limits<streamsize>::max(), '\n');
//				cout << "Invalid input. Please try again." << endl;
//			}else{
//				validInput = true;
//			}
//		}
//	}else if(input == 0){
//	}
	cout << "Selected value: " << selectToolsNum << endl;
	for(i=0; i<sharedFPFilenames.size(); i++){
		if(i > 0){
			entriesStartIndex[i] = entriesStartIndex[i - 1] + fp_vec[i - 1].size();
		}else{
			entriesStartIndex[i] = 0;
		}
		totalEntries += fp_vec[i].size();
		for(j=0; j<fp_vec[i].size(); j++){
			addEntryMap(entryIndex, fp_vec[i][j]);
//			connectedComponents[entryIndex].push_back(entryMap[entryIndex]);
			entryIndex++;
		}
	}
//	connectedComponents.reserve(totalEntries);
//	connectedComponents.resize(totalEntries, vector<SV_item*>());
//	for(j=0; j<totalEntries; j++){
//		connectedComponents[j].push_back(entryMap[j]);
//	}
	connectedComponents.clear();
	connectedComponents.resize(totalEntries, vector<int>());
	for(j=0; j<totalEntries; j++){
		connectedComponents[j].push_back(j);
	}
	set<string> chrname_set1, chrname_set2;
	set<string> chrname_set_tmp;
	vector<string> chrname_vec_sorted;
	vector<vector<vector<SV_item*>>> subsets;
	vector<SV_item*> subset;
	string chr;
	vector<vector<int>> subsetStartIndex;
	chrname_set1 = getChrnames(fp_vec[0]);
	chrname_set2 = getChrnames(fp_vec[1]);
	chrname_set_tmp = getChrUnion(chrname_set1, chrname_set2);
	chrname_vec_sorted = sortChrnames(chrname_set_tmp); // sort chromosome names
	subsets.resize(chrname_vec_sorted.size(), vector<vector<SV_item*>>());
	subsetStartIndex.resize(chrname_vec_sorted.size(), vector<int>(sharedFPFilenames.size(), 0));
	for(i=0; i<chrname_vec_sorted.size(); i++){
		//subsets[i].resize(sharedFPFilenames.size(), vector<SV_item*>());
		chr = chrname_vec_sorted.at(i);
		for(l=0; l<sharedFPFilenames.size(); l++){
			subset = getItemsByChr(chr, fp_vec[l]);
			subsets[i].push_back(subset);
		}
	}
	for(i=1; i<chrname_vec_sorted.size(); i++){
		for(l=0; l<sharedFPFilenames.size(); l++){
			subsetStartIndex[i][l] = subsetStartIndex[i-1][l] + subsets[i-1][l].size();
		}
	}
	vector<thread> threads;
	for(l=0; l<chrname_vec_sorted.size(); l++){
		//vector<vector<SV_item*>> sets = subsets[l];
		threads.emplace_back(match, ref(subsets[l]), l, ref(subsetStartIndex), fai, totalEntries, ref(entriesStartIndex), ref(connectedComponents));

	}
	for(vector<thread>::iterator it=threads.begin(); it!=threads.end(); it++){
		if (it->joinable()){
			it->join();
		}
	}

	//Write matched variants
	vector<size_t> connectedIndex;
//	unordered_set<size_t> connectedIndex;
	vector<int> differentToolsIndex;
	vector<int> mostRepresentativeVariantIndices;
	unordered_map<int, int> sizeDifferenceCount;
	ofstream outputFile1(newFilename);
	ofstream outputFile2(newFilename1);
//	ofstream outputFile3(newFilename2);
	//ofstream outputFile4(newFilename3);
	vector<string> fields;
	vector<SV_item*> shared_fp_vec;
	string infoField, customInfo, newInfoField, newLineInfo, tools;
	SV_item* currentVariant;
	SV_item* mostRepresentativeVariant;
	SV_item* otherVariant;
	bool flag;
	size_t size;
	int totalOverlapForCurrentVariant, maxTotalOverlap, mostRepresentativeVariantIndex, toolsNum, count, maxCount, len;
//	outputFile1 << VCF_HEADER << endl;
//	outputFile2 << VCF_HEADER << endl;
	for (const auto& header : headers){
		outputFile1 << header << endl;
		outputFile2 << header << endl;
	}
//	outputFile3 << VCF_HEADER << endl;
	for(i=0; i<connectedComponents.size(); i++){
		size = connectedComponents[i].size();
		if(size >= selectToolsNum){
			differentToolsIndex.clear();
			for(j=0; j<connectedComponents[i].size(); j++){
				for(l=0; l<sharedFPFilenames.size(); l++){
					if(l != (sharedFPFilenames.size() - 1)){
						if(connectedComponents[i][j] < entriesStartIndex[l + 1]){
							if(find(differentToolsIndex.begin(), differentToolsIndex.end(), l) == differentToolsIndex.end()){
								differentToolsIndex.push_back(l);
							}
							break;
						}
					}else{
						if(find(differentToolsIndex.begin(), differentToolsIndex.end(), l) == differentToolsIndex.end()){
							differentToolsIndex.push_back(l);
						}
						break;
					}
				}
			}

			if(differentToolsIndex.size() >= selectToolsNum){
			flag = true;

			if(!connectedIndex.empty()){
				for(m=0; m<connectedIndex.size(); m++){
					if(i == connectedIndex[m]){
						flag = false;
						break;
					}
				}
			}
			for(j=0; j<connectedComponents[i].size(); j++){
				for(m=0; m<connectedComponents[connectedComponents[i][j]].size(); m++){
					if(find(connectedIndex.begin(), connectedIndex.end(), connectedComponents[connectedComponents[i][j]][m]) == connectedIndex.end()){
						connectedIndex.push_back(connectedComponents[connectedComponents[i][j]][m]);
					}
				}
			}


				if(flag){
					if(entryMap[i] != nullptr){
						maxCount = 0;
						len = 0;
						while(maxCount == 0){
							len = len + 10;
							for(k=0; k<connectedComponents[i].size(); k++){
								differentToolsIndex.clear();
								for(j=0; j<connectedComponents[connectedComponents[i][k]].size(); j++){
									for(l=0; l<sharedFPFilenames.size(); l++){
										if(l != (sharedFPFilenames.size() - 1)){
											if(connectedComponents[connectedComponents[i][k]][j] < entriesStartIndex[l + 1]){
												if(find(differentToolsIndex.begin(), differentToolsIndex.end(), l) == differentToolsIndex.end()){
													differentToolsIndex.push_back(l);
												}
												break;
											}
										}else{
											if(find(differentToolsIndex.begin(), differentToolsIndex.end(), l) == differentToolsIndex.end()){
												differentToolsIndex.push_back(l);
											}
										}
									}
								}
								if(differentToolsIndex.size() >= selectToolsNum){
	//							if (connectedComponents[connectedComponents[i][k].second].size() >= (sharedFPFilenames.size()/2)+1){
									currentVariant = entryMap[connectedComponents[i][k]];
									count=0;
									for(o=0; o<connectedComponents[i].size(); o++){
										if(k != o){
											otherVariant = entryMap[connectedComponents[i][o]];
											if(abs(currentVariant->sv_len - otherVariant->sv_len) <= len){
												count++;
											}
										}
									}
									sizeDifferenceCount.insert(make_pair(connectedComponents[i][k], count));
								}
							}
							for(const auto& kv : sizeDifferenceCount){
								if(kv.second > maxCount){
									maxCount = kv.second;
									mostRepresentativeVariantIndices.clear();
									mostRepresentativeVariantIndices.push_back(kv.first);
								}else if(kv.second == maxCount){
									mostRepresentativeVariantIndices.push_back(kv.first);
								}
							}
							sizeDifferenceCount.clear();
							if(maxCount == 0) mostRepresentativeVariantIndices.clear();
						}

						maxTotalOverlap = -1;
						for(k=0; k<mostRepresentativeVariantIndices.size(); k++){
							currentVariant = entryMap[mostRepresentativeVariantIndices[k]];
							totalOverlapForCurrentVariant = 0;
							for(o=0; o<mostRepresentativeVariantIndices.size(); o++){
								if(k != o){
									otherVariant = entryMap[mostRepresentativeVariantIndices[o]];
									totalOverlapForCurrentVariant += calculateOverlap(currentVariant,otherVariant);
								}
							}
							if(totalOverlapForCurrentVariant >= maxTotalOverlap){
								maxTotalOverlap = totalOverlapForCurrentVariant;
								mostRepresentativeVariant = currentVariant;
								mostRepresentativeVariantIndex = mostRepresentativeVariantIndices[k];
							}
						}
						mostRepresentativeVariantIndices.clear();


						toolsNum = 0;
						tools = "";
						for(k=0; k<connectedComponents[mostRepresentativeVariantIndex].size(); k++){
							for(l=0; l<sharedFPFilenames.size(); l++){
								if(l != (sharedFPFilenames.size() - 1)){
									if(connectedComponents[mostRepresentativeVariantIndex][k] < entriesStartIndex[l + 1]){
										if(k!= 0){
											if(tools.find(SVcallernames[l]) == string::npos) { // 检查tools是否已包含SVcallernames[l]
												tools += "," + SVcallernames[l];
												toolsNum++;
												}
											break;
										}else{
											if(tools.find(SVcallernames[l]) == string::npos) { // 检查tools是否已包含SVcallernames[l]
												tools += SVcallernames[l];
												toolsNum++;
											}
											break;
										}
									}
								}else if(k != 0){
									if (tools.find(SVcallernames[l]) == string::npos) { // 检查tools是否已包含SVcallernames[l]
										tools += "," + SVcallernames[l];
										toolsNum++;
									}
								}else{
									if(tools.find(SVcallernames[l]) == string::npos) { // 检查tools是否已包含SVcallernames[l]
										tools += SVcallernames[l];
										toolsNum++;
									}
								}
							}
						}
//						shared_fp_vec.push_back(mostRepresentativeVariant);
						fields = split(mostRepresentativeVariant->lineInfo, '\t');
						//outputFile1 << entryMap[i]->lineInfo << endl;
						infoField = fields[7];
						customInfo = "TOOLSNUM=" + to_string(toolsNum) + ";TOOLS=" + tools;
						newInfoField = appendCustomInfo(infoField, customInfo);
						fields[7] = newInfoField;
						newLineInfo = "";
						for(n=0; n<fields.size(); n++){
							newLineInfo += fields[n] + (n < fields.size() - 1 ? "\t" : "");
						}
						sharedFPPair = make_pair(newLineInfo,1);
						sharedFPlineMap.insert(sharedFPPair);
//						for(j=0;j<connectedComponents[i].size();j++){
//							if(entryMap[connectedComponents[i][j]] == mostRepresentativeVariant){
//								connectedComponents[i].erase(connectedComponents[i].begin() + j);
//								break;
//							}
//						}

						//outputFile1 << newLineInfo << endl;
					}else{
						//outputFile1 << "NULL pointer" << endl;
					}
				}
			}
//			size--;
		}
	}

//	sortSharedFPVec(shared_fp_vec);
	vector<pair<string, int>> sortedVec1(sharedFPlineMap.begin(), sharedFPlineMap.end());
	sort(sortedVec1.begin(), sortedVec1.end(), numericSort);
	for(const auto& pair : sortedVec1){
			outputFile1 << pair.first << endl;
	}
	string mate_filename = "", snv_filename = "";
	string label1 = "user";
	string sharedFP = outputCommonFN + "/" + "shared_fp.bed";
	convertVcf(newFilename, sharedFP, ref_file, mate_filename, snv_filename, label1);
	shared_fp_vec = loadData(sharedFP);
//	vector<vector<vector<SV_item*>>> subsets2;
//	vector<int> subset2StartIndex;
//	subsets2.resize(chrname_vec_sorted.size(), vector<vector<SV_item*>>());
//	subset2StartIndex.resize(chrname_vec_sorted.size(), 0);
//	for(i=0; i<chrname_vec_sorted.size(); i++){
//		chr = chrname_vec_sorted.at(i);
//		subsets2[i].push_back(getItemsByChr(chr, shared_fp_vec));
//		subsets2[i].push_back(getItemsByChr(chr, sv_item_bench_vec));
//		subsets2[i].push_back(getItemsByChr(chr, tp_vec));
//	}
//	for(i=1; i<chrname_vec_sorted.size(); i++){
//		subset2StartIndex[i] = subset2StartIndex[i-1] + subsets2[i-1][0].size();
//	}
//	vector<thread> threads2;
//	for(l=0; l<chrname_vec_sorted.size(); l++){
//		threads2.emplace_back(filter, ref(subsets2[l]), l, ref(subset2StartIndex), ref(shared_fp_vec));
//	}
//	for(vector<thread>::iterator it=threads2.begin(); it!=threads2.end(); it++){
//		if (it->joinable()){
//			it->join();
//		}
//	}
//	for(i=0; i<shared_fp_vec.size(); i++){
//		if(shared_fp_vec[i] == nullptr) continue;
//		refinedsharedFPlineMap.insert(make_pair(shared_fp_vec[i]->lineInfo,1));
//	}
	for(const auto& pair : benchmarklineMap){
		if(pair.second != int(SVcallernames.size())){
			benchmarkMap_nonFNS.insert(pair);
		}
	}
	for(const auto& pair : sharedFPlineMap){
		benchmarkMap_nonFNS.insert(pair);
	}
	vector<pair<string, int>> sortedVec(benchmarkMap_nonFNS.begin(), benchmarkMap_nonFNS.end());
//	vector<pair<string, int>> sortedVec2(refinedsharedFPlineMap.begin(), refinedsharedFPlineMap.end());
	sort(sortedVec.begin(), sortedVec.end(), numericSort);
//	sort(sortedVec2.begin(), sortedVec2.end(), numericSort);
	for(const auto& pair : sortedVec){
		outputFile2 << pair.first << endl;
	}
//	for(const auto& pair : sortedVec2){
//			outputFile3 << pair.first << endl;
//	}
	outputFile1.close();
	outputFile2.close();
//	outputFile3.close();
	fai_destroy(fai);
	releaseSVItems(entryMap);
}
