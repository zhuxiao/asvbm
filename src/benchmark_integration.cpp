/*
 * benchmark_integration.cpp
 *
 *  Created on: Jan 5, 2025
 *      Author: node
 */
#include "benchmark_integration.h"
// 安全地将字符串转换为整数
int safeStoi(const string& str){
	if(str.empty()) return -1; // 空字符串返回 -1
	try{
		return stoi(str);
	}catch(const invalid_argument& e){
		return -1; // 无效输入返回 -1
	}catch(const out_of_range& e){
		return -1; // 超出范围返回 -1
	}
}

// 安全地将字符串转换为浮点数
float safeStof(const string& str){
	if(str.empty()) return -1.0f; // 空字符串返回 -1.0
	try{
		return stof(str);
	}catch(const invalid_argument& e){
		return -1.0f; // 无效输入返回 -1.0
	}catch(const out_of_range& e){
		return -1.0f; // 超出范围返回 -1.0
	}
}
string extractDP(string& line){
	istringstream iss(line);
	string field;
	int fieldCount = 0;
	while(getline(iss, field, '\t')){
		if(++fieldCount == 8){
			size_t svlenStart = field.find("DP=");
			if(svlenStart != string::npos){
				size_t svlenEnd = field.find(';', svlenStart);
				return (svlenEnd != string::npos) ? field.substr(svlenStart + strlen("DP="), svlenEnd - svlenStart - strlen("DP=")) : "";
			}
			break;
		}
	}
	return "";
}
string extractAF(string& line){
	istringstream iss(line);
	string field;
	int fieldCount = 0;
	while(getline(iss, field, '\t')){
		if(++fieldCount == 8){
			size_t svlenStart = field.find("AF=");
			if(svlenStart != string::npos){
				size_t svlenEnd = field.find(';', svlenStart);
				return (svlenEnd != string::npos) ? field.substr(svlenStart + strlen("AF="), svlenEnd - svlenStart - strlen("AF=")) : "";
			}
			break;
		}
	}
	return "";
}
//void sort(vector<SV_item*>& vec){
//	map<string, int> chrname_to_index;
//	int i;
//	for(i=1; i<=22; i++){
//		chrname_to_index[to_string(i)] = i;
//	}
//	chrname_to_index["X"] = 23;
//	chrname_to_index["Y"] = 24;
//	sort(vec.begin(), vec.end(),[&chrname_to_index](const SV_item* a, const SV_item* b){
//		int indexA = chrname_to_index.at(a->chrname);
//		int indexB = chrname_to_index.at(b->chrname);
//		if(indexA != indexB){
//			return indexA < indexB;
//		}
//		return a->startPos < b->startPos;
//	});
//}
vector<string> split1(const string& str, char delimiter){
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
string appendCustomInfo1(const string& originalInfo, const string& customInfo){
	if(originalInfo.empty()){
		return customInfo;
	}else{
		return originalInfo + ";" + customInfo;
	}
}
void releaseSVItems1(vector<SV_item*>& items){
	for(SV_item* item : items){
		delete item;
	}
	items.clear();
}
int calculateOverlap1(SV_item* item1, SV_item* item2){
	if(item1->sv_type==VAR_INS){
		item1->endPos=item1->startPos+item1->sv_len;
		item2->endPos=item2->startPos+item2->sv_len;
	}
	return (max(item1->startPos, item2->startPos) < min(item1->endPos, item2->endPos)) ? (min(item1->endPos, item2->endPos) - max(item1->startPos, item2->startPos)) : 0;
}
bool hasMatchingPairs1(const unordered_map<int, set<pair<int, int>>>& map, int index1, int index2){
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
void match1(vector<vector<SV_item*>>& sets, size_t l, vector<vector<int>>& subsetStartIndex, faidx_t *fai, size_t totalEntries, vector<int>& entriesStartIndex, vector<vector<int>>& connectedComponents){
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
	for(m=0; m<ToolsFilenames.size(); m++){
		for(n=m+1; n<ToolsFilenames.size(); n++){
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
						if(hasMatchingPairs1(map, index1, index2)){
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
											consistency = computeVarseqConsistency(item1, item2, fai, percentAlleleSeqIdentity);
											if(consistency >= percentAlleleSeqIdentity){
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
								if(entryMap1[index1] != nullptr && entryMap1[index2] != nullptr){
									if(index1 < entryMap1.size() && index2 < entryMap1.size()){
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
void addEntryMap1(size_t index, SV_item *entry){
	if(index >= 0 && entry != nullptr){
		entryMap1.push_back(entry);
	}else{
		cerr<< "Index out of bounds or null pointer in addEntryMap function: "<< index << endl;
	}
}
void integrationBenchmark(string &outputFile, string &ref_file){
	faidx_t *fai;
	fai = fai_load(ref_file.c_str());
	string newFilename = outputFile + "Benchmark.vcf";
	vector<int> entriesStartIndex(ToolsFilenames.size(), 0);
	vector<vector<int>> connectedComponents;
	pair<string, int> Pair;
	size_t i, j, m, n, l, k, o;
	size_t totalEntries = 0;
	int entryIndex = 0;
	size_t selectToolsNum;
	selectToolsNum = ceil(ToolsFilenames.size() * 0.55);

	for(i=0; i<ToolsFilenames.size(); i++){
		if(i > 0){
			entriesStartIndex[i] = entriesStartIndex[i - 1] + benchmark_vec[i - 1].size();
		}else{
			entriesStartIndex[i] = 0;
		}
		totalEntries += benchmark_vec[i].size();
		for(j=0; j<benchmark_vec[i].size(); j++){
			addEntryMap1(entryIndex, benchmark_vec[i][j]);
//			connectedComponents[entryIndex].push_back(entryMap[entryIndex]);
			entryIndex++;
		}
	}

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
	chrname_set1 = getChrnames(benchmark_vec[0]);
	chrname_set2 = getChrnames(benchmark_vec[1]);
	chrname_set_tmp = getChrUnion(chrname_set1, chrname_set2);
	chrname_vec_sorted = sortChrnames(chrname_set_tmp); // sort chromosome names
	subsets.resize(chrname_vec_sorted.size(), vector<vector<SV_item*>>());
	subsetStartIndex.resize(chrname_vec_sorted.size(), vector<int>(ToolsFilenames.size(), 0));

	for(i=0; i<chrname_vec_sorted.size(); i++){
		chr = chrname_vec_sorted.at(i);
		for(l=0; l<ToolsFilenames.size(); l++){
			subset = getItemsByChr(chr, benchmark_vec[l]);
			subsets[i].push_back(subset);
		}
	}
	for(i=1; i<chrname_vec_sorted.size(); i++){
		for(l=0; l<ToolsFilenames.size(); l++){
			subsetStartIndex[i][l] = subsetStartIndex[i-1][l] + subsets[i-1][l].size();
		}
	}
	vector<thread> threads2;
	for(l=0; l<chrname_vec_sorted.size(); l++){
		threads2.emplace_back(match1, ref(subsets[l]), l, ref(subsetStartIndex), fai, totalEntries, ref(entriesStartIndex), ref(connectedComponents));

	}
	for(vector<thread>::iterator it=threads2.begin(); it!=threads2.end(); it++){
		if (it->joinable()){
			it->join();
		}
	}

	//Write matched variants
		vector<size_t> connectedIndex;
		vector<int> differentToolsIndex;
		vector<int> mostRepresentativeVariantIndices;
		unordered_map<int, int> sizeDifferenceCount;
		ofstream outputFile1(newFilename);
		vector<string> fields;
//		vector<SV_item*> shared_fp_vec;
		string infoField, customInfo, newInfoField, newLineInfo, tools;
		SV_item* currentVariant;
		SV_item* mostRepresentativeVariant;
		SV_item* otherVariant;
		bool flag;
		int totalOverlapForCurrentVariant, maxTotalOverlap, mostRepresentativeVariantIndex, toolsNum, count, maxCount, len;
		outputFile1 << VCF_HEADER << endl;
		for(i=0; i<connectedComponents.size(); i++){
			if(connectedComponents[i].size() >= selectToolsNum){
//				if((extractDP(entryMap1[connectedComponents[i][0]]->lineInfo) == "" && extractAF(entryMap1[connectedComponents[i][0]]->lineInfo) == "") || (extractDP(entryMap1[connectedComponents[i][0]]->lineInfo) == "" && safeStof(extractAF(entryMap1[connectedComponents[i][0]]->lineInfo)) >= 0.2) || (safeStoi(extractDP(entryMap1[connectedComponents[i][0]]->lineInfo)) >= 25 && safeStoi(extractDP(entryMap1[connectedComponents[i][0]]->lineInfo)) <= 75 && extractAF(entryMap1[connectedComponents[i][0]]->lineInfo) == "")|| (safeStoi(extractDP(entryMap1[connectedComponents[i][0]]->lineInfo)) >= 25 && safeStoi(extractDP(entryMap1[connectedComponents[i][0]]->lineInfo)) <= 75 && safeStof(extractAF(entryMap1[connectedComponents[i][0]]->lineInfo)) >= 0.2)){
					differentToolsIndex.clear();
					for(j=0; j<connectedComponents[i].size(); j++){
						for(l=0; l<ToolsFilenames.size(); l++){
							if(l != (ToolsFilenames.size() - 1)){
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
						if(find(connectedIndex.begin(), connectedIndex.end(), connectedComponents[i][j]) == connectedIndex.end()){
							connectedIndex.push_back(connectedComponents[i][j]);
						}
					}


						if(flag){
							if(entryMap1[i] != nullptr){
								maxCount = 0;
								len = 0;
								while(maxCount == 0){
									len = len + 10;
									for(k=0; k<connectedComponents[i].size(); k++){
										differentToolsIndex.clear();
										for(j=0; j<connectedComponents[connectedComponents[i][k]].size(); j++){
											for(l=0; l<ToolsFilenames.size(); l++){
												if(l != (ToolsFilenames.size() - 1)){
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
											currentVariant = entryMap1[connectedComponents[i][k]];
											count=0;
											for(o=0; o<connectedComponents[i].size(); o++){
												if(k != o){
													otherVariant = entryMap1[connectedComponents[i][o]];
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
									currentVariant = entryMap1[mostRepresentativeVariantIndices[k]];
									totalOverlapForCurrentVariant = 0;
									for(o=0; o<mostRepresentativeVariantIndices.size(); o++){
										if(k != o){
											otherVariant = entryMap1[mostRepresentativeVariantIndices[o]];
											totalOverlapForCurrentVariant += calculateOverlap1(currentVariant,otherVariant);
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
									for(l=0; l<ToolsFilenames.size(); l++){
										if(l != (ToolsFilenames.size() - 1)){
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
								fields = split1(mostRepresentativeVariant->lineInfo, '\t');
								infoField = fields[7];
								customInfo = "TOOLSNUM=" + to_string(toolsNum) + ";TOOLS=" + tools;
								newInfoField = appendCustomInfo1(infoField, customInfo);
								fields[7] = newInfoField;
								newLineInfo = "";
								for(n=0; n<fields.size(); n++){
									newLineInfo += fields[n] + (n < fields.size() - 1 ? "\t" : "");
								}
								Pair = make_pair(newLineInfo,1);
								benchmark_hg002_Map.insert(Pair);
							}else{

							}
						}
					}
//				}
			}
		}

		vector<pair<string, int>> sortedVec(benchmark_hg002_Map.begin(), benchmark_hg002_Map.end());
//		sort(sortedVec.begin(), sortedVec.end(), numericSort);
		try{
			sort(sortedVec.begin(), sortedVec.end(), numericSort);
		}catch(const exception& e){
			cerr << "Error: " << e.what() << endl;
		}
		for(const auto& pair : sortedVec){
			outputFile1 << pair.first << endl;
		}
		outputFile1.close();
		fai_destroy(fai);
		releaseSVItems1(entryMap1);
}
