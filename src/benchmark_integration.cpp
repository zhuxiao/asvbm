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
vector<int> extractAlleles(const string& genotype){
	vector<int> alleles;
	stringstream ss(genotype);
	string part;

	while(getline(ss, part, '/') || getline(ss, part, '|')){
		if(part == "."){
			alleles.push_back(-1);
		}else{
			alleles.push_back(std::stoi(part));
		}
	}

	return alleles;
}
vector<vector<int>> parseGenotypes(const string& vcfLine){
	vector<vector<int>> genotypes;

	if(vcfLine.empty() || vcfLine[0] == '#'){
		return genotypes;
	}
	stringstream ss(vcfLine);
	string field;
	int i, currentIndex;
	size_t gtIndex;
	for(i = 0; i < 8; i++){
		if(!getline(ss, field, '\t')){
			cerr << "Error: Incomplete VCF line, missing required fields" << std::endl;
			return genotypes;
		}
	}

	gtIndex = SIZE_MAX;
	currentIndex = 0;
	if(getline(ss, field, '\t')){
		stringstream formatSS(field);
		string formatField;
		while(getline(formatSS, formatField, ':')){
			if(formatField == "GT"){
				gtIndex = currentIndex;
				break;
			}
			currentIndex++;
		}
	}

	if(gtIndex == SIZE_MAX){
		cerr << "Warning: No GT field found in FORMAT" << std::endl;
		return genotypes;
	}

	while(getline(ss, field, '\t')){
		vector<string> sampleValues;
		stringstream sampleSS(field);
		string value;
		while(getline(sampleSS, value, ':')){
			sampleValues.push_back(value);
		}
		if(gtIndex < sampleValues.size()){
			genotypes.push_back(extractAlleles(sampleValues[gtIndex]));
		}else{
			genotypes.push_back({-1, -1});
		}
	}
	return genotypes;
}
bool checkMendelianLaw(const vector<int>& sonAlleles, const vector<int>& fatherAlleles, const vector<int>& motherAlleles){
	if(sonAlleles.size() != 2 || fatherAlleles.size() != 2 || motherAlleles.size() != 2){
		cerr << "Error: Invalid genotype format, expected diploid (2 alleles)" << endl;
		return false;
	}
	bool sonHasMissing, fatherHasMissing, motherHasMissing, allele1FromFather, allele1FromMother, allele2FromFather, allele2FromMother, validInheritance;
	sonHasMissing = (sonAlleles[0] == -1 || sonAlleles[1] == -1);
	fatherHasMissing = (fatherAlleles[0] == -1 || fatherAlleles[1] == -1);
	motherHasMissing = (motherAlleles[0] == -1 || motherAlleles[1] == -1);

	if(!sonHasMissing && !fatherHasMissing && !motherHasMissing){
		allele1FromFather = false, allele1FromMother = false, allele2FromFather = false, allele2FromMother = false;
		for(int fatherAllele : fatherAlleles){
			if(fatherAllele == sonAlleles[0]){
				allele1FromFather = true;
				break;
			}
		}
		for(int motherAllele : motherAlleles){
			if(motherAllele == sonAlleles[0]){
				allele1FromMother = true;
				break;
			}
		}

		for(int fatherAllele : fatherAlleles){
			if(fatherAllele == sonAlleles[1]){
				allele2FromFather = true;
				break;
			}
		}
		for(int motherAllele : motherAlleles){
			if(motherAllele == sonAlleles[1]){
				allele2FromMother = true;
				break;
			}
		}
		validInheritance = (allele1FromFather && allele2FromMother) || (allele1FromMother && allele2FromFather);
		return validInheritance;
	}else{
		return true;
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
	size_t i, j, m, n, k;
	vector<size_t> overlap_type_vec;
	vector<SV_item*> subset1, subset2, merge_item1_vec, merge_item2_vec;
	SV_item* Definite_Merge_item;
	int64_t start_idx, end_idx, new_start_idx, start_search_pos, end_search_pos;
	int64_t start_idx_tmp0, end_idx_tmp0, sv_len1, sv_len2, feature_len;
	float svlenRatio_tmp1, svlenRatio;
	bool flag, foundIndex;
	double consistency;
	size_t index1, index2;
	pair<int, int> connectedPair;
	unordered_map<int, set<pair<int, int>>> map;
	for(m=0; m<ToolsFilenames.size(); m++){
		for(n=m+1; n<ToolsFilenames.size(); n++){
			if(m == n) continue;
			subset1=sets[m];
			subset2=sets[n];
			for(i=0; i<subset1.size(); i++) {subset1.at(i)->overlapped = false; subset1.at(i)->seqcons = "-"; subset1.at(i)->mergeMark = false; subset1.at(i)->mergeMark2 = false;subset1.at(i)-> LpFlag = false;}
			for(i=0; i<subset2.size(); i++) {subset2.at(i)->overlapped = false;	subset2.at(i)->seqcons = "-"; subset2.at(i)->mergeMark = false; subset2.at(i)->mergeMark2 = false;subset2.at(i)-> LpFlag = false;}
			SV_item *item1, *item2;
			start_idx = end_idx = -1;
			for(i=0; i<subset1.size(); i++){
				item1 = subset1.at(i);
//				item1_merge_flag = false;
//				item2_merge_flag = false;
				start_idx_tmp0 = start_idx;
				end_idx_tmp0 = end_idx;
//				if(abs(item1->sv_len) <= 100){
//					start_search_pos = item1->startPos - tolerance;
//					end_search_pos = item1->endPos + tolerance;
//				}else{
//					start_search_pos = item1->startPos - abs(item1->sv_len);
//					end_search_pos = item1->endPos + abs(item1->sv_len);
//				}
				start_search_pos = item1->startPos - 1000;
				end_search_pos = item1->endPos + 1000;
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

//					if(item1->mergeMark == true){
//						continue;
//					}
					for(j=start_idx; j<subset2.size(); j++){
						item2 = subset2.at(j);
						if(item2->startPos < start_search_pos)
							continue;
//						if(item1_merge_flag){
//							continue;
//						}
//						if(item2_merge_flag){
//							continue;
//						}
						overlap_type_vec = computeOverlapType(item1, item2);
						if((overlap_type_vec.size() > 0 && overlap_type_vec.at(0) != NO_OVERLAP)){
							if(item1->sv_type == item2->sv_type){
								IsmergeJudge2(i, item1, item2, subset1, fai, merge_item1_vec);
							}

							if(!merge_item1_vec.empty()){
								index2 = entriesStartIndex[n] + subsetStartIndex[l][n] + j;
								for(k=0; k < merge_item1_vec.size(); k++){
									Definite_Merge_item = merge_item1_vec.at(k);
									index1 = itemToIndex[Definite_Merge_item];
									connectedComponents[index1].push_back(index2);
									connectedComponents[index2].push_back(index1);
								}
								merge_item1_vec.clear();
//								if(item1->mergeMark == true) {item1_merge_flag = true; continue;}
							}
							if(item1->sv_type == item2->sv_type){
								IsmergeJudge4(j, item2, item1, subset2, fai, merge_item2_vec);
							}

							if(!merge_item2_vec.empty()){
								index1 = entriesStartIndex[m] + subsetStartIndex[l][m] + i;
								for(k=0; k < merge_item2_vec.size(); k++){
									Definite_Merge_item = merge_item2_vec.at(k);
									index2 = itemToIndex[Definite_Merge_item];
									connectedComponents[index1].push_back(index2);
									connectedComponents[index2].push_back(index1);
								}
								merge_item2_vec.clear();
//								if(item2->mergeMark2 == true) {item2_merge_flag = true; continue;}
							}
						}
						flag = false;
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
								feature_len = max(sv_len1, sv_len2);
								if(feature_len < 50){
									svlenRatio = 0;
								}else if(feature_len >= 1000){
									svlenRatio= 0.7;
								}else{
									svlenRatio = 0.7 * (1 - exp(-(feature_len - 50) / 100));
//									svlenRatio = 0.7 * (1 - exp(-0.006 * (feature_len - 50)));
//									svlenRatio = (0.7 / 950) * (feature_len - 50);
								}

								if(svlenRatio_tmp1 >= svlenRatio){
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
//											consistency = computeVarseqConsistency(item1, item2, fai, percentAlleleSeqSim);
//											if(consistency >= percentAlleleSeqSim){
//											if(overlap_type_vec.size() <= 0 || overlap_type_vec.at(0) == NO_OVERLAP){
//												consistency = 0;
//											}else{
//											consistency = computeVarseqConsistency(item1, item2, fai, AlleleSeqSim);
//											}
											consistency = computeVarseqConsistency(item1, item2, fai, 0.7);
											if(consistency >= 0.7){
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
							index1 = entriesStartIndex[m] + subsetStartIndex[l][m] + i;
							index2 = entriesStartIndex[n] + subsetStartIndex[l][n] + j;
							if(index1 < totalEntries && index2 < totalEntries){
								if(entryMap1[index1] != nullptr && entryMap1[index2] != nullptr){
									if(index1 < entryMap1.size() && index2 < entryMap1.size()){
										if(find(connectedComponents[index1].begin(),connectedComponents[index1].end(),index2) == connectedComponents[index1].end()){
											connectedComponents[index1].push_back(index2);
										}
										if(find(connectedComponents[index2].begin(),connectedComponents[index2].end(),index1) == connectedComponents[index2].end()){
											connectedComponents[index2].push_back(index1);
										}

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
void writeMatchedVariant(size_t subsetSize, size_t selectToolsNum, vector<vector<SV_item*>>& sets, vector<int>& entriesStartIndex, vector<vector<int>>& connectedComponents){
	size_t i, j, m, n, l, k, o;
	unordered_set<size_t> connectedIndex;
	vector<int> differentToolsIndex;
	vector<int> mostRepresentativeVariantIndices, connectedComponentsToAdd;
	unordered_map<int, int> sizeDifferenceCount;
	vector<string> fields;
	pair<string, int> Pair;
//		vector<SV_item*> shared_fp_vec;
	string infoField, customInfo, newInfoField, newLineInfo, tools;
	SV_item* currentVariant;
	SV_item* mostRepresentativeVariant;
	SV_item* otherVariant;
	bool flag;
	int totalOverlapForCurrentVariant, maxTotalOverlap, mostRepresentativeVariantIndex, toolsNum, count, maxCount, len, max_len, maxSVLen;
	for(i=0; i<subsetSize; i++){
			if(connectedComponents[i].size() >= selectToolsNum){
				max_len = 0;
	//				mergeFlagExist = false;
			////	if((extractDP(entryMap1[connectedComponents[i][0]]->lineInfo) == "" && extractAF(entryMap1[connectedComponents[i][0]]->lineInfo) == "") || (extractDP(entryMap1[connectedComponents[i][0]]->lineInfo) == "" && safeStof(extractAF(entryMap1[connectedComponents[i][0]]->lineInfo)) >= 0.2) || (safeStoi(extractDP(entryMap1[connectedComponents[i][0]]->lineInfo)) >= 25 && safeStoi(extractDP(entryMap1[connectedComponents[i][0]]->lineInfo)) <= 75 && extractAF(entryMap1[connectedComponents[i][0]]->lineInfo) == "")|| (safeStoi(extractDP(entryMap1[connectedComponents[i][0]]->lineInfo)) >= 25 && safeStoi(extractDP(entryMap1[connectedComponents[i][0]]->lineInfo)) <= 75 && safeStof(extractAF(entryMap1[connectedComponents[i][0]]->lineInfo)) >= 0.2)){
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
						if(differentToolsIndex.size() >= selectToolsNum) break;
					}

					if(differentToolsIndex.size() >= selectToolsNum){
						toolsNum = 0;
						tools = "";

					flag = true;

	//					if(!connectedIndex.empty()){
	//						for(m=0; m<connectedIndex.size(); m++){
	//							if(i == connectedIndex[m]){
	//								flag = false;
	//								break;
	//							}
	//						}
	//					}
					if(!connectedIndex.empty() && connectedIndex.find(i) != connectedIndex.end()){
						flag = false;
					}
	//					for(j=0; j<connectedComponents[i].size(); j++){
	//						if(find(connectedIndex.begin(), connectedIndex.end(), connectedComponents[i][j]) == connectedIndex.end()){
	//							connectedIndex.push_back(connectedComponents[i][j]);
	//						}
	//					}
	//					for(j=0; j<connectedComponents[i].size(); j++){
	//						for(m=0; m<connectedComponents[connectedComponents[i][j]].size(); m++){
	//							if(find(connectedIndex.begin(), connectedIndex.end(), connectedComponents[connectedComponents[i][j]][m]) == connectedIndex.end()){
	//								connectedIndex.push_back(connectedComponents[connectedComponents[i][j]][m]);
	//							}
	//						}
	//					}

					for(j=0; j<connectedComponents[i].size(); j++){
						connectedComponentsToAdd.clear();
						for(m=0; m<connectedComponents[connectedComponents[i][j]].size(); m++){
							if(connectedIndex.find(connectedComponents[connectedComponents[i][j]][m]) == connectedIndex.end()){
								connectedIndex.insert(connectedComponents[connectedComponents[i][j]][m]);
							}

							if(j >= 1 && m >= 1){
								if(find(connectedComponents[i].begin(),connectedComponents[i].end(),connectedComponents[connectedComponents[i][j]][m]) == connectedComponents[i].end()){
									connectedComponentsToAdd.push_back(connectedComponents[connectedComponents[i][j]][m]);
								}
							}

						}
	//						if(entryMap1[connectedComponents[i][j]]->mergeFlag){mergeFlagExist = true;}
						connectedComponents[i].insert(connectedComponents[i].end(), connectedComponentsToAdd.begin(), connectedComponentsToAdd.end());
					}


						if(flag){
							if(entryMap1[i] != nullptr){
								for(k=0; k<connectedComponents[i].size(); k++){
									if(entryMap1[connectedComponents[i][k]]->sv_len > max_len)
										max_len = entryMap1[connectedComponents[i][k]]->sv_len;
								}

									maxCount = 0;
									len = 0;
									while(maxCount == 0){
										len = len + 10;
										for(k=0; k<connectedComponents[i].size(); k++){
	/*											differentToolsIndex.clear();
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
											}*/
	//											if(differentToolsIndex.size() >= selectToolsNum){
	//												if(entryMap1[connectedComponents[i][k]]->mergeFlag) continue;
												if(entryMap1[connectedComponents[i][k]]->sv_len > (max_len * 0.7)){

												currentVariant = entryMap1[connectedComponents[i][k]];
												count=0;
												for(o=0; o<connectedComponents[i].size(); o++){
	//													if(k != o && !entryMap1[connectedComponents[i][o]]->mergeFlag && entryMap1[connectedComponents[i][o]]->sv_len > (max_len * 0.7)){
													if(k != o && entryMap1[connectedComponents[i][o]]->sv_len > (max_len * 0.7)){
														otherVariant = entryMap1[connectedComponents[i][o]];
														if(abs(currentVariant->sv_len - otherVariant->sv_len) <= len){
															count++;
														}
													}
												}
												sizeDifferenceCount.insert(make_pair(connectedComponents[i][k], count));
												}
	//											}
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
	//									sizeDifferenceCount.clear();
	//									if(connectedComponents[i].size() == 1) break;
										if(sizeDifferenceCount.size() == 1 || sizeDifferenceCount.size() == 0){sizeDifferenceCount.clear();break;}
										else{sizeDifferenceCount.clear();}
										if(maxCount == 0) mostRepresentativeVariantIndices.clear();

									}



								maxTotalOverlap = -1, maxSVLen = 0;
								for(k=0; k<mostRepresentativeVariantIndices.size(); k++){
									currentVariant = entryMap1[mostRepresentativeVariantIndices[k]];
									totalOverlapForCurrentVariant = 0;

									for(o=0; o<mostRepresentativeVariantIndices.size(); o++){
										if(k != o){
											otherVariant = entryMap1[mostRepresentativeVariantIndices[o]];
											totalOverlapForCurrentVariant += calculateOverlap1(currentVariant,otherVariant);
										}

									}
									if(totalOverlapForCurrentVariant > maxTotalOverlap){
										maxSVLen = currentVariant->sv_len;
										maxTotalOverlap = totalOverlapForCurrentVariant;
										mostRepresentativeVariant = currentVariant;
										mostRepresentativeVariantIndex = mostRepresentativeVariantIndices[k];
									}else if(totalOverlapForCurrentVariant == maxTotalOverlap && currentVariant->sv_len >= maxSVLen){
										maxSVLen = currentVariant->sv_len;
										maxTotalOverlap = totalOverlapForCurrentVariant;
										mostRepresentativeVariant = currentVariant;
										mostRepresentativeVariantIndex = mostRepresentativeVariantIndices[k];
									}
								}
								mostRepresentativeVariantIndices.clear();

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
								for(k=0; k<connectedComponents[i].size(); k++){
									for(l=0; l<ToolsFilenames.size(); l++){
										if(l != (ToolsFilenames.size() - 1)){
											if(connectedComponents[i][k] < entriesStartIndex[l + 1]){
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
				//}
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
set<string> readVCFHeader(const string& filename){
	set<string> headers;
	ifstream file(filename);
	string line;

	if(!file.is_open()){
		cerr << "无法打开文件: " << filename << endl;
		return headers;
	}
	while(getline(file, line)){
		if(line.empty() || line[0] != '#'){
			break;
		}
		if(line.substr(0, 2) == "##"){
			headers.insert(line);
		}
	}

	file.close();
	return headers;
}
string getCurrentTime(){
	auto now = chrono::system_clock::now();
	auto in_time_t = chrono::system_clock::to_time_t(now);
	stringstream ss;
	ss << put_time(localtime(&in_time_t), "%Y-%m-%dT%H:%M:%S%z");
	return ss.str();
}
vector<string> mergeHeaders(const vector<string>& filenames){
	set<string> mergedHeaders;
//	map<string, string> uniqueHeaders;
	map<string, map<string, string>> uniqueHeadersByType;
	set<string> sources, headers, sampleNames;
	string reference, combinedSource, sv_caller, source, id, currentFileDate, line, sample, sampleHeader, toolsNumInfo, toolsInfo;
	size_t pos, eq_pos, end;
	vector<string> orderedHeaders;
	toolsNumInfo = "##INFO=<ID=TOOLSNUM,Number=1,Type=Integer,Description=\"The number of tools for identifying variants that match this variant\">";
	toolsInfo = "##INFO=<ID=TOOLS,Number=1,Type=String,Description=\"Tools for identifying variants that match the variant\">";
	for(const auto& filename : filenames){
		headers = readVCFHeader(filename);
		for(const auto& header : headers){
			if(header.find("##fileformat=") != string::npos){
				mergedHeaders.insert(header);
			}else if(header.find("##reference=") != string::npos){
				if(reference.empty()){
					reference = header;
				}
//				mergedHeaders.insert(header);
			}else if(header.find("##contig=") != string::npos){
				mergedHeaders.insert(header);
			}else if(header.find("##fileDate=") != string::npos){
				continue;
			}else if(header.find("##source=") != string::npos){
				pos = header.find("##source=");
				if(pos != string::npos){
					eq_pos = header.find('=', pos);
					if(eq_pos != string::npos){
						source = header.substr(eq_pos + 1);
						sources.insert(source);
					}
				}
			}else if(header.find("##SVCaller=") != string::npos){
				pos = header.find("##SVCaller=");
				if(pos != string::npos){
					eq_pos = header.find('=', pos);
					if(eq_pos != string::npos){
						sv_caller = header.substr(eq_pos + 1);
						sources.insert(sv_caller);
					}
				}
			}else if (header.find("##INFO=") != string::npos || header.find("##FORMAT=") != string::npos || header.find("##FILTER=") != string::npos || header.find("##ALT=") != string::npos){
				string type;
				if(header.find("##INFO=") != string::npos) type = "INFO";
				else if(header.find("##FORMAT=") != string::npos) type = "FORMAT";
				else if(header.find("##FILTER=") != string::npos) type = "FILTER";
				else if(header.find("##ALT=") != string::npos) type = "ALT";

				pos = header.find("ID=");
				if(pos != string::npos){
					end = header.find(',', pos);
					if(end == string::npos){
						end = header.size();
					}
					id = header.substr(pos + 3, end - pos - 3);
					if(uniqueHeadersByType[type].find(id) == uniqueHeadersByType[type].end()){
						uniqueHeadersByType[type][id] = header;
					}
				}
			}
//			else if(header.find("##INFO=") != string::npos || header.find("##FORMAT=") != string::npos || header.find("##FILTER=") != string::npos || header.find("##ALT=") != string::npos){
//				pos = header.find("ID=");
//				if(pos != string::npos){
//					end = header.find(',', pos);
//					if(end == string::npos){
//						end = header.size();
//					}
//					id = header.substr(pos + 3, end - pos - 3);
//					if(uniqueHeaders.find(id) == uniqueHeaders.end()){
//						uniqueHeaders[id] = header;
//					}
//				}
//			}
		}
		ifstream file(filename);
		while(getline(file, line)){
			if(line.empty() || line[0] != '#'){
				break;
			}
			if(line[0] == '#' && line[1] != '#'){
				istringstream iss(line);
				for(int i = 0; i < 9; ++i){
					getline(iss, sample, '\t');
				}
				while(getline(iss, sample, '\t')){
					sampleNames.insert(sample);
				}
				break;
			}
		}
		file.close();
	}
//	for(const auto& pair : uniqueHeaders){
//		mergedHeaders.insert(pair.second);
//	}

	for(const auto& typePair : uniqueHeadersByType){
		for(const auto& idPair : typePair.second){
			mergedHeaders.insert(idPair.second);
		}
	}
	currentFileDate = getCurrentTime();
	headers.insert("##fileDate=" + currentFileDate);
	if(!reference.empty()){
		mergedHeaders.insert(reference);
	}
	if(!sources.empty()){
		combinedSource = "##source=";
		for(const auto& source : sources){
			combinedSource += source + ";";
		}
		combinedSource.pop_back();
		mergedHeaders.insert(combinedSource);
	}
	mergedHeaders.insert(toolsNumInfo);
	mergedHeaders.insert(toolsInfo);
	for(const auto& header : mergedHeaders){
		if(header.find("##fileformat=") != string::npos){
			orderedHeaders.push_back(header);
		}
	}
	for(const auto& header : mergedHeaders){
		if(header.find("##fileDate=") != string::npos){
			orderedHeaders.push_back(header);
		}
	}
	for(const auto& header : mergedHeaders){
		if(header.find("##source=") != string::npos){
			orderedHeaders.push_back(header);
		}
	}
	for(const auto& header : mergedHeaders){
		if(header.find("##reference=") != string::npos){
			orderedHeaders.push_back(header);
		}
	}
	for(const auto& header : mergedHeaders){
		if(header.find("##contig=") != string::npos){
			orderedHeaders.push_back(header);
		}
	}
	for(const auto& header : mergedHeaders){
		if(header.find("##INFO=") != string::npos){
			orderedHeaders.push_back(header);
		}
	}
	for(const auto& header : mergedHeaders){
		if(header.find("##FILTER=") != string::npos){
			orderedHeaders.push_back(header);
		}
	}
	for(const auto& header : mergedHeaders){
		if(header.find("##FORMAT=") != string::npos){
			orderedHeaders.push_back(header);
		}
	}
	for(const auto& header : mergedHeaders){
		if(header.find("##ALT=") != string::npos){
			orderedHeaders.push_back(header);
		}
	}

//	for(const auto& sample : sampleNames){
//		if (!sampleHeader.empty()){
//			sampleHeader += "_";
//		}
//		sampleHeader += sample;
//	}
	sampleHeader="sample";
	orderedHeaders.push_back("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sampleHeader);
	return orderedHeaders;
}

void integrationBenchmark(string &outputFile, string &ref_file, vector<string> &sv_files1){
	vector<string> headers = mergeHeaders(sv_files1);
	faidx_t *fai;
	string newFilename;
	fai = fai_load(ref_file.c_str());

	vector<int> entriesStartIndex(ToolsFilenames.size(), 0);
	vector<vector<int>> connectedComponents;
	pair<string, int> Pair;
	size_t i, j, l, m ,n ,k, o;
	size_t totalEntries = 0;
	int entryIndex = 0;
	size_t selectToolsNum;
	size_t input, input2;
	bool validInput;
//	selectToolsNum = ceil(ToolsFilenames.size() * 0.6);
	cout << "Do you want to set the number of selected tools? (1 for Yes, 0 for No): " << endl;
	cout << "If you input 0, the selectToolsNum = " << ceil(sharedFPFilenames.size() * 0.6) << endl;

	while(true){
		cin >> input;
		if(cin.fail() || input < 0 || input > 1){
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cout << "Invalid choice. Please enter 1 for Yes or 0 for No." << endl;
		}else{
			break;
		}
	}
	if(input == 1){
		validInput = false;
		while(!validInput){
			cout << "Set your number (between 1 and " << sharedFPFilenames.size() << "): ";
			cin >> input2;
			if(cin.fail() || input2 < 1 || input2 > sharedFPFilenames.size()){
				cin.clear();
				cin.ignore(numeric_limits<streamsize>::max(), '\n');
				cout << "Invalid input. Please try again." << endl;
			}else{
				validInput = true;
				selectToolsNum = input2;
			}
		}
	}else if(input == 0){
		selectToolsNum = ceil(sharedFPFilenames.size() * 0.6);
	}
	cout << "Selected value: " << selectToolsNum << endl;

	if(!genotype_filterBench_flag){
		newFilename = outputFile + to_string(ToolsFilenames.size()) + "_" + to_string(selectToolsNum) + "_" +"Benchmark.vcf";
	}else{
		newFilename = outputFile + "genotype_filteredBenchmark.vcf";
	}

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
			itemToIndex[benchmark_vec[i][j]] = entryIndex;
			entryIndex++;
//			benchmark_vec[i][j]->mergeFlag = false;
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
	size_t size;
	vector<size_t> subsetsize;
	vector<vector<int>> subsetStartIndex;
	chrname_set1 = getChrnames(benchmark_vec[0]);
	chrname_set2 = getChrnames(benchmark_vec[1]);
	chrname_set_tmp = getChrUnion(chrname_set1, chrname_set2);
	chrname_vec_sorted = sortChrnames(chrname_set_tmp); // sort chromosome names
	subsets.resize(chrname_vec_sorted.size(), vector<vector<SV_item*>>());
	subsetStartIndex.resize(chrname_vec_sorted.size(), vector<int>(ToolsFilenames.size(), 0));

	for(i=0; i<chrname_vec_sorted.size(); i++){
		chr = chrname_vec_sorted.at(i);
		size = 0;
		for(l=0; l<ToolsFilenames.size(); l++){
			subset = getItemsByChr(chr, benchmark_vec[l]);
			subsets[i].push_back(subset);
			size = size + subset.size();
		}
		subsetsize.push_back(size);
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
//		vector<size_t> connectedIndex;
		unordered_set<size_t> connectedIndex;
		vector<int> differentToolsIndex;
		vector<int> mostRepresentativeVariantIndices, connectedComponentsToAdd;
		unordered_map<int, int> sizeDifferenceCount;
		ofstream outputFile1(newFilename);
		vector<string> fields;
//		vector<SV_item*> shared_fp_vec;
		string infoField, customInfo, newInfoField, newLineInfo, tools;
		SV_item* currentVariant;
		SV_item* mostRepresentativeVariant;
		SV_item* otherVariant;
		bool flag;
		int totalOverlapForCurrentVariant, maxTotalOverlap, mostRepresentativeVariantIndex, toolsNum, count, maxCount, len, max_len, maxSVLen;
//		outputFile1 << VCF_HEADER << endl;

		for(const auto& header : headers){
			outputFile1 << header << endl;
		}
		if(!genotype_filterBench_flag){

/*			vector<thread> threads3;
			for(l=0; l<chrname_vec_sorted.size(); l++){
				threads3.emplace_back(writeMatchedVariant, ref(subsetsize[l]), selectToolsNum, ref(subsets[l]),ref(entriesStartIndex), ref(connectedComponents));

			}
			for(vector<thread>::iterator it=threads3.begin(); it!=threads3.end(); it++){
				if (it->joinable()){
					it->join();
				}
			}*/
			for(i=0; i<connectedComponents.size(); i++){
					if(connectedComponents[i].size() >= selectToolsNum){
						max_len = 0;
			//				mergeFlagExist = false;
					////	if((extractDP(entryMap1[connectedComponents[i][0]]->lineInfo) == "" && extractAF(entryMap1[connectedComponents[i][0]]->lineInfo) == "") || (extractDP(entryMap1[connectedComponents[i][0]]->lineInfo) == "" && safeStof(extractAF(entryMap1[connectedComponents[i][0]]->lineInfo)) >= 0.2) || (safeStoi(extractDP(entryMap1[connectedComponents[i][0]]->lineInfo)) >= 25 && safeStoi(extractDP(entryMap1[connectedComponents[i][0]]->lineInfo)) <= 75 && extractAF(entryMap1[connectedComponents[i][0]]->lineInfo) == "")|| (safeStoi(extractDP(entryMap1[connectedComponents[i][0]]->lineInfo)) >= 25 && safeStoi(extractDP(entryMap1[connectedComponents[i][0]]->lineInfo)) <= 75 && safeStof(extractAF(entryMap1[connectedComponents[i][0]]->lineInfo)) >= 0.2)){
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
								if(differentToolsIndex.size() >= selectToolsNum) break;
							}

							if(differentToolsIndex.size() >= selectToolsNum){
								toolsNum = 0;
								tools = "";

							flag = true;

			//					if(!connectedIndex.empty()){
			//						for(m=0; m<connectedIndex.size(); m++){
			//							if(i == connectedIndex[m]){
			//								flag = false;
			//								break;
			//							}
			//						}
			//					}
							if(!connectedIndex.empty() && connectedIndex.find(i) != connectedIndex.end()){
								flag = false;
							}
			//					for(j=0; j<connectedComponents[i].size(); j++){
			//						if(find(connectedIndex.begin(), connectedIndex.end(), connectedComponents[i][j]) == connectedIndex.end()){
			//							connectedIndex.push_back(connectedComponents[i][j]);
			//						}
			//					}
			//					for(j=0; j<connectedComponents[i].size(); j++){
			//						for(m=0; m<connectedComponents[connectedComponents[i][j]].size(); m++){
			//							if(find(connectedIndex.begin(), connectedIndex.end(), connectedComponents[connectedComponents[i][j]][m]) == connectedIndex.end()){
			//								connectedIndex.push_back(connectedComponents[connectedComponents[i][j]][m]);
			//							}
			//						}
			//					}

							for(j=0; j<connectedComponents[i].size(); j++){
								connectedComponentsToAdd.clear();
								for(m=0; m<connectedComponents[connectedComponents[i][j]].size(); m++){
									if(connectedIndex.find(connectedComponents[connectedComponents[i][j]][m]) == connectedIndex.end()){
										connectedIndex.insert(connectedComponents[connectedComponents[i][j]][m]);
									}

									if(j >= 1 && m >= 1){
										if(find(connectedComponents[i].begin(),connectedComponents[i].end(),connectedComponents[connectedComponents[i][j]][m]) == connectedComponents[i].end()){
											connectedComponentsToAdd.push_back(connectedComponents[connectedComponents[i][j]][m]);
										}
									}

								}
			//						if(entryMap1[connectedComponents[i][j]]->mergeFlag){mergeFlagExist = true;}
								connectedComponents[i].insert(connectedComponents[i].end(), connectedComponentsToAdd.begin(), connectedComponentsToAdd.end());
							}


								if(flag){
									if(entryMap1[i] != nullptr){
										for(k=0; k<connectedComponents[i].size(); k++){
											if(entryMap1[connectedComponents[i][k]]->sv_len > max_len)
												max_len = entryMap1[connectedComponents[i][k]]->sv_len;
										}

											maxCount = 0;
											len = 0;
											while(maxCount == 0){
												len = len + 10;
												for(k=0; k<connectedComponents[i].size(); k++){
			/*											differentToolsIndex.clear();
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
													}*/
			//											if(differentToolsIndex.size() >= selectToolsNum){
			//												if(entryMap1[connectedComponents[i][k]]->mergeFlag) continue;
														if(entryMap1[connectedComponents[i][k]]->sv_len > (max_len * 0.7)){

														currentVariant = entryMap1[connectedComponents[i][k]];
														count=0;
														for(o=0; o<connectedComponents[i].size(); o++){
			//													if(k != o && !entryMap1[connectedComponents[i][o]]->mergeFlag && entryMap1[connectedComponents[i][o]]->sv_len > (max_len * 0.7)){
															if(k != o && entryMap1[connectedComponents[i][o]]->sv_len > (max_len * 0.7)){
																otherVariant = entryMap1[connectedComponents[i][o]];
																if(abs(currentVariant->sv_len - otherVariant->sv_len) <= len){
																	count++;
																}
															}
														}
														sizeDifferenceCount.insert(make_pair(connectedComponents[i][k], count));
														}
			//											}
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
			//									sizeDifferenceCount.clear();
			//									if(connectedComponents[i].size() == 1) break;
												if(sizeDifferenceCount.size() == 1 || sizeDifferenceCount.size() == 0){sizeDifferenceCount.clear();break;}
												else{sizeDifferenceCount.clear();}
												if(maxCount == 0) mostRepresentativeVariantIndices.clear();

											}



										maxTotalOverlap = -1, maxSVLen = 0;
										for(k=0; k<mostRepresentativeVariantIndices.size(); k++){
											currentVariant = entryMap1[mostRepresentativeVariantIndices[k]];
											totalOverlapForCurrentVariant = 0;

											for(o=0; o<mostRepresentativeVariantIndices.size(); o++){
												if(k != o){
													otherVariant = entryMap1[mostRepresentativeVariantIndices[o]];
													totalOverlapForCurrentVariant += calculateOverlap1(currentVariant,otherVariant);
												}

											}
											if(totalOverlapForCurrentVariant > maxTotalOverlap){
												maxSVLen = currentVariant->sv_len;
												maxTotalOverlap = totalOverlapForCurrentVariant;
												mostRepresentativeVariant = currentVariant;
												mostRepresentativeVariantIndex = mostRepresentativeVariantIndices[k];
											}else if(totalOverlapForCurrentVariant == maxTotalOverlap && currentVariant->sv_len >= maxSVLen){
												maxSVLen = currentVariant->sv_len;
												maxTotalOverlap = totalOverlapForCurrentVariant;
												mostRepresentativeVariant = currentVariant;
												mostRepresentativeVariantIndex = mostRepresentativeVariantIndices[k];
											}
										}
										mostRepresentativeVariantIndices.clear();

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
										for(k=0; k<connectedComponents[i].size(); k++){
											for(l=0; l<ToolsFilenames.size(); l++){
												if(l != (ToolsFilenames.size() - 1)){
													if(connectedComponents[i][k] < entriesStartIndex[l + 1]){
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
						//}
					}
				}
		}else{
			vector<vector<int>> sonGts, fatherGts, motherGts;
			bool satisfiesLaw;
			for(i=0; i<benchmark_vec[0].size(); i++){
				if(connectedComponents[i].size() >= selectToolsNum){
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
						sonGts = parseGenotypes(entryMap1[connectedComponents[i][0]]->lineInfo);
						for(j=1; j<connectedComponents[i].size(); j++){
							if(connectedComponents[i][j] >= entriesStartIndex[1] && connectedComponents[i][j] < entriesStartIndex[2]){
								fatherGts = parseGenotypes(entryMap1[connectedComponents[i][1]]->lineInfo);
							}
							if(connectedComponents[i][j] > entriesStartIndex[2]){
								motherGts = parseGenotypes(entryMap1[connectedComponents[i][2]]->lineInfo);
							}
						}


						if(sonGts.empty() || fatherGts.empty() || motherGts.empty()){
							continue;
						}
//						if(fatherGts == motherGts) continue;
						satisfiesLaw = checkMendelianLaw(sonGts[0], fatherGts[0], motherGts[0]);
						if(!satisfiesLaw){
							Pair = make_pair(entryMap1[connectedComponents[i][0]]->lineInfo,1);
							benchmark_hg002_Map.insert(Pair);
						}
					}
				}
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
