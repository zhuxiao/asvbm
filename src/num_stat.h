#ifndef SRC_NUM_STAT_H_
#define SRC_NUM_STAT_H_

#include <iostream>
#include <cstring>
#include <vector>
#include <limits.h>
#include <sys/stat.h>
#include <set>
#include <algorithm>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <algorithm>
//#include <chrono>

#include "constants.h"
#include "structure.h"
#include "util.h"
#include "meminfo.h"

using namespace std;

void upperSeq(string &seq);
void SVNumStat(string &user_file, string &benchmark_file, string &ref_file, int32_t max_valid_reg_thres, string &outputPathname, vector<string> &sv_files1);
void SVNumStatOp(string &user_file, string &benchmark_file, string &ref_file, int32_t max_valid_reg_thres, string &dirname);
void computeNumStat(vector<SV_item*> &sv_data1, vector<SV_item*> &sv_data2, string &file_prefix, faidx_t *fai, int Markers);
void computeNumStatFromFile(vector<SV_item*> TPbench_data, vector<SV_item*> TPuser_data, vector<SV_item*> FP_data, vector<SV_item*> FN_data, string &file_prefix, int32_t endpos, int Markers);

void CollectData(float recall, float precision_user, float F1_score_user, double seqcons, vector<float> &Data, size_t num, int Markers);
void CollectData(int LP, int TP_benchmark, int FP, int FN, vector<int> &Data, size_t num, int Markers);
void computeLenStat(vector<SV_item*> &data, string &description_str);

vector<vector<SV_item*>> intersect(vector<SV_item*> &data1, vector<SV_item*> &data2, faidx_t *fai);
vector<vector<SV_item*>> constructSubsetByChr(vector<SV_item*> &user_data, vector<SV_item*> &benchmark_data);
set<string> getChrnames(vector<SV_item*> &dataset);
set<string> getChrUnion(set<string> &chrname_set1, set<string> &chrname_set2);
vector<string> sortChrnames(set<string> &chrname_set);
vector<vector<SV_item*>> constructSubsetByChrOp(vector<SV_item*> &user_data, vector<SV_item*> &benchmark_data, vector<string> &chrname_vec);
void sortSubsets(vector<vector<SV_item*>> &subsets);
void* sortSubsetOp(void *arg);
bool sortFunSameChr(const SV_item *item1, const SV_item *item2);
void checkOrder(vector<vector<SV_item*>> &subsets);
vector<vector<SV_item*>> intersectOp(vector<vector<SV_item*>> &subsets, faidx_t *fai);
void* intersectSubset(void *arg);
vector<SV_item*> getItemsByChr(string &chrname, vector<SV_item*> &dataset);
SV_item* itemdup(SV_item* item);
bool IsSameChrname(string &chrname1, string &chrname2);
vector<size_t> computeOverlapType(SV_item* item1, SV_item* item2);
bool isOverlappedPos(size_t startPos1, size_t endPos1, size_t startPos2, size_t endPos2);
pair<int, vector<SV_item*>> findClosestSubarray(vector<SV_item*> items, int32_t target);
void IsmergeJudge(size_t user_pos, SV_item* user, SV_item* bench, vector<SV_item*> subset1, faidx_t *fai);
//int MINDistance(string seq1, string seq2);
int32_t minDistance(const string &seq1, const string &seq2);
double computeVarseqConsistency(SV_item *item1, SV_item *item2, faidx_t *fai, double SeqIdentity);
void SeqConsNumStat(double consistency);
void extractRefSeq(SV_item* item1, SV_item* item2, string &seq_new1, string &seq_new2, faidx_t *fai);
size_t customHashFunction(const string& kmer, size_t kmerSize);
double MinimizerMethodOp(string& seq1, string& seq2, string& aln_seq1, string& aln_seq2, string& AlignSeq, string& AlignSeq1);
vector<Minimizer> findMinimizers(const string& sequence, size_t windowSize, size_t kmerSize);
void retainCommonHashValuesMinimizers(vector<Minimizer>& containerA, vector<Minimizer>& containerB);
void findSimilarityPosMinimizers(vector<Minimizer>& minimizers, vector<Minimizer>& minimizers1);
void FilterDissimilaPosMinimizers(vector<Minimizer>& minimizers, vector<Minimizer>& minimizers1);
void needleman_wunsch(const string &seq1, const string &seq2, int32_t match_score, int32_t mismatch_score, int32_t gap_penalty, string &seq1_new, string &seq2_new);
void needleman_wunschOp(const string &seq1, const string &seq2, int32_t match_score, int32_t mismatch_score, int32_t gap_penalty, string &seq1_new, string &seq2_new);
int32_t max(int32_t a, int32_t b, int32_t c);
int32_t min(int32_t a, int32_t b, int32_t c);
void LongSequenceSplitAlignment(vector<Minimizer>& minimizers, vector<Minimizer>& minimizers1, int k, string& AlignSeq, string& AlignSeq1, string sequence, string sequence1);
double calculate_consistency(const string& seq1, const string& seq2);

void destroyResultData(vector<vector<SV_item*>> &result);

void SVNumStatTraOp(string &user_file, string &benchmark_file, string &dirname);
void computeNumStatTra(vector<SV_item*> &sv_data1, vector<SV_item*> &sv_data2, string &file_prefix);
void computeOverlapTra(vector<SV_item*> &data1, vector<SV_item*> &data2);
void computeBPNumStatTra(vector<Breakpoint_t*> &bp_vec_user, vector<Breakpoint_t*> &bp_vec_benchmark, string &file_prefix);
Breakpoint_t *allocateBPItem(string &chrname, size_t bp_loc);
void destroyBPData(vector<Breakpoint_t*> &bp_vec);
vector<Breakpoint_t*> constructBPVec(vector<SV_item*> &sv_data);
void outputBP2File(string &filename, vector<Breakpoint_t*> &bp_vec);

#endif /* SRC_NUM_STAT_H_ */
