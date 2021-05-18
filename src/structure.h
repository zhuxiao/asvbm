#ifndef SRC_STRUCTURE_H_
#define SRC_STRUCTURE_H_

#include <iostream>
#include <cstring>

using namespace std;

#define VAR_UNC				0	// uncertain
#define VAR_INS				1	// insertion
#define VAR_DEL				2	// deletion
#define VAR_DUP				3	// duplication
#define VAR_INV				4	// inversion
#define VAR_TRA				5	// translocation
#define VAR_BND				6	// BND
#define VAR_INV_TRA			7	// inverted translocation
#define VAR_MIX				10	// mixed variation
#define VAR_MNP				11	// MNP
#define VAR_SNV				12	// SNV
#define VAR_CNV				13	// CNV

#define BALANCED_TRA_STR	"BALANCED"
#define UNBALANCED_TRA_STR	"UNBALANCED"

typedef struct {
	string chrname;
	int64_t startPos, endPos;  // 1-based
	string chrname2;
	int64_t startPos2, endPos2;  // 1-based
	bool overlapped;	// default: false
	bool validFlag;		// default: true
	//bool balancedTraFlag;	// TRA, default: false
	//bool isTraBreakpointFlag, isTraBreakpointFlag2;	// TRA, default: false
	bool traOverlappedArr[4];	// four TRA breakpoints, default: false
	size_t sv_type;
	int32_t sv_len;
}SV_item;

typedef struct {
	SV_item *sv_item1, *sv_item2;
	double dif_size, size_ratio, dif_rmse;
}SV_pair;

typedef struct {
	string chrname;
	int64_t bp_loc;
	bool overlappedFlag;
}Breakpoint_t;

typedef struct {
	vector<SV_item*> subset1, subset2;
	vector<SV_item*> *intersect_vec_user, *intersect_vec_benchmark, *private_vec_user, *private_vec_benchmark;
	vector<SV_pair*> *sv_pair_vec;
}overlapWork_opt;

#endif /* SRC_STRUCTURE_H_ */
