#ifndef SRC_CONSTANTS_H_
#define SRC_CONSTANTS_H_

//decoy
#define DECOY_PREFIX				"hs37d"
#define DECOY_PREFIX2				"hs38d"

// maximal valid SVLEN
#define MAX_VALID_SVLEN				99999999	// 99M

//Scoring matrix
#define MATCH_SCORE					2
#define MISMATCH_SCORE				-1
#define GAP_PENALTY				    -1

//Relieving factor of continuous gap
#define RELIEF_FACTOR				0.1

//Sequence consistency condition setting
#define SEQ_CONSISTENCY				0.7f
#define ALLELE_SEQ_CONSISTENCY		0.7f

#define MAX_VALID_REG_THRES			50000
#define EXTEND_SIZE					200  // 100

#define MIN_SIZE_LARGE_SV			100
#define EXTEND_SIZE_LARGE_SV		1000

#define MAXLEN  					15
#define EXTEND_NUM					5
#define WINDOWSIZE					10
#define KMERSIZE					15
#define EDLIBLEN					20000
#define MAX_SEQ_LEN					10000

#define SVLEN_RATIO					0.7f

#define minSizeLittleSV			    100

#define NO_OVERLAP					0
#define WHOLE_REG_OVERLAP			1

//MatchLevel
#define MATCHLEVEL_S				"strict"
#define MATCHLEVEL_L				"loose"

#define ALLELIC_DISTANCE			50

// breakpoints
#define ONE_OVERLAP_ONE_BP			2
#define ONE_OVERLAP_TWO_BP			3
#define ONE_OVERLAP_THREE_BP		4
#define ONE_OVERLAP_FOUR_BP			5
#define TWO_OVERLAP_ONE_BP			6
#define TWO_OVERLAP_TWO_BP			7
#define TWO_OVERLAP_THREE_BP		8
#define TWO_OVERLAP_FOUR_BP			9
#define THREE_OVERLAP_ONE_BP		10
#define THREE_OVERLAP_TWO_BP		11
#define THREE_OVERLAP_THREE_BP		12
#define THREE_OVERLAP_FOUR_BP		13
#define FOUR_OVERLAP_ONE_BP			14
#define FOUR_OVERLAP_TWO_BP			15
#define FOUR_OVERLAP_THREE_BP		16
#define FOUR_OVERLAP_FOUR_BP		17



#endif /* SRC_CONSTANTS_H_ */
