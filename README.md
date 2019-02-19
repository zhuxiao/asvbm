# SV_STAT
A tool for Structural Variation Statistics Evaluation

------------------------------------------------------------------------------
SV_STAT is a tool for Structural Variation (SV) Statistics Evaluation. SV_STAT uses a user-called SV set and a gold standard set as input. It first computes the tranditional metrics, such as the number of true positives (TPs), false positives (TPs), false negatives (FNs), recall, precision and F1 score; and it then compute the difference size between the user-called variants and the one in gold standard set by computing the distance bwtween the centers of the variant regions, and the region size ratio for the two variant regions; and finally, it also computes statistics for different variant types with varies variant size.

Before using SV_STAT, both the user-called SV set and the gold standard set are in bed file format with first five columns are below:
>    chromosome,	start_ref_pos,	end_ref_pos,	SV_type,	SV_len

For translocations, they should be in bedpe format before using SV_STAT, and the first 8 columns are listed as below:
>    chromosome1,	start_ref_pos1,	end_ref_pos1,	chromosome2,	start_ref_pos2,	end_ref_pos2,	SV_type,	SV_len

The SV_type can be TRA or BND, and the SV_len will be 0.

SV_STAT can convert VCF file format to bed or bedpe file format by typing:
>   sv_stat convert -f vcf var.vcf var.bed 


SV_STAT can remove invalid long user-called regions using -m option as they are two long to be a valid variant region, the command can be:
>   sv_stat stat -m 10000 sv_user.bed standard.bed 

------------------------------------------------------------------------------
If you have problems or some suggestions, please contact: xzhu@hrbnu.edu.cn  

---- Enjoy!!! -----

