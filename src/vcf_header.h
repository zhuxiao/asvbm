/*
 * vcf_header.h
 *
 *  Created on: Nov 18, 2024
 *      Author: node
 */

#ifndef SRC_VCF_HEADER_H_
#define SRC_VCF_HEADER_H_
#define VCF_HEADER \
			"##fileformat=VCFv4.2\n" \
			"##fileDate=20241110\n" \
			"##source=ASVCLR 1.4.4,cuteSV 2.1.0,DeBreak,Kled 1.2.9,Sniffles 2.2.2,svdss,SVIM 2.0.0,pbsv 2.9.0\n" \
			"##reference=hs37d5.fa\n" \
			"##phasing=None\n" \
			"##contig=<ID=1,length=249250621>\n" \
			"##contig=<ID=2,length=243199373>\n" \
			"##contig=<ID=3,length=198022430>\n" \
			"##contig=<ID=4,length=191154276>\n" \
			"##contig=<ID=5,length=180915260>\n" \
			"##contig=<ID=6,length=171115067>\n" \
			"##contig=<ID=7,length=159138663>\n" \
			"##contig=<ID=8,length=146364022>\n" \
			"##contig=<ID=9,length=141213431>\n" \
			"##contig=<ID=10,length=135534747>\n" \
			"##contig=<ID=11,length=135006516>\n" \
			"##contig=<ID=12,length=133851895>\n" \
			"##contig=<ID=13,length=115169878>\n" \
			"##contig=<ID=14,length=107349540>\n" \
			"##contig=<ID=15,length=102531392>\n" \
			"##contig=<ID=16,length=90354753>\n" \
			"##contig=<ID=17,length=81195210>\n" \
			"##contig=<ID=18,length=78077248>\n" \
			"##contig=<ID=19,length=59128983>\n" \
			"##contig=<ID=20,length=63025520>\n" \
			"##contig=<ID=21,length=48129895>\n" \
			"##contig=<ID=22,length=51304566>\n" \
			"##contig=<ID=X,length=155270560>\n" \
			"##contig=<ID=Y,length=59373566>\n" \
			"##contig=<ID=MT,length=16569>\n" \
			"##INFO=<ID=TOOLSNUM,Number=1,Type=Integer,Description=\"The number of tools for identifying variants that match this variant\">\n" \
			"##INFO=<ID=TOOLS,Number=1,Type=String,Description=\"Tools for identifying variants that match the variant\">\n" \
			"##INFO=<ID=SS,Number=1,Type=Integer,Description=\"Number of supported signatures\">\n" \
			"##INFO=<ID=ST,Number=1,Type=Integer,Description=\"Number of supported templates\">\n" \
			"##INFO=<ID=SS2,Number=1,Type=Integer,Description=\"Number of supported NGS signatures\">\n" \
			"##INFO=<ID=ST2,Number=1,Type=Integer,Description=\"Number of supported NGS templates\">\n" \
			"##INFO=<ID=LS,Number=1,Type=Integer,Description=\"Length SD ratio score (100-ratio*100) of the cluster.\">\n" \
			"##INFO=<ID=LR,Number=1,Type=String,Description=\"L,R,LR\">\n" \
			"##INFO=<ID=CV,Number=1,Type=Float,Description=\"Nearby coverage for genotyping.\">\n" \
			"##INFO=<ID=CC,Number=1,Type=Float,Description=\"Average coverage of this contig.\">\n" \
			"##INFO=<ID=CN,Number=1,Type=Float,Description=\"Copy number for duplication alleles(normal=1).\">\n" \
			"##INFO=<ID=CR,Number=1,Type=Float,Description=\"Core ratio.\">\n" \
			"##INFO=<ID=M3L,Number=3,Type=Integer,Description=\"(min length, medium length, max length).\">\n" \
			"##INFO=<ID=STDEV_POS,Number=1,Type=Float,Description=\"Standard deviation of structural variation start position\">\n" \
			"##INFO=<ID=STDEV_LEN,Number=1,Type=Float,Description=\"Standard deviation of structural variation length\">\n" \
			"##INFO=<ID=COVERAGE,Number=.,Type=Float,Description=\"Coverages near upstream, start, center, end, downstream of structural variation\">\n" \
			"##INFO=<ID=PSTD,Number=1,Type=Float,Description=\"Position STD, -1 if not calculated.\">\n" \
			"##INFO=<ID=MULTI,Number=0,Type=Flag,Description=\"If the SV is multi-allelic SV\">\n" \
			"##INFO=<ID=MOSAIC,Number=0,Type=Flag,Description=\"Structural variation classified as putative mosaic\">\n" \
			"##INFO=<ID=LARGEINS,Number=0,Type=Flag,Description=\"Large insertion identified from local assembly\">\n" \
			"##INFO=<ID=START2,Number=1,Type=Integer,Description=\"SV start position on the second haplotype of multi-allelic SV\">\n" \
			"##INFO=<ID=END2,Number=1,Type=Integer,Description=\"SV end position on the second haplotype of multi-allelic SV\">\n" \
			"##INFO=<ID=SVANN,Number=.,Type=String,Description=\"Repeat annotation of structural variant\">\n" \
			"##INFO=<ID=SVLEN2,Number=1,Type=Integer,Description=\"SV length on the second haplotype of multi-allelic SV\">\n" \
			"##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Type of approach used to detect SV\">\n" \
			"##INFO=<ID=MAPQ,Number=1,Type=Integer,Description=\"Mean mapping quality of supporting reads\">\n" \
			"##INFO=<ID=SUPPREAD,Number=1,Type=Integer,Description=\"Number of supporting reads\">\n" \
			"##INFO=<ID=SUPPORT_INLINE,Number=1,Type=Integer,Description=\"Number of reads supporting an INS/DEL SV (non-split events only)\">\n" \
			"##INFO=<ID=SUPPORT_LONG,Number=1,Type=Integer,Description=\"Number of soft-clipped reads putatively supporting the long insertion SV\">\n" \
			"##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n" \
			"##INFO=<ID=CILEN,Number=2,Type=Integer,Description=\"Confidence interval around inserted/deleted material between breakends\">\n" \
			"##INFO=<ID=RE,Number=1,Type=Integer,Description=\"Number of read support this record\">\n" \
			"##INFO=<ID=STRAND,Number=A,Type=String,Description=\"Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)\">\n" \
			"##INFO=<ID=RNAMES,Number=.,Type=String,Description=\"Supporting read names of SVs (comma separated)\">\n" \
			"##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">\n" \
			"##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variant\">\n" \
			"##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variant\">\n" \
			"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant described in this record\">\n" \
			"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n" \
			"##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n" \
			"##INFO=<ID=DUPNUM,Number=1,Type=Integer,Description=\"Copy number of DUP\">\n" \
			"##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">\n" \
			"##INFO=<ID=MATEDIST,Number=1,Type=Integer,Description=\"Distance to the mate breakend for mates on the same contig\">\n" \
			"##INFO=<ID=LDISCOV,Number=1,Type=String,Description=\"Variant discover level: ALN for variants are discovered from consensus alignments, RES-ALN for variants are discovered from rescued consensus alignments around variant regions, READS for variants are discovered from reads alignments directly\">\n" \
			"##INFO=<ID=VARTYPE,Number=A,Type=String,Description=\"Variant class\">\n" \
			"##INFO=<ID=WEIGHT,Number=1,Type=Integer,Description=\"Number of alignments supporting this record\">\n" \
			"##INFO=<ID=COV,Number=1,Type=Integer,Description=\"Total number of alignments covering this locus\">\n" \
			"##INFO=<ID=COV0,Number=1,Type=Integer,Description=\"Total number of alignments covering this locus (no HP)\">\n" \
			"##INFO=<ID=COV1,Number=1,Type=Integer,Description=\"Total number of alignments covering this locus (HP=1)\">\n" \
			"##INFO=<ID=COV2,Number=1,Type=Integer,Description=\"Total number of alignments covering this locus (HP=2)\">\n" \
			"##INFO=<ID=AS,Number=1,Type=Integer,Description=\"Alignment score\">\n" \
			"##INFO=<ID=NV,Number=1,Type=Integer,Description=\"Number of variations on same consensus\">\n" \
			"##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n" \
			"##INFO=<ID=CIGAR,Number=A,Type=String,Description=\"CIGAR of consensus\">\n" \
			"##INFO=<ID=READS,Number=.,Type=String,Description=\"Reads identifiers supporting the call\">\n" \
			"##INFO=<ID=RVEC,Number=.,Type=String,Description=\"Reads vector used by genotyper\">\n" \
			"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n" \
			"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n" \
			"##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Allele count, summed up over all samples\">\n" \
			"##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description=\"List of read support for all samples\">\n" \
			"##INFO=<ID=CONSENSUS_SUPPORT,Number=1,Type=Integer,Description=\"Number of reads that support the generated insertion (INS) consensus sequence\">\n" \
			"##INFO=<ID=NM,Number=.,Type=Float,Description=\"Mean number of query alignment length adjusted mismatches of supporting reads\">\n" \
			"##INFO=<ID=PHASE,Number=.,Type=String,Description=\"Phasing information derived from supporting reads, represented as list of: HAPLOTYPE,PHASESET,HAPLOTYPE_SUPPORT,PHASESET_SUPPORT,HAPLOTYPE_FILTER,PHASESET_FILTER\">\n" \
			"##INFO=<ID=CUTPASTE,Number=0,Type=Flag,Description=\"Genomic origin of interspersed duplication seems to be deleted\">\n" \
			"##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"Number of reads supporting this variant\">\n" \
			"##INFO=<ID=STD_SPAN,Number=1,Type=Float,Description=\"Standard deviation in span of merged SV signatures\">\n" \
			"##INFO=<ID=STD_POS,Number=1,Type=Float,Description=\"Standard deviation in position of merged SV signatures\">\n" \
			"##INFO=<ID=STD_POS1,Number=1,Type=Float,Description=\"Standard deviation of breakend 1 position\">\n" \
			"##INFO=<ID=STD_POS2,Number=1,Type=Float,Description=\"Standard deviation of breakend 2 position\">\n" \
			"##INFO=<ID=NotFullySpanned,Number=0,Type=Flag,Description=\"Duplication variant does not have any fully spanning reads\">\n" \
			"##INFO=<ID=ClusterIDs,Number=1,Type=String,Description=\"IDs of SVs that cluster with this SV\">\n" \
			"##INFO=<ID=NumClusterSVs,Number=1,Type=Integer,Description=\"Total number of SV calls in this cluster\">\n" \
			"##INFO=<ID=ExactMatchIDs,Number=1,Type=String,Description=\"IDs of SVs that are exactly the same call as this SV\">\n" \
			"##INFO=<ID=NumExactMatchSVs,Number=1,Type=Integer,Description=\"Total number of SVs in this exact cluster\">\n" \
			"##INFO=<ID=ClusterMaxShiftDist,Number=1,Type=Float,Description=\"Maximum relative shift distance between two SVs in this cluster\">\n" \
			"##INFO=<ID=ClusterMaxSizeDiff,Number=1,Type=Float,Description=\"Maximum relative size difference between two SVs in this cluster\">\n" \
			"##INFO=<ID=ClusterMaxEditDist,Number=1,Type=Float,Description=\"Maximum relative edit distance between two SVs in this cluster\">\n" \
			"##INFO=<ID=PBcalls,Number=1,Type=Integer,Description=\"Number of PacBio calls in this cluster\">\n" \
			"##INFO=<ID=Illcalls,Number=1,Type=Integer,Description=\"Number of Illumina calls in this cluster\">\n" \
			"##INFO=<ID=TenXcalls,Number=1,Type=Integer,Description=\"Number of 10X Genomics calls in this cluster\">\n" \
			"##INFO=<ID=CGcalls,Number=1,Type=Integer,Description=\"Number of Complete Genomics calls in this cluster\">\n" \
			"##INFO=<ID=PBexactcalls,Number=1,Type=Integer,Description=\"Number of PacBio calls exactly matching the call output for this cluster\">\n" \
			"##INFO=<ID=Illexactcalls,Number=1,Type=Integer,Description=\"Number of Illumina calls exactly matching the call output for this cluster\">\n" \
			"##INFO=<ID=TenXexactcalls,Number=1,Type=Integer,Description=\"Number of 10X Genomics calls exactly matching the call output for this cluster\">\n" \
			"##INFO=<ID=CGexactcalls,Number=1,Type=Integer,Description=\"Number of Complete Genomics calls exactly matching the call output for this cluster\">\n" \
			"##INFO=<ID=HG2count,Number=1,Type=Integer,Description=\"Number of calls discovered in HG002 in this cluster\">\n" \
			"##INFO=<ID=HG3count,Number=1,Type=Integer,Description=\"Number of calls discovered in HG003 in this cluster\">\n" \
			"##INFO=<ID=HG4count,Number=1,Type=Integer,Description=\"Number of calls discovered in HG004 in this cluster\">\n" \
			"##INFO=<ID=NumTechs,Number=1,Type=Integer,Description=\"Number of technologies from which calls were discovered in this cluster\">\n" \
			"##INFO=<ID=NumTechsExact,Number=1,Type=Integer,Description=\"Number of technologies from which calls were discovered that exactly match the call output for this cluster\">\n" \
			"##INFO=<ID=DistBack,Number=1,Type=Integer,Description=\"Distance to the closest non-matching variant before this variant\">\n" \
			"##INFO=<ID=DistForward,Number=1,Type=Integer,Description=\"Distance to the closest non-matching variant after this variant\">\n" \
			"##INFO=<ID=DistMin,Number=1,Type=Integer,Description=\"Distance to the closest non-matching variant in either direction\">\n" \
			"##INFO=<ID=DistMinlt1000,Number=1,Type=String,Description=\"TRUE if Distance to the closest non-matching variant in either direction is less than 1000bp, suggesting possible complex or compound heterozygous variant\">\n" \
			"##INFO=<ID=MultiTech,Number=1,Type=String,Description=\"TRUE if callsets from more than one technology are in this cluster, i.e., NumTechs>1\">\n" \
			"##INFO=<ID=MultiTechExact,Number=1,Type=String,Description=\"TRUE if callsets from more than one technology exactly matches the call output for this cluster, i.e., NumTechsExact>1\">\n" \
			"##INFO=<ID=sizecat,Number=1,Type=String,Description=\"Size category in bp: 20to49, 50to99, 100to299, 300to999, or gt1000\">\n" \
			"##INFO=<ID=DistPASSHG2gt49Minlt1000,Number=1,Type=String,Description=\"TRUE if Distance to the closest non-matching PASS variant >49bp in HG002 in either direction is less than 1000bp, suggesting possible complex or compound heterozygous variant or inaccurate call\">\n" \
			"##INFO=<ID=DistPASSMinlt1000,Number=1,Type=String,Description=\"TRUE if Distance to the closest non-matching PASS variant in either direction is less than 1000bp, suggesting possible complex or compound heterozygous variant or inaccurate call\">\n" \
			"##INFO=<ID=MendelianError,Number=1,Type=String,Description=\"TRUE if all individuals have a consensus GT and they are not consistent with Mendelian inheritance\">\n" \
			"##INFO=<ID=HG003_GT,Number=1,Type=String,Description=\"Consensus GT for HG003/Father\">\n" \
			"##INFO=<ID=HG004_GT,Number=1,Type=String,Description=\"Consensus GT for HG004/Mother\">\n" \
			"##INFO=<ID=BREAKSIMLENGTH,Number=1,Type=Integer,Description=\"Length of alignable similarity at event breakpoints as determined by the aligner\">\n" \
			"##INFO=<ID=REFWIDENED,Number=1,Type=String,Description=\"Widened boundaries of the event in the reference allele\">\n" \
			"##INFO=<ID=REPTYPE,Number=1,Type=String,Description=\"Type of SV, with designation of uniqueness of new or deleted sequence:SIMPLEDEL=Deletion of at least some unique sequence, SIMPLEINS=Insertion of at least some unique sequence, CONTRAC=Contraction, or deletion of sequence entirely similar to remaining sequence, DUP=Duplication, or insertion of sequence entirely similar to pre-existing sequence, INV=Inversion, SUBSINS=Insertion of new sequence with alteration of some pre-existing sequence, SUBSDEL=Deletion of sequence with alteration of some remaining sequence\">\n" \
			"##INFO=<ID=TRall,Number=1,Type=String,Description=\"TRUE if at least 20% of the REF bases are tandem repeats of any length\">\n" \
			"##INFO=<ID=TRgt100,Number=1,Type=String,Description=\"TRUE if at least 20% of the REF bases are tandem repeats at least 100bp long\">\n" \
			"##INFO=<ID=TRgt10k,Number=1,Type=String,Description=\"TRUE if at least 20% of the REF bases are tandem repeats at least 10kbp long\">\n" \
			"##INFO=<ID=segdup,Number=1,Type=String,Description=\"TRUE if at least 20% of the REF bases are segmental duplications at least 10kbp long\">\n" \
			"##INFO=<ID=BREAKSIMLENGTH,Number=1,Type=Integer,Description=\"Length of alignable similarity at event breakpoints as determined by the aligner\">\n" \
			"##INFO=<ID=REFWIDENED,Number=1,Type=String,Description=\"Widened boundaries of the event in the reference allele\">\n" \
			"##INFO=<ID=REPTYPE,Number=1,Type=String,Description=\"Type of SV, with designation of uniqueness of new or deleted sequence:SIMPLEDEL=Deletion of at least some unique sequence, SIMPLEINS=Insertion of at least some unique sequence, CONTRAC=Contraction, or deletion of sequence entirely similar to remaining sequence, DUP=Duplication, or insertion of sequence entirely similar to pre-existing sequence, INV=Inversion, SUBSINS=Insertion of new sequence with alteration of some pre-existing sequence, SUBSDEL=Deletion of sequence with alteration of some remaining sequence\">\n" \
			"##FILTER=<ID=hom_ref,Description=\"Genotype is homozygous reference\">\n" \
			"##FILTER=<ID=NoConsensusGT,Description=\"No individual had genotypes from svviz agree across all datasets with confident genotypes\">\n" \
			"##FILTER=<ID=LongReadHomRef,Description=\"Long reads supported homozygous reference for all individuals\">\n" \
			"##FILTER=<ID=ClusteredCalls,Description=\"Supported SVs that are within 1kb of a different supported SV\">\n" \
			"##FILTER=<ID=lt50bp,Description=\"Supported variant but smaller than 50bp\">\n" \
			"##FILTER=<ID=not_fully_covered,Description=\"Tandem duplication is not fully covered by a single read\">\n" \
			"##FILTER=<ID=q5,Description=\"Quality below 5\">\n" \
			"##FILTER=<ID=PASS,Description=\"All filters passed\">\n" \
			"##FILTER=<ID=GT,Description=\"Genotype filter\">\n" \
			"##FILTER=<ID=SUPPORT_MIN,Description=\"Minimum read support filter\">\n" \
			"##FILTER=<ID=STDEV_POS,Description=\"SV Breakpoint standard deviation filter\">\n" \
			"##FILTER=<ID=STDEV_LEN,Description=\"SV length standard deviation filter\">\n" \
			"##FILTER=<ID=COV_MIN,Description=\"Minimum coverage filter\">\n" \
			"##FILTER=<ID=COV_CHANGE,Description=\"Coverage change filter\">\n" \
			"##FILTER=<ID=COV_CHANGE_FRAC,Description=\"Coverage fractional change filter\">\n" \
			"##FILTER=<ID=MOSAIC_AF,Description=\"Mosaic maximum allele frequency filter\">\n" \
			"##FILTER=<ID=ALN_NM,Description=\"Length adjusted mismatch filter\">\n" \
			"##FILTER=<ID=STRAND,Description=\"Strand support filter\">\n" \
			"##FILTER=<ID=SVLEN_MIN,Description=\"SV length filter\">\n" \
			"##FILTER=<ID=Decoy,Description=\"Variant involves a decoy sequence\">\n" \
			"##FILTER=<ID=NearReferenceGap,Description=\"Variant is near (< 1000 bp) from a gap (run of >= 50 Ns) in the reference assembly\">\n" \
			"##FILTER=<ID=NearContigEnd,Description=\"Variant is near (< 1000 bp) from the end of a contig\">\n" \
			"##FILTER=<ID=InsufficientStrandEvidence,Description=\"Variant has insufficient number of reads per strand (< 0).\">\n" \
			"##FORMAT=<ID=ID,Number=1,Type=String,Description=\"Individual sample SV ID for multi-sample output\">\n" \
			"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" \
			"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read depth per allele\">\n" \
			"##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# High-quality reference reads\">\n" \
			"##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# High-quality variant reads\">\n" \
			"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"# Phred-scaled genotype likelihoods rounded to the closest integer\">\n" \
			"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"# Genotype quality\">\n" \
			"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth at this position for this sample\">\n" \
			"##FORMAT=<ID=SAC,Number=.,Type=Integer,Description=\"Number of reads on the forward and reverse strand supporting each allele including reference\">\n" \
			"##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number genotype for imprecise events\">\n" \
			"##FORMAT=<ID=GTcons1,Number=1,Type=String,Description=\"Consensus Genotype using the GT from svviz2 rather than ref and alt allele counts, which is sometimes inaccurate for large variants\">\n" \
			"##FORMAT=<ID=PB_GT,Number=1,Type=String,Description=\"Genotype predicted by svviz from PacBio\">\n" \
			"##FORMAT=<ID=PB_REF,Number=1,Type=Integer,Description=\"Number of PacBio reads supporting the REF allele as predicted by svviz\">\n" \
			"##FORMAT=<ID=PB_ALT,Number=1,Type=Integer,Description=\"Number of PacBio reads supporting the ALT allele as predicted by svviz\">\n" \
			"##FORMAT=<ID=PBHP_GT,Number=1,Type=String,Description=\"Genotype predicted by svviz from PacBio haplotype-separated by 10X\">\n" \
			"##FORMAT=<ID=PB_REF_HP1,Number=1,Type=Integer,Description=\"Number of PacBio reads on haplotype 1 supporting the REF allele as predicted by svviz\">\n" \
			"##FORMAT=<ID=PB_ALT_HP1,Number=1,Type=Integer,Description=\"Number of PacBio reads on haplotype 1 supporting the ALT allele as predicted by svviz\">\n" \
			"##FORMAT=<ID=PB_REF_HP2,Number=1,Type=Integer,Description=\"Number of PacBio reads on haplotype 2 supporting the REF allele as predicted by svviz\">\n" \
			"##FORMAT=<ID=PB_ALT_HP2,Number=1,Type=Integer,Description=\"Number of PacBio reads on haplotype 2 supporting the ALT allele as predicted by svviz\">\n" \
			"##FORMAT=<ID=TenX_GT,Number=1,Type=String,Description=\"Genotype predicted by svviz from 10X by combining the GT's from each haplotype\">\n" \
			"##FORMAT=<ID=TenX_REF_HP1,Number=1,Type=Integer,Description=\"Number of 10X reads on haplotype 1 supporting the REF allele as predicted by svviz\">\n" \
			"##FORMAT=<ID=TenX_ALT_HP1,Number=1,Type=Integer,Description=\"Number of 10X reads on haplotype 1 supporting the ALT allele as predicted by svviz\">\n" \
			"##FORMAT=<ID=TenX_REF_HP2,Number=1,Type=Integer,Description=\"Number of 10X reads on haplotype 2 supporting the REF allele as predicted by svviz\">\n" \
			"##FORMAT=<ID=TenX_ALT_HP2,Number=1,Type=Integer,Description=\"Number of 10X reads on haplotype 2 supporting the ALT allele as predicted by svviz\">\n" \
			"##FORMAT=<ID=ILL250bp_GT,Number=1,Type=String,Description=\"Genotype predicted by svviz from Illumina 250bp reads\">\n" \
			"##FORMAT=<ID=ILL250bp_REF,Number=1,Type=Integer,Description=\"Number of Illumina 250bp reads supporting the REF allele as predicted by svviz\">\n" \
			"##FORMAT=<ID=ILL250bp_ALT,Number=1,Type=Integer,Description=\"Number of Illumina 250bp reads supporting the ALT allele as predicted by svviz\">\n" \
			"##FORMAT=<ID=ILLMP_GT,Number=1,Type=String,Description=\"Genotype predicted by svviz from Illumina mate-pair reads\">\n" \
			"##FORMAT=<ID=ILLMP_REF,Number=1,Type=Integer,Description=\"Number of Illumina mate-pair reads supporting the REF allele as predicted by svviz\">\n" \
			"##FORMAT=<ID=ILLMP_ALT,Number=1,Type=Integer,Description=\"Number of Illumina mate-pair reads supporting the ALT allele as predicted by svviz\">\n" \
			"##FORMAT=<ID=BNG_LEN_DEL,Number=1,Type=Integer,Description=\"Length of a deletion predicted by BioNano in a region overlapping this variant\">\n" \
			"##FORMAT=<ID=BNG_LEN_INS,Number=1,Type=Integer,Description=\"Length of an insertion predicted by BioNano in a region overlapping this variant\">\n" \
			"##FORMAT=<ID=nabsys_svm,Number=1,Type=Float,Description=\"Nabsys SVM score for this variant if it was evaluated\">\n" \
			"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample" \



#endif /* SRC_VCF_HEADER_H_ */
