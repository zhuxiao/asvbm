#include "paras.h"
#include "global.h"

// get the program version
string getProgramVersion(const string &cmd_str){
	FILE *stream;
	char tmp[256], info[256] = {0};
	string pg_version_str = "";

	sprintf(tmp, "%s", cmd_str.c_str());
	stream = popen(tmp, "r");
	if(fread(info, 1, sizeof(info), stream)>0){
		pg_version_str = info;
		if(pg_version_str.at(pg_version_str.size()-1)=='\n'){
			pg_version_str.erase(pg_version_str.size()-1); // Remove the newline character
		}
	}
	pclose(stream);

	return pg_version_str;
}
bool check_bam_average_read_length() {
    const int REGION_LENGTH = 1000;
    const int MAX_ATTEMPTS = 1000;
    hts_set_log_level(HTS_LOG_OFF);
    samFile* in = sam_open(BamFilePath.c_str(), "r");
    if (!in) return false;

    bam_hdr_t* header = sam_hdr_read(in);
    if (!header) {
        sam_close(in);
        return false;
    }

    hts_idx_t* idx = sam_index_load(in, BamFilePath.c_str());
    if (!idx) {
        bam_hdr_destroy(header);
        sam_close(in);
        return false;
    }
    // Collect valid reference sequences (length >= REGION_LENGTH)
    vector<int> valid_tids;
    for (int i = 0; i < header->n_targets; ++i) {
        if (header->target_len[i] >= REGION_LENGTH) {
            valid_tids.push_back(i);
        }
    }
    if (valid_tids.empty()) {
        hts_idx_destroy(idx);
        bam_hdr_destroy(header);
        sam_close(in);
        return false;
    }

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> tid_dist(0, valid_tids.size() - 1);
    bam1_t* read = bam_init1();
    vector<double> averages;

    for (int i = 0; i < 3; ++i) {
        bool found = false;
        int attempts = 0;

        while (!found && attempts++ < MAX_ATTEMPTS) {
            int tid = valid_tids[tid_dist(gen)];
            int contig_len = header->target_len[tid];
            uniform_int_distribution<> pos_dist(0, contig_len - REGION_LENGTH);
            int start = pos_dist(gen);
            int end = start + REGION_LENGTH;

            hts_itr_t* iter = sam_itr_queryi(idx, tid, start, end);
            if (!iter) continue;

            uint64_t total_len = 0;
            int count = 0;
            int ret;
            while ((ret = sam_itr_next(in, iter, read)) >= 0) {
                count++;
                total_len += read->core.l_qseq;
            }

            hts_itr_destroy(iter);

            if (count > 0) {
                averages.push_back(static_cast<double>(total_len) / count);
                found = true;
            }
        }

        if (!found) {
            bam_destroy1(read);
            hts_idx_destroy(idx);
            bam_hdr_destroy(header);
            sam_close(in);
            return false;
        }
    }

    bam_destroy1(read);
    hts_idx_destroy(idx);
    bam_hdr_destroy(header);
    sam_close(in);

    for (double avg : averages) {
        if (avg >= 2000.0) {
            return true;
        }
    }
    return false;
}
//SVregion
vector<string> generateSVregion(){
	vector<string> scenarios;

	if(size_div_vec.empty()){
		scenarios.push_back(">0");
		return scenarios;
	}

	scenarios.push_back("1-" + to_string(size_div_vec[0]) + "bp");

	for(size_t i = 1; i < size_div_vec.size(); ++i){
		scenarios.push_back(to_string(size_div_vec[i - 1] + 1) + "-" + to_string(size_div_vec[i]) + "bp");
	}

	scenarios.push_back(">=" + to_string(size_div_vec.back() + 1) + "bp");

	return scenarios;
}
//show version
void show_version(){
	cout << PROG_VERSION << endl;
}

// show usage
void showUsage(){

	cout << "Program: " << PROG_NAME << " (" << PROG_DESC << ")" << endl;
	cout << "Version: " << PROG_VERSION << endl << endl;

	cout << "Usage:  sv_stat  <command> [options]" << endl << endl;

	cout << "Commands:" << endl;
	cout << "     convert      convert SV results to 5-column BED (or 8-column BEDPE) file format" << endl;
	cout << "     stat         compute the SV statistics" << endl;
}

// show usage for convert command

void showUsageConvert(){

	cout << "Program: " << PROG_NAME << " (" << PROG_DESC << ")" << endl;
	cout << "Version: " << PROG_VERSION << endl << endl;

	cout << "Usage:  sv_stat convert [options] <Ref> <infile> <outfile>" << endl << endl;

	cout << "Options:" << endl;
	cout << "     -f STR       SV input file format (required):" << endl;
	cout << "                  bed: BED/BEDPE format" << endl;
	cout << "                  vcf: VCF format" << endl;
	cout << "                  csv: CSV format" << endl;
	cout << "     -o FILE      output directory: [output]" << endl;
	cout << "     -r INT       remove mate variant items [1]: 1 for yes, 0 for no" << endl;
	cout << "     -m FILE      duplicated mate variant items file:" << endl;
	cout << "                  [" << mateItemFilename << "]" << endl;
	cout << "     -s INT       remove snv items [1]: 1 for yes, 0 for no" << endl;
	cout << "     -S FILE      snv items file: [" << snvFilename << "]" << endl;
	cout << "     -h           show this help message and exit" << endl;
}
void showUsageCreate(){

	cout << "Program: " << PROG_NAME << " (" << PROG_DESC << ")" << endl;
	cout << "Version: " << PROG_VERSION << endl << endl;

	cout << "Usage:  asvbm create [options] <USER_FILE> [<USER_FILE1>...] <BENCH_FILE> <REF_FILE>" << endl << endl;

	cout << "Description:" << endl;
	cout << "   USER_FILE      User called SV result file." << endl;
	cout << "   BENCH_FILE     Benchmark SV file." << endl;
	cout << "   REF_FILE       Reference file." << endl << endl;

	cout << "Options:" << endl;
	cout << "   -M INT         valid maximal region size for statistics: [" << MAX_VALID_REG_THRES << "]" << endl;
	cout << "                  0 is for all variant size are valid, and while positive" << endl;
	cout << "                  values are for the valid maximal region size, then longer" << endl;
	cout << "                  regions are omitted and saved to the file specified with '-l' option" << endl;
	cout << "   -m INT         valid minimal region size for statistics: [" << MIN_SVLEN << "]" << endl;
	cout << "                  values are for the valid minimal region size, then shorter" << endl;
	cout << "                  regions are omitted and saved to the file specified with '-l' option" << endl;
	cout << "   -S             enable the strict type match mode which is disabled by default." << endl;
	cout << "                  There are two variation type match modes:" << endl;
	cout << "                  " << MATCHLEVEL_L << ": allow type match between DUP and INS, which takes effect by '-S' option" << endl;
	cout << "                  " << MATCHLEVEL_S << ": strict type match which is disabled by default" << endl;
	cout << "                  The default enabled match mode is 'loose' to allow the type match between DUP and INS." << endl;
	cout << "   -C STR         Chromosomes to be processed: [null]" << endl;
	cout << "                  no decoy indicates not specifying the chromosome set for benchmarking." << endl;
	cout << "                  This parameter is used to specify the chromosomes to be benchmarked." << endl;
	cout << "                  Chromosome names should match the format within the VCF file. " <<endl;
	cout << "                  Chromosome names are separated by ';'. Example: -C \"1;2;3\" " << endl;
	cout << "   -s INT         overlap extend size: [" << EXTEND_SIZE << "]" << endl;
	cout << "   -i FLOAT       minimal sequence similarly for variant match: [" << SEQ_CONSISTENCY << "]" << endl;
	cout << "   -a FLOAT       minimal sequence similarly for allelic variants match: [" << ALLELE_SEQ_CONSISTENCY << "]" << endl;
	cout << "   -p FLOAT       minimal percent size ratio for variant match: [" << SVLEN_RATIO << "]" << endl;
	cout << "   -t INT         number of threads [0]. 0 for the maximal number of threads" << endl;
	cout << "                  in machine" << endl;
	cout << "   -T STR         Tool names [null]." << endl;
	cout << "                  This parameter is used for comparing multiple datasets. The number" << endl;
	cout << "                  of inputs should be consistent with the data set. Tool names are " <<endl;
	cout << "                  separated by ';'. Example: -T \"tool1;tool2;tool3\" " << endl;
	cout << "   -B STR         BAM file [null]." << endl;
	cout << "                  This parameter is used for local joint validation analysis." << endl;
	cout << "                  Recommend for achieving more accurate local variant joint analysis." <<endl;
	cout << "                  BAM file that supports long-read sequencing data as input." <<endl;
	cout << "   -R STR         Custom region sizes [null]." << endl;
	cout << "                  This parameter allows you to specify custom region sizes for analysis." << endl;
	cout << "                  Custom region sizes are separated by ';'. Example: -R \"100;500;1000\" " <<endl;
	cout << "   -o FILE        output directory: [" << outputPathname << "]" << endl;
	cout << "   -l FILE        file name of long SV regions: [" << longSVFilename << "]" << endl;
	cout << "                  file name of short SV regions: [" << shortSVFilename << "]" << endl;
	cout << "   -r FILE        file name of benchmarking results to report: [" << htmlFilename << "]" << endl;
	cout << "                  Ensure that the filename extension is '.html'." << endl;
	cout << "   -v,--version   show version information" << endl;
	cout << "   -h,--help      show this help message and exit" << endl << endl;

	cout << "Example:" << endl;
	cout << "   # run the benchmarking on the user-called set (method) for a single sample to allow match between DUPs as INSs" << endl;
	cout << "   $ asvbm create -T method user_sv.vcf benchmark_sv.vcf ref.fa" << endl << endl;

	cout << "   # run the benchmarking on the user-called set (method) for a single sample to perform the strict type matching by '-S' option" << endl;
	cout << "   $ asvbm create -T method -S user_sv.vcf benchmark_sv.vcf ref.fa" << endl << endl;

	cout << "   # run the benchmarking on the user-called sets (tool1, tool2 and tool3) for multiple user callsets" << endl;
	cout << "   $ asvbm create -T \"tool1;tool2;tool3\" user_sv1.vcf user_sv2.vcf user_sv3.vcf benchmark_sv.vcf ref.fa" << endl;
}
// show usage for stat command
void showUsageStat(){

	cout << "Program: " << PROG_NAME << " (" << PROG_DESC << ")" << endl;
	cout << "Version: " << PROG_VERSION << endl << endl;

	cout << "Usage:  asvbm stat [options] <USER_FILE> [<USER_FILE1>...] <BENCH_FILE> <REF_FILE>" << endl << endl;

	cout << "Description:" << endl;
	cout << "   USER_FILE      User called SV result file." << endl;
	cout << "   BENCH_FILE     Benchmark SV file." << endl;
	cout << "   REF_FILE       Reference file." << endl << endl;

	cout << "Options:" << endl;
	cout << "   -M INT         valid maximal region size for statistics: [" << MAX_VALID_REG_THRES << "]" << endl;
	cout << "                  0 is for all variant size are valid, and while positive" << endl;
	cout << "                  values are for the valid maximal region size, then longer" << endl;
	cout << "                  regions are omitted and saved to the file specified with '-l' option" << endl;
	cout << "   -m INT         valid minimal region size for statistics: [" << MIN_SVLEN << "]" << endl;
	cout << "                  values are for the valid minimal region size, then shorter" << endl;
	cout << "                  regions are omitted and saved to the file specified with '-l' option" << endl;
	cout << "   -S             enable the strict type match mode which is disabled by default." << endl;
	cout << "                  There are two variation type match modes:" << endl;
	cout << "                  " << MATCHLEVEL_L << ": allow type match between DUP and INS, which takes effect by '-S' option" << endl;
	cout << "                  " << MATCHLEVEL_S << ": strict type match which is disabled by default" << endl;
	cout << "                  The default enabled match mode is 'loose' to allow the type match between DUP and INS." << endl;
	cout << "   -C STR         Chromosomes to be processed: [null]" << endl;
	cout << "                  no decoy indicates not specifying the chromosome set for benchmarking." << endl;
	cout << "                  This parameter is used to specify the chromosomes to be benchmarked." << endl;
	cout << "                  Chromosome names should match the format within the VCF file. " <<endl;
	cout << "                  Chromosome names are separated by ';'. Example: -C \"1;2;3\" " << endl;
	cout << "   -s INT         overlap extend size: [" << EXTEND_SIZE << "]" << endl;
	cout << "   -i FLOAT       minimal sequence similarly for variant match: [" << SEQ_CONSISTENCY << "]" << endl;
	cout << "   -a FLOAT       minimal sequence similarly for allelic variants match: [" << ALLELE_SEQ_CONSISTENCY << "]" << endl;
	cout << "   -p FLOAT       minimal percent size ratio for variant match: [" << SVLEN_RATIO << "]" << endl;
	cout << "   -t INT         number of threads [0]. 0 for the maximal number of threads" << endl;
	cout << "                  in machine" << endl;
	cout << "   -T STR         Tool names [null]." << endl;
	cout << "                  This parameter is used for comparing multiple datasets. The number" << endl;
	cout << "                  of inputs should be consistent with the data set. Tool names are " <<endl;
	cout << "                  separated by ';'. Example: -T \"tool1;tool2;tool3\" " << endl;
	cout << "   -B STR         BAM file [null]." << endl;
	cout << "                  This parameter is used for local joint validation analysis." << endl;
	cout << "                  Recommend for achieving more accurate local variant joint analysis." <<endl;
	cout << "                  BAM file that supports long-read sequencing data as input." <<endl;
	cout << "   -R STR         Custom region sizes [null]." << endl;
	cout << "                  This parameter allows you to specify custom region sizes for analysis." << endl;
	cout << "                  Custom region sizes are separated by ';'. Example: -R \"100;500;1000\" " <<endl;
	cout << "   -b,--bench_refine"<< endl;
	cout << "                  perform the refinement of benchmark set along with the benchmarking process" << endl;
	cout << "   -o FILE        output directory: [" << outputPathname << "]" << endl;
	cout << "   -l FILE        file name of long SV regions: [" << longSVFilename << "]" << endl;
	cout << "                  file name of short SV regions: [" << shortSVFilename << "]" << endl;
	cout << "   -r FILE        file name of benchmarking results to report: [" << htmlFilename << "]" << endl;
	cout << "                  Ensure that the filename extension is '.html'." << endl;
	cout << "   -v,--version   show version information" << endl;
	cout << "   -h,--help      show this help message and exit" << endl << endl;

	cout << "Example:" << endl;
	cout << "   # run the benchmarking on the user-called set (method) for a single sample to allow match between DUPs as INSs" << endl;
	cout << "   $ asvbm stat -T method user_sv.vcf benchmark_sv.vcf ref.fa" << endl << endl;

	cout << "   # run the benchmarking on the user-called set (method) for a single sample to perform the strict type matching by '-S' option" << endl;
	cout << "   $ asvbm stat -b -T method -S user_sv.vcf benchmark_sv.vcf ref.fa" << endl << endl;

	cout << "   # run the benchmarking on the user-called sets (tool1, tool2 and tool3) for multiple user callsets" << endl;
	cout << "   $ asvbm stat --bench_refine -T \"tool1;tool2;tool3\" user_sv1.vcf user_sv2.vcf user_sv3.vcf benchmark_sv.vcf ref.fa" << endl;
}

void show(){

	cout << "Program: " << PROG_NAME << " (" << PROG_DESC << ")" << endl;
	cout << "Version: " << PROG_VERSION << endl << endl;

	cout << "Usage:  asvbm <command> [options] <USER_FILE> [<USER_FILE1>...] <BENCH_FILE> <REF_FILE>" << endl << endl;

	cout << "Description:" << endl;
	cout << "   USER_FILE    User called SV result file." << endl;
	cout << "   BENCH_FILE   Benchmark SV file." << endl;
	cout << "   REF_FILE     Reference file." << endl << endl;

	cout << "Commands:" << endl;
	cout << "   stat         revise the existing benchmark set" << endl;
	cout << "   create       create a new benchmark set" << endl<< endl;

//	cout << "Options:" << endl;
//	cout << "   -m INT    valid maximal region size for statistics: [" << MAX_VALID_REG_THRES << "]" << endl;
//	cout << "             0 is for all variant size are valid, and while positive" << endl;
//	cout << "             values are for the valid maximal region size, then longer" << endl;
//	cout << "             regions are omitted and saved to the file specified with '-l' option" << endl;
//	cout << "   -S        enable the strict type match mode which is disabled by default." << endl;
//	cout << "             There are two variation type match modes:" << endl;
//	cout << "             " << MATCHLEVEL_L << ": allow type match between DUP and INS, which takes effect by '-S' option" << endl;
//	cout << "             " << MATCHLEVEL_S << ": strict type match which is disabled by default" << endl;
//	cout << "             The default enabled match mode is 'loose' to allow the type match between DUP and INS." << endl;
//	cout << "   -C STR    Chromosomes to be processed: [null]" << endl;
//	cout << "             no decoy indicates not specifying the chromosome set for benchmarking." << endl;
//	cout << "             This parameter is used to specify the chromosomes to be benchmarked." << endl;
//	cout << "             Chromosome names should match the format within the VCF file. " <<endl;
//	cout << "             Chromosome names are separated by ';'. Example: -C \"1;2;3\" " << endl;
//	cout << "   -s INT    overlap extend size: [" << EXTEND_SIZE << "]" << endl;
//	cout << "   -i FLOAT  minimal sequence similarly for variant match: [" << SEQ_CONSISTENCY << "]" << endl;
//	cout << "   -a FLOAT  minimal sequence similarly for allelic variants match: [" << ALLELE_SEQ_CONSISTENCY << "]" << endl;
//	cout << "   -p FLOAT  minimal percent size ratio for variant match: [" << SVLEN_RATIO << "]" << endl;
//	cout << "   -t INT    number of threads [0]. 0 for the maximal number of threads" << endl;
//	cout << "             in machine" << endl;
//	cout << "   -T STR    Tool names [null]." << endl;
//	cout << "             This parameter is used for comparing multiple datasets. The number" << endl;
//	cout << "             of inputs should be consistent with the data set. Tool names are " <<endl;
//	cout << "             separated by ';'. Example: -T \"tool1;tool2;tool3\" " << endl;
//	cout << "   -o FILE   output directory: [" << outputPathname << "]" << endl;
//	cout << "   -l FILE   file name of long SV regions: [" << longSVFilename << "]" << endl;
//	cout << "   -r FILE   file name of benchmarking results to report: [" << htmlFilename << "]" << endl;
//	cout << "             Ensure that the filename extension is '.html'." << endl;
//	cout << "   -v        show version information" << endl;
//	cout << "   -h        show this help message and exit" << endl << endl;

	cout << "Example:" << endl;
	cout << "   # run the benchmarking on the user-called set (method) for a single sample to allow match between DUPs as INSs" << endl;
	cout << "   $ asvbm stat -T method user_sv.vcf benchmark_sv.vcf ref.fa" << endl << endl;

	cout << "   # run the benchmarking on the user-called set (method) for a single sample to perform the strict type matching by '-S' option" << endl;
	cout << "   $ asvbm stat -T method -S user_sv.vcf benchmark_sv.vcf ref.fa" << endl << endl;

	cout << "   # run the benchmarking on the user-called sets (tool1, tool2 and tool3) for multiple user callsets" << endl;
	cout << "   $ asvbm create -T \"tool1;tool2;tool3\" user_sv1.vcf user_sv2.vcf user_sv3.vcf benchmark_sv.vcf ref.fa" << endl;
}

// parse the parameters for convert command
int parseConvert(int argc, char **argv)
{
	int opt;
	string sv_format, ref_file, in_file, out_file, remove_mate_str = "1", mate_filename = "", remove_snv_str = "0", snv_filename = "";

	while( (opt = getopt(argc, argv, ":f:o:r:R:s:S:h")) != -1 ){
		switch(opt){
			case 'f': sv_format = optarg; break;
			case 'o': outputPathname = optarg; break;
			case 'r': remove_mate_str = optarg; break;
			case 'm': mate_filename = optarg; break;
			case 's': remove_snv_str = optarg; break;
			case 'S': snv_filename = optarg; break;
			case 'h': showUsageConvert(); exit(0);
			case '?': cout << "unknown option -" << (char)optopt << endl; exit(1);
			case ':': cout << "the option -" << (char)optopt << " needs a value" << endl; exit(1);
		}
	}

	if(remove_mate_str.compare("1")==0) {
		if(mate_filename.size()==0) mate_filename = mateItemFilename;
	}else if(remove_mate_str.compare("0")==0) mate_filename = "";
	else{
		cout << "Error: Please specify the correct value for '-r' option" << endl << endl;
		showUsageConvert();
		return 1;
	}

	if(remove_snv_str.compare("1")==0) {
		if(snv_filename.size()==0) snv_filename = snvFilename;
	}else if(remove_snv_str.compare("0")==0) snv_filename = "";
	else{
		cout << "Error: Please specify the correct value for '-s' option" << endl << endl;
		showUsageConvert();
		return 1;
	}
	if(outputPathname.at(outputPathname.size()-1)!='/') outputPathname += "/";

	opt = argc - optind; // the number of SAMs on the command line
	if(opt==3) {
		ref_file = argv[optind];
		in_file = argv[optind+1];
		out_file = argv[optind+2];
	}else { showUsageConvert(); return 1; }

	if(ref_file.size()>0 and (ref_file.compare(in_file)==0 or ref_file.compare(out_file)==0)){
		cout << "Please specify different reference file name" << endl << endl;
		showUsageConvert();
		return 1;
	}

	if(in_file.size()>0 and in_file.compare(out_file)==0){
		cout << "Please specify different SV file names" << endl << endl;
		showUsageConvert();
		return 1;
	}

	if(sv_format.compare("bed")==0 or sv_format.compare("vcf")==0 or sv_format.compare("csv")==0 or sv_format.compare("nm")==0) // nm: private usage
		convert(in_file, out_file, ref_file, mate_filename, snv_filename, sv_format);
	else{
		cout << "Error: Please specify the correct SV file format" << endl << endl;
		showUsageConvert();
		return 1;
	}

	return 0;
}

void convert(string &infilename, string &outfilename, string &reffilename, string &mate_filename, string &snv_filename, string &sv_format){
	string label = "benchmark";
	mkdir(outputPathname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	convertScreenFilename = outputPathname + convertScreenFilename;
	outConvertScreenFile.open(convertScreenFilename);
	if(!outConvertScreenFile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << convertScreenFilename << endl;
		exit(1);
	}
	outConvertScreenFile << "Program command: " << program_cmd_str << endl << endl;

	printConvertParas(infilename, outfilename, reffilename, mate_filename, snv_filename, sv_format); // print parameters

	cout << "############# Convert variants: #############" << endl;
	outConvertScreenFile << "############# Convert variants: #############" << endl;

	outfilename = outputPathname + outfilename;
	if(sv_format.compare("bed")==0)
		convertBed(infilename, outfilename, reffilename, mate_filename, snv_filename , label);
	else if(sv_format.compare("vcf")==0)
		convertVcf(infilename, outfilename, reffilename, mate_filename, snv_filename, label);
	else if(sv_format.compare("csv")==0)
		convertCsv(infilename, outfilename, reffilename, mate_filename, snv_filename, label);
	else if(sv_format.compare("nm")==0)  // private usage: nm
		convertNm(infilename, outfilename, reffilename, mate_filename, snv_filename, label);

	outConvertScreenFile.close();
}
int parseStat(int argc, char **argv){
	int32_t opt, threadNum_tmp, longindex = 0;
	string sv_file1, ref_file, sv_file2, ToolNameStore, ChromosomeNames;
	vector<string> sv_files1;
	vector<string> tool_names;

	if (argc < 2) { showUsageStat(); return 1; }
	minsvlen = MIN_SVLEN;
	maxValidRegThres = MAX_VALID_REG_THRES;
	percentSeqSim = SEQ_CONSISTENCY;
	percentAlleleSeqSim = ALLELE_SEQ_CONSISTENCY;
	extendSize = EXTEND_SIZE;
	threadNum_tmp = 0;
	minSizeLargeSV = MIN_SIZE_LARGE_SV;
	extendSizeLargeSV = EXTEND_SIZE_LARGE_SV;
	svlenRatio = SVLEN_RATIO;
	typeMatchLevel = MATCHLEVEL_L;
	ToolNameStore = ChromosomeNames = BamFilePath = "";
	BamFileSign = false;
	static struct option long_options[] = {
		{"bench_refine", no_argument, 0, 'b'},
		{"help", no_argument, 0, 'h'},
		{"version", no_argument, 0, 'v'},
		{"minsvlen",no_argument, 0, 'm'},
		{0, 0, 0, 0}
	};
	while( (opt = getopt_long(argc, argv, ":m:M:s:t:r:o:i:a:p:l:T:B:R:hvC:bS", long_options, &longindex)) != -1 ){
		switch(opt){
			case 'm': minsvlen = stoi(optarg); break;
			case 'M': maxValidRegThres = stoi(optarg); break;
			case 's': extendSize = stoi(optarg); break;
			case 't': threadNum_tmp = stoi(optarg); break;
			case 'o': outputPathname = optarg; break;
			case 'i': percentSeqSim = stod(optarg); break;
			case 'p': svlenRatio = stod(optarg); break;
			case 'a': percentAlleleSeqSim = stod(optarg); break;
			case 'l': longSVFilename = optarg; break;
			case 'T': ToolNameStore = optarg; break;
			case 'B': BamFilePath = optarg;	BamFileSign = true;	break;
			case 'R': size_div_vec = parseSemicolonSeparatedValues(optarg); break;
			case 'r': htmlFilename = optarg; break;
			case 'C': ChromosomeNames = optarg; break;
			case 'h': showUsageStat(); exit(0);
			case 'v': show_version();; exit(0);
			case 'b': refine_bench_flag = true; break;
			case 'S': typeMatchLevel = MATCHLEVEL_S; break;
			case '?': cout << "unknown option -" << (char)optopt << endl; exit(1);
			case ':': cout << "the option -" << (char)optopt << " needs a value" << endl; exit(1);
		}
	}

	if(maxValidRegThres<0){
		cout << "Error: Please specify the correct positive value for option: -m" << endl;
		showUsageStat();
		return 1;
	}
	if(minsvlen<=0){
		cout << "Error: Please specify the correct positive value for option: -M" << endl;
		showUsageStat();
		return 1;
	}
	if(extendSize<=0){
		cout << "Error: Please specify the correct positive value for option: -s" << endl;
		showUsageStat();
		return 1;
	}
	if(ToolNameStore.size()<0){
		cout << "Error: Specify the correct tool name for the option, separated by ';' :-T" << endl;
		showUsageStat();
		return 1;
	}
	if(ChromosomeNames.size()<0){
		cout << "Error: Specify the correct set of chromosomes for the option, separated by ';' :-C" << endl;
		showUsageStat();
		return 1;
	}

	/*if(typeMatchLevel.compare("strict")==0 or typeMatchLevel.compare("loose")==0);
	else {
		cout << "Error: Specify a type matching mode for the option: -L" << endl;
		showUsageStat();
		return 1;
	}*/
	if(threadNum_tmp==0) num_threads = sysconf(_SC_NPROCESSORS_ONLN);
	else num_threads = (threadNum_tmp>=sysconf(_SC_NPROCESSORS_ONLN)) ? sysconf(_SC_NPROCESSORS_ONLN) : threadNum_tmp;

	if(outputPathname.at(outputPathname.size()-1)!='/') outputPathname += "/";
	MeticsValues4_num.resize(size_div_vec.size()+1);
	region = generateSVregion();

	tool_names = split(ToolNameStore,";");
	chromosomeSet = split(ChromosomeNames,";");

	opt = argc - optind; // the number of SAMs on the command line
	if(tool_names.size()>=1){ //use the -T parameter to enter multiple tool names
		for(int i = optind; i<argc; i++){
			if(i == argc-2) sv_file2 = argv[i];
			else if(i == argc-1) ref_file = argv[i];
			else sv_files1.push_back(argv[i]);
		}
		if(tool_names.size() != sv_files1.size()){
			cout << "Error: The number of tool names does not match the data set :-T" << endl;
			showUsageStat();
			return 1;
		}
	}else if(opt==3){ //evaluate a data set
		if(opt==3) {
			sv_file1 = argv[optind];
			sv_file2 = argv[optind+1];
			ref_file = argv[optind+2];
		}else { showUsageStat(); return 1; }
	}else{ //evaluate multiple data sets without using the "-T" parameter
		for(int i = optind ; i<argc; i++){
			if(i == argc-2) sv_file2 = argv[i];
			else if(i == argc-1) ref_file = argv[i];
			else sv_files1.push_back(argv[i]);
		}
	}
	if(BamFileSign){
		if(!check_bam_average_read_length()){
			cout << "Error: The reads length are too short to support local joint analysis verification." << endl;
			showUsageStat();
			return 1;
		}
	}

	mem_total = getMemInfo("MemTotal", 2);
	swap_total = getMemInfo("SwapTotal", 2);
	if(mem_total<0 or swap_total<0){
		cerr << "line=" << __LINE__ << ", mem_total=" << mem_total << ", swap_total=" << swap_total << ", cannot get the supplied memory information, error." << endl;
		exit(1);
	}
	extend_total = mem_total<swap_total ? mem_total : swap_total;

	//cout << "mem_total=" << mem_total << ", swap_total=" << swap_total << ", extend_total=" << extend_total << endl;

	SVStatOp(ref_file, sv_file1, sv_file2, sv_files1, tool_names);

	return 0;
}

//int parseState(int argc, char **argv){
//	int32_t opt, threadNum_tmp;
//	string sv_file1, ref_file, sv_file2, ToolNameStore, ChromosomeNames;
//	vector<string> sv_files1;
//	vector<string> tool_names;
//
//	if (argc < 2) { showUsageStat(); return 1; }
//
//	maxValidRegThres = MAX_VALID_REG_THRES;
//	percentSeqSim = SEQ_CONSISTENCY;
//	percentAlleleSeqSim = ALLELE_SEQ_CONSISTENCY;
//	extendSize = EXTEND_SIZE;
//	threadNum_tmp = 0;
//	minSizeLargeSV = MIN_SIZE_LARGE_SV;
//	extendSizeLargeSV = EXTEND_SIZE_LARGE_SV;
//	svlenRatio = SVLEN_RATIO;
//	typeMatchLevel = MATCHLEVEL_L;
//	ToolNameStore = ChromosomeNames = "";
//
//	while( (opt = getopt(argc, argv, ":m:s:t:r:o:i:a:p:l:T:h:C:S")) != -1 ){
//		switch(opt){
//			case 'm': maxValidRegThres = stoi(optarg); break;
//			case 's': extendSize = stoi(optarg); break;
//			case 't': threadNum_tmp = stoi(optarg); break;
//			case 'o': outputPathname = optarg; break;
//			case 'i': percentSeqSim = stod(optarg); break;
//			case 'p': svlenRatio = stod(optarg); break;
//			case 'a': percentAlleleSeqSim = stod(optarg); break;
//			case 'l': longSVFilename = optarg; break;
//			case 'T': ToolNameStore = optarg; break;
//			case 'r': htmlFilename = optarg; break;
//			case 'C': ChromosomeNames = optarg; break;
//			case 'h': showUsageStat(); exit(0);
//			case 'S': typeMatchLevel = MATCHLEVEL_S; break;
//			case '?': cout << "unknown option -" << (char)optopt << endl; exit(1);
//			case ':': cout << "the option -" << (char)optopt << " needs a value" << endl; exit(1);
//		}
//	}
//
//	if(maxValidRegThres<0){
//		cout << "Error: Please specify the correct positive value for option: -m" << endl;
//		showUsageStat();
//		return 1;
//	}
//	if(extendSize<=0){
//		cout << "Error: Please specify the correct positive value for option: -s" << endl;
//		showUsageStat();
//		return 1;
//	}
//	if(ToolNameStore.size()<0){
//		cout << "Error: Specify the correct tool name for the option, separated by ';' :-T" << endl;
//		showUsageStat();
//		return 1;
//	}
//	if(ChromosomeNames.size()<0){
//		cout << "Error: Specify the correct set of chromosomes for the option, separated by ';' :-C" << endl;
//		showUsageStat();
//		return 1;
//	}
//
//	/*if(typeMatchLevel.compare("strict")==0 or typeMatchLevel.compare("loose")==0);
//	else {
//		cout << "Error: Specify a type matching mode for the option: -L" << endl;
//		showUsageStat();
//		return 1;
//	}*/
//	if(threadNum_tmp==0) num_threads = sysconf(_SC_NPROCESSORS_ONLN);
//	else num_threads = (threadNum_tmp>=sysconf(_SC_NPROCESSORS_ONLN)) ? sysconf(_SC_NPROCESSORS_ONLN) : threadNum_tmp;
//
//	if(outputPathname.at(outputPathname.size()-1)!='/') outputPathname += "/";
//
//	tool_names = split(ToolNameStore,";");
//	chromosomeSet = split(ChromosomeNames,";");
//
//	opt = argc - optind; // the number of SAMs on the command line
//	if(tool_names.size()>=1){ //use the -T parameter to enter multiple tool names
//		for(int i = optind; i<argc; i++){
//			if(i == argc-2) sv_file2 = argv[i];
//			else if(i == argc-1) ref_file = argv[i];
//			else sv_files1.push_back(argv[i]);
//		}
//		if(tool_names.size() != sv_files1.size()){
//			cout << "Error: The number of tool names does not match the data set :-T" << endl;
//			showUsageStat();
//			return 1;
//		}
//	}else if(opt==3){ //evaluate a data set
//		if(opt==3) {
//			sv_file1 = argv[optind];
//			sv_file2 = argv[optind+1];
//			ref_file = argv[optind+2];
//		}else { showUsageStat(); return 1; }
//	}else{ //evaluate multiple data sets without using the "-T" parameter
//		for(int i = optind ; i<argc; i++){
//			if(i == argc-2) sv_file2 = argv[i];
//			else if(i == argc-1) ref_file = argv[i];
//			else sv_files1.push_back(argv[i]);
//		}
//	}
//
//	mem_total = getMemInfo("MemTotal", 2);
//	swap_total = getMemInfo("SwapTotal", 2);
//	if(mem_total<0 or swap_total<0){
//		cerr << "line=" << __LINE__ << ", mem_total=" << mem_total << ", swap_total=" << swap_total << ", cannot get the supplied memory information, error." << endl;
//		exit(1);
//	}
//	extend_total = mem_total<swap_total ? mem_total : swap_total;
//
//	//cout << "mem_total=" << mem_total << ", swap_total=" << swap_total << ", extend_total=" << extend_total << endl;
//
//	SVStatOp(ref_file, sv_file1, sv_file2, sv_files1, tool_names);
//
//	return 0;
//}
// parse the parameters for create command
int parseCreate(int argc, char **argv){
	int32_t opt, threadNum_tmp, longindex = 0;;
	string sv_file1, ref_file, sv_file2, ToolNameStore, ChromosomeNames;
	vector<string> sv_files1;
	vector<string> tool_names;

	if (argc < 2) { showUsageCreate(); return 1; }
	minsvlen = MIN_SVLEN;
	maxValidRegThres = MAX_VALID_REG_THRES;
	percentSeqSim = SEQ_CONSISTENCY;
	percentAlleleSeqSim = ALLELE_SEQ_CONSISTENCY;
	extendSize = EXTEND_SIZE;
	threadNum_tmp = 0;
	minSizeLargeSV = MIN_SIZE_LARGE_SV;
	extendSizeLargeSV = EXTEND_SIZE_LARGE_SV;
	svlenRatio = SVLEN_RATIO;
	typeMatchLevel = MATCHLEVEL_L;
	ToolNameStore = ChromosomeNames = BamFilePath = "";
	BamFileSign = false;
	static struct option long_options[] = {
		{"help", no_argument, 0, 'h'},
		{"version", no_argument, 0, 'v'},
		{"minsvlen",no_argument, 0, 'm'},
		{0, 0, 0, 0}
	};
	while( (opt = getopt_long(argc, argv, ":m:M:s:t:r:o:i:a:p:l:T:B:R:hvC:S", long_options, &longindex)) != -1 ){
		switch(opt){
			case 'm': minsvlen = stoi(optarg); break;
			case 'M': maxValidRegThres = stoi(optarg); break;
			case 's': extendSize = stoi(optarg); break;
			case 't': threadNum_tmp = stoi(optarg); break;
			case 'o': outputPathname = optarg; break;
			case 'i': percentSeqSim = stod(optarg); break;
			case 'p': svlenRatio = stod(optarg); break;
			case 'a': percentAlleleSeqSim = stod(optarg); break;
			case 'l': longSVFilename = optarg; break;
			case 'T': ToolNameStore = optarg; break;
			case 'B': BamFilePath = optarg;	BamFileSign = true;	break;
			case 'R': size_div_vec = parseSemicolonSeparatedValues(optarg); break;
			case 'r': htmlFilename = optarg; break;
			case 'C': ChromosomeNames = optarg; break;
			case 'h': showUsageCreate(); exit(0);
			case 'v': show_version();; exit(0);
			case 'S': typeMatchLevel = MATCHLEVEL_S; break;
			case '?': cout << "unknown option -" << (char)optopt << endl; exit(1);
			case ':': cout << "the option -" << (char)optopt << " needs a value" << endl; exit(1);
		}
	}

	if(maxValidRegThres<0){
		cout << "Error: Please specify the correct positive value for option: -m" << endl;
		showUsageCreate();
		return 1;
	}
	if(minsvlen<=0){
		cout << "Error: Please specify the correct positive value for option: -M" << endl;
		showUsageStat();
		return 1;
	}
	if(extendSize<=0){
		cout << "Error: Please specify the correct positive value for option: -s" << endl;
		showUsageCreate();
		return 1;
	}
	if(ToolNameStore.size()<0){
		cout << "Error: Specify the correct tool name for the option, separated by ';' :-T" << endl;
		showUsageCreate();
		return 1;
	}
	if(ChromosomeNames.size()<0){
		cout << "Error: Specify the correct set of chromosomes for the option, separated by ';' :-C" << endl;
		showUsageCreate();
		return 1;
	}

	/*if(typeMatchLevel.compare("strict")==0 or typeMatchLevel.compare("loose")==0);
	else {
		cout << "Error: Specify a type matching mode for the option: -L" << endl;
		showUsageStat();
		return 1;
	}*/
	if(threadNum_tmp==0) num_threads = sysconf(_SC_NPROCESSORS_ONLN);
	else num_threads = (threadNum_tmp>=sysconf(_SC_NPROCESSORS_ONLN)) ? sysconf(_SC_NPROCESSORS_ONLN) : threadNum_tmp;

	if(outputPathname.at(outputPathname.size()-1)!='/') outputPathname += "/";
	MeticsValues4_num.resize(size_div_vec.size()+1);
	region = generateSVregion();

	tool_names = split(ToolNameStore,";");
	chromosomeSet = split(ChromosomeNames,";");

	opt = argc - optind; // the number of SAMs on the command line
	if(tool_names.size()>=1){ //use the -T parameter to enter multiple tool names
		for(int i = optind; i<argc; i++){
			if(i == argc-2) sv_file2 = argv[i];
			else if(i == argc-1) ref_file = argv[i];
			else sv_files1.push_back(argv[i]);
		}
		if(tool_names.size() != sv_files1.size()){
			cout << "Error: The number of tool names does not match the data set :-T" << endl;
			showUsageCreate();
			return 1;
		}
	}else if(opt==3){ //evaluate a data set
		if(opt==3) {
			sv_file1 = argv[optind];
			sv_file2 = argv[optind+1];
			ref_file = argv[optind+2];
		}else { showUsageCreate(); return 1; }
	}else{ //evaluate multiple data sets without using the "-T" parameter
		for(int i = optind ; i<argc; i++){
			if(i == argc-2) sv_file2 = argv[i];
			else if(i == argc-1) ref_file = argv[i];
			else sv_files1.push_back(argv[i]);
		}
	}
	if(BamFileSign){
		if(!check_bam_average_read_length()){
			cout << "Error: The reads length are too short to support local joint analysis verification." << endl;
			showUsageStat();
			return 1;
		}
	}

	mem_total = getMemInfo("MemTotal", 2);
	swap_total = getMemInfo("SwapTotal", 2);
	if(mem_total<0 or swap_total<0){
		cerr << "line=" << __LINE__ << ", mem_total=" << mem_total << ", swap_total=" << swap_total << ", cannot get the supplied memory information, error." << endl;
		exit(1);
	}
	extend_total = mem_total<swap_total ? mem_total : swap_total;

	//cout << "mem_total=" << mem_total << ", swap_total=" << swap_total << ", extend_total=" << extend_total << endl;

	SVStatOp(ref_file, sv_file1, sv_file2, sv_files1, tool_names);

	return 0;
}

//SVStatOp :Determine whether a single data set evaluation or a multi-data set evaluation
void SVStatOp(string &ref_file, string &sv_file1, string &sv_file2, vector<string> &sv_files1, vector<string> tool_names){
	string sv_file_name, mate_filename, snv_filename, convert_sv_file2;
	string sv_file_name_Path, convert_sv_file2_Path;
	string label = "benchmark", label1 = "user";
	vector<string> convert_sv_files;
	vector<SV_item*> sv_item_vec;
	mkdir(outputPathname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	suboutputDirnamePath = outputPathname + suboutputDirname;
	mkdir(suboutputDirnamePath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	//convert
	if(sv_files1.size()>=1){
		for(size_t i=0; i< sv_files1.size(); i++){
			sv_file_name = sv_files1[i];
			sv_file_name = PathqueryDot(sv_file_name);
			//path
			sv_file_name = Pathquerybackslash(sv_file_name);
			mate_filename = to_string(i)+ "_" + mateItemFilename;
//			snv_filename = to_string(i)+ "_" + snvFilename;
			sv_file_name_Path = suboutputDirnamePath + "/" + sv_file_name;
			convertVcf(sv_files1[i], sv_file_name_Path, ref_file, mate_filename, snv_filename, label1);
			convert_sv_files.push_back(sv_file_name_Path);

			sv_item_vec = loadDataWithoutTra(sv_file_name_Path);
			benchmark_vec.push_back(sv_item_vec);
			ToolsFilenames.push_back(sv_files1[i]);
		}
		//benchmark
		convert_sv_file2 = PathqueryDot(sv_file2);
		mate_filename = "Benckmark_" + mateItemFilename;
//		snv_filename = "Benckmark_" + snvFilename;
		convert_sv_file2 = Pathquerybackslash(convert_sv_file2);
		convert_sv_file2_Path = suboutputDirnamePath + "/" + convert_sv_file2;
		convertVcf(sv_file2, convert_sv_file2_Path, ref_file, mate_filename, snv_filename, label);
	}else{
		sv_file_name = sv_file1;
		sv_file_name = PathqueryDot(sv_file_name);
		mate_filename = mateItemFilename;
		//path
		sv_file_name = Pathquerybackslash(sv_file_name);
		sv_file_name_Path = suboutputDirnamePath + "/" + sv_file_name;
		convertVcf(sv_file1, sv_file_name_Path, ref_file, mate_filename, snv_filename, label1);
		//benchmark
		convert_sv_file2 = PathqueryDot(sv_file2);
		mate_filename = "Benckmark_" + mateItemFilename;
		convert_sv_file2 = Pathquerybackslash(convert_sv_file2);
		convert_sv_file2_Path = suboutputDirnamePath + "/" + convert_sv_file2;
		convertVcf(sv_file2, convert_sv_file2_Path, ref_file, mate_filename, snv_filename, label);
	}
	outConvertScreenFile.close();
	usersets_num = 0;

	if(convert_sv_files.size()>=1){	//multiple data sets are evaluated			sv_files1
		for(string& sv_file1 : convert_sv_files){
			SVStat(ref_file, sv_file1, convert_sv_file2_Path, convert_sv_files, tool_names);
		}
	}else
		SVStat(ref_file, sv_file_name_Path, convert_sv_file2_Path, sv_files1, tool_names);

	ResultPresentation(sv_files1, outputPathname, tool_names, outputBasicMetricschart, MeticsValues, MeticsValues1, ref_file);

}

// SV stat
void SVStat(string &ref_file, string &user_file, string &benchmark_file, vector<string> &sv_files1, vector<string> &tool_names){

//	mkdir(outputPathname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if(sv_files1.size()>=1){
		for(size_t i=1; i<=sv_files1.size();i++){
			if(user_file == sv_files1.at(i-1)){
				if(tool_names.size()>0) {outputInsideToolDirname = tool_names.at(i-1);	SVcallernames.push_back(tool_names.at(i-1));}
				else{
					//Instead of entering the tool name, use the data set name as the name of the output file
					size_t lastSlashPos = user_file.find_last_of('/');
					if (lastSlashPos == string::npos) {
						// If '/' is not found, the entire string is printed
						outputInsideToolDirname = sv_files1.at(i-1);
					} else {
						// When '/' is reached, the part after '/' is printed
						outputInsideToolDirname = user_file.substr(lastSlashPos + 1);
					}
					SVcallernames.push_back(outputInsideToolDirname);
				}
				alltoolnames.push_back(outputInsideToolDirname);
				method_name = outputInsideToolDirname;
				allmetric.push_back(outputInsideToolDirname);
				centerDistance.push_back(outputInsideToolDirname);
				sizeratio.push_back(outputInsideToolDirname);
				outputInsideToolDirname = outputPathname + outputInsideToolDirname;
				mkdir(outputInsideToolDirname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
				statScreenFilename = outputInsideToolDirname +'/' + statScreenFilename;
				outStatScreenFile.open(statScreenFilename);
				break;
			}
		}
	}else{
		alltoolnames.push_back("method");
		allmetric.push_back(acquiesce_count);
		centerDistance.push_back(acquiesce_count);
		sizeratio.push_back(acquiesce_count);
		statScreenFilename = outputPathname + statScreenFilename;
		outStatScreenFile.open(statScreenFilename);
	}

	if(!outStatScreenFile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << statScreenFilename << endl;
		exit(1);
	}
	outStatScreenFile << "Program command: " << program_cmd_str << endl << endl;

	printStatParas(user_file, benchmark_file, ref_file); // print parameters

	cout << "############# Stage 1: SV size statistics: #############" << endl;
	outStatScreenFile << "############# Stage 1: SV size statistics: #############" << endl;
	refRegSizeStat(user_file, benchmark_file, maxValidRegThres, minsvlen, sv_files1);

	cout << "\n\n############# Stage 2: Num statistics: #############" << endl;
	outStatScreenFile << "\n\n############# Stage 2: Num statistics: #############" << endl;
	SVNumStat(user_file, benchmark_file, ref_file,  maxValidRegThres, minsvlen, outputPathname, sv_files1);

	// compute SV size difference statistics
	cout << "\n\n############# Stage 3: SV size difference statistics: #############" << endl;
	outStatScreenFile << "\n\n############# Stage 3: SV size difference statistics: #############" << endl;
	SVSizeDifStat(user_file, benchmark_file, maxValidRegThres, minsvlen, sv_files1);

	cout << "\n\n############# Stage 4: SV size and num statistics: #############" << endl;
	outStatScreenFile << "\n\n############# Stage 4: SV size and num statistics: #############" << endl;
	SVSizeNumStat(user_file, benchmark_file, ref_file, maxValidRegThres, minsvlen, sv_files1);

	outStatScreenFile.close();
	//Evaluate multiple data sets and initialize the path for the next data save
	if(sv_files1.size()>=1)	FolderInit();
}

// print parameters for 'convert' command
void printConvertParas(string &infilename, string &outfilename, string &reffilename, string &mate_filename, string &snv_filename, string &sv_format){

	cout << "Program: " << PROG_NAME << " (" << PROG_DESC << ")" << endl;
	cout << "Version: " << PROG_VERSION << endl << endl;

	cout << "############# Parameters for 'convert' command: #############" << endl;
	cout << "          Input SV file: " << infilename << endl;
	cout << "         Output SV file: " << outfilename << endl;
	cout << "         Input Ref file: " << reffilename << endl;

	if(mate_filename.size())
		cout << "mate variant items: " << mate_filename << endl;
	if(snv_filename.size())
		cout << "              snv items: " << snv_filename << endl;

	cout << "       Output directory: " << outputPathname.substr(0, outputPathname.size()-1) << endl;

	if(sv_format.compare("bed")==0)
		cout << "      Input file format: BED" << endl;
	else if(sv_format.compare("vcf")==0)
		cout << "      Input file format: VCF" << endl;
	else if(sv_format.compare("csv")==0)
		cout << "      Input file format: CSV" << endl;
	else if(sv_format.compare("nm")==0)  // private usage: nm
		cout << "      Input file format: nm, [private usage]" << endl;
	cout << endl;

	// print to file
	outConvertScreenFile << "Program: " << PROG_NAME << " (" << PROG_DESC << ")" << endl;
	outConvertScreenFile << "Version: " << PROG_VERSION << endl << endl;

	outConvertScreenFile << "############# Parameters for 'convert' command: #############" << endl;

	outConvertScreenFile << "          Input SV file: " << infilename << endl;
	outConvertScreenFile << "         Output SV file: " << outfilename << endl;

	if(mate_filename.size())
		outConvertScreenFile << "mate variant items: " << mate_filename << endl;

	if(snv_filename.size())
			outConvertScreenFile << "snv items: " << snv_filename << endl;

	outConvertScreenFile << "       Output directory: " << outputPathname.substr(0, outputPathname.size()-1) << endl;

	if(sv_format.compare("bed")==0)
		outConvertScreenFile << "      Input file format: BED" << endl;
	else if(sv_format.compare("vcf")==0)
		outConvertScreenFile << "      Input file format: VCF" << endl;
	else if(sv_format.compare("csv")==0)
		outConvertScreenFile << "      Input file format: CSV" << endl;
	else if(sv_format.compare("nm")==0)  // private usage: nm
		outConvertScreenFile << "      Input file format: nm, [private usage]" << endl;
	outConvertScreenFile << endl;
}

// print parameters for 'stat' command
void printStatParas(string &user_file, string &benchmark_file ,string &ref_file){

	cout << "Program: " << PROG_NAME << " (" << PROG_DESC << ")" << endl;
	cout << "Version: " << PROG_VERSION << endl << endl;

	cout << "############# Parameters for 'ASVBM' command: #############" << endl;
	cout << "      Reference file: " << ref_file << endl;
	cout << " User-called SV file: " << user_file << endl;
	cout << "   Benchmark SV file: " << benchmark_file << endl;

	cout << " Maximal region size: " << maxValidRegThres << endl;
	cout << " Overlap extend size: " << extendSize << endl;
	cout << "  Percent size ratio: " << svlenRatio << endl;
	cout << "   Sequence similarly: " << percentSeqSim << endl;
	cout << "Allelic seq similarly: " << percentAlleleSeqSim << endl;
	cout << "   Number of threads: " << num_threads << endl;
	if(typeMatchLevel.compare(MATCHLEVEL_L)==0)
	cout << "    Type match level: " << typeMatchLevel << " (allow type match between DUPLICTION and INSERTION)" << endl;
	else
	cout << "    Type match level: " << typeMatchLevel << " (strict variant type match)" << endl;
	cout << "    Output directory: " << outputPathname.substr(0, outputPathname.size()-1) << endl;
	cout << "Long SV regions file: " << longSVFilename << endl << endl;

	// print to file
	outStatScreenFile << "Program: " << PROG_NAME << " (" << PROG_DESC << ")" << endl;
	outStatScreenFile << "Version: " << PROG_VERSION << endl << endl;

	outStatScreenFile << "############# Parameters for 'stat' command: #############" << endl;
	outStatScreenFile << "      Reference file: " << ref_file << endl;
	outStatScreenFile << " User-called SV file: " << user_file << endl;
	outStatScreenFile << "   Benchmark SV file: " << benchmark_file << endl;

	outStatScreenFile << " Maximal region size: " << maxValidRegThres << endl;
	outStatScreenFile << " Overlap extend size: " << extendSize << endl;
	outStatScreenFile << "  Percent size ratio: " << svlenRatio << endl;
	outStatScreenFile << "   Sequence similarly: " << percentSeqSim << endl;
	outStatScreenFile << "Allelic seq similarly: " << percentAlleleSeqSim << endl;
	outStatScreenFile << "   Number of threads: " << num_threads << endl;
	if(typeMatchLevel.compare(MATCHLEVEL_L)==0)
	outStatScreenFile << "    Type match level: " << typeMatchLevel << " (treat DUPLTCATION as INSERTION)" << endl;
	else
	outStatScreenFile << "    Type match level: " << typeMatchLevel << " (strict variant type matching)" << endl;
	outStatScreenFile << "    Output directory: " << outputPathname.substr(0, outputPathname.size()-1) << endl;
	outStatScreenFile << "Long SV regions file: " << longSVFilename << endl << endl;
}
