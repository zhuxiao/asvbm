#include "paras.h"
#include "global.h"

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

// show usage for stat command
void showUsageStat(){

	cout << "Program: " << PROG_NAME << " (" << PROG_DESC << ")" << endl;
	cout << "Version: " << PROG_VERSION << endl << endl;

	cout << "Usage:  sv_stat [options] <USER_FILE> [<USER_FILE1>...] <BENCH_FILE> <REF_FILE>" << endl << endl;

	cout << "Description:" << endl;
	cout << "     USER_FILE    User called SV result file." << endl;
	cout << "     BENCH_FILE   Benchmark SV file." << endl;
	cout << "     REF_FILE     Reference file." << endl << endl;

	cout << "Options:" << endl;
	cout << "     -m INT       valid maximal region size for statistics: [" << MAX_VALID_REG_THRES << "]" << endl;
	cout << "                  0 is for all variant size are valid, and while positive" << endl;
	cout << "                  values are for the valid maximal region size, then longer" << endl;
	cout << "                  regions are omitted and saved to the file specified with '-l' option" << endl;
	cout << "     -L STR       Variant type matching pattern: [" << MATCHLEVEL_S << "]" << endl;
	cout << "                  " << MATCHLEVEL_S << ": strict variant type matching" << endl;
	cout << "                  " << MATCHLEVEL_L << ": treating duplications as insertions" << endl;
	cout << "     -C STR       Chromosome set to be evaluated: [null]" << endl;
	cout << "                  no decoy indicates not specifying the chromosome set for evaluation." << endl;
	cout << "                  This parameter is used to specify the chromosomes to be evaluated." << endl;
	cout << "                  Chromosome names should match the format within the VCF file. " <<endl;
	cout << "                  Chromosome names are separated by ';'. Example: -C \"1;2;3\" " << endl;
	cout << "     -s INT       overlap extend size: [" << EXTEND_SIZE << "]" << endl;
	cout << "     -t INT       number of threads [0]. 0 for the maximal number of threads" << endl;
	cout << "                  in machine" << endl;
	cout << "     -T STR       Tool names [null]. 0 indicates that the tool name is not entered." << endl;
	cout << "                  This parameter is used for comparing multiple datasets. The number" << endl;
	cout << "                  of inputs should be consistent with the data set. Tool names are " <<endl;
	cout << "                  separated by ';'. Example: -T \"tool1;tool2;tool3\" " << endl;
	cout << "     -o FILE      output directory: [" << outputPathname << "]" << endl;
	cout << "     -l FILE      file name for long SV regions: [" << longSVFilename << "]" << endl;
	cout << "     -v           show version" << endl;
	cout << "     -h           show this help message and exit" << endl;
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
		convertBed(infilename, outfilename, reffilename, mate_filename, snv_filename);
	else if(sv_format.compare("vcf")==0)
		convertVcf(infilename, outfilename, reffilename, mate_filename, snv_filename);
	else if(sv_format.compare("csv")==0)
		convertCsv(infilename, outfilename, reffilename, mate_filename, snv_filename);
	else if(sv_format.compare("nm")==0)  // private usage: nm
		convertNm(infilename, outfilename, reffilename, mate_filename, snv_filename);

	outConvertScreenFile.close();
}

// parse the parameters for stat command
int parseStat(int argc, char **argv){
	int32_t opt, threadNum_tmp;
	string sv_file1, ref_file, sv_file2, ToolNameStore, ChromosomeNames;
	vector<string> sv_files1;
	vector<string> tool_names;

	if (argc < 2) { showUsageStat(); return 1; }

	maxValidRegThres = MAX_VALID_REG_THRES;
	extendSize = EXTEND_SIZE;
	threadNum_tmp = 0;
	minSizeLargeSV = MIN_SIZE_LARGE_SV;
	extendSizeLargeSV = EXTEND_SIZE_LARGE_SV;
	svlenRatio = SVLEN_RATIO;
	typeMatchLevel = MATCHLEVEL_S;
	ToolNameStore = ChromosomeNames = "";

	while( (opt = getopt(argc, argv, ":m:s:t:o:l:T:h:L:C:")) != -1 ){
		switch(opt){
			case 'm': maxValidRegThres = stoi(optarg); break;
			case 's': extendSize = stoi(optarg); break;
			case 't': threadNum_tmp = stoi(optarg); break;
			case 'o': outputPathname = optarg; break;
			case 'l': longSVFilename = optarg; break;
			case 'T': ToolNameStore = optarg; break;
			case 'C': ChromosomeNames = optarg; break;
			case 'h': showUsageStat(); exit(0);
			case 'L': typeMatchLevel = optarg; break;
			case '?': cout << "unknown option -" << (char)optopt << endl; exit(1);
			case ':': cout << "the option -" << (char)optopt << " needs a value" << endl; exit(1);
		}
	}

	if(maxValidRegThres<0){
		cout << "Error: Please specify the correct positive value for option: -m" << endl;
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
	if(typeMatchLevel.compare("strict")==0 or typeMatchLevel.compare("loose")==0);
	else {
		cout << "Error: Specify a type matching mode for the option: -L" << endl;
		showUsageStat();
		return 1;
	}
	if(threadNum_tmp==0) num_threads = sysconf(_SC_NPROCESSORS_ONLN);
	else num_threads = (threadNum_tmp>=sysconf(_SC_NPROCESSORS_ONLN)) ? sysconf(_SC_NPROCESSORS_ONLN) : threadNum_tmp;

	if(outputPathname.at(outputPathname.size()-1)!='/') outputPathname += "/";

	tool_names = split(ToolNameStore,";");
	chromosomeSet = split(ChromosomeNames,";");

	opt = argc - optind; // the number of SAMs on the command line
	if(tool_names.size()>1){ //use the -T parameter to enter multiple tool names
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

	mem_total = getMemInfo("MemTotal", 2);
	swap_total = getMemInfo("SwapTotal", 2);
	if(mem_total<0 or swap_total<0){
		cerr << "line=" << __LINE__ << ", mem_total=" << mem_total << ", swap_total=" << swap_total << ", cannot get the supplied memory information, error." << endl;
		exit(1);
	}
	extend_total = mem_total<swap_total ? mem_total : swap_total;

	cout << "mem_total=" << mem_total << ", swap_total=" << swap_total << ", extend_total=" << extend_total << endl;

	SVStatOp(ref_file, sv_file1, sv_file2, sv_files1, tool_names);

	return 0;
}

//SVStatOp :Determine whether a single data set evaluation or a multi-data set evaluation
void SVStatOp(string &ref_file, string &sv_file1, string &sv_file2, vector<string> &sv_files1, vector<string> tool_names){
	string sv_file_name, mate_filename, snv_filename, convert_sv_file2;
	string sv_file_name_Path, convert_sv_file2_Path;
	vector<string> convert_sv_files;
	mkdir(outputPathname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	suboutputDirnamePath = outputPathname + suboutputDirname;
	mkdir(suboutputDirnamePath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	//convert
	if(sv_files1.size()>1){
		for(size_t i=0; i< sv_files1.size(); i++){
			sv_file_name = sv_files1[i];
			sv_file_name = PathqueryDot(sv_file_name);
			//path
			sv_file_name = Pathquerybackslash(sv_file_name);
			mate_filename = to_string(i)+ "_" + mateItemFilename;
			sv_file_name_Path = suboutputDirnamePath + "/" + sv_file_name;
			convertVcf(sv_files1[i], sv_file_name_Path, ref_file, mate_filename, snv_filename);
			convert_sv_files.push_back(sv_file_name_Path);
		}
		//benchmark
		convert_sv_file2 = PathqueryDot(sv_file2);
		mate_filename = "Benckmark_" + mateItemFilename;
		convert_sv_file2 = Pathquerybackslash(convert_sv_file2);
		convert_sv_file2_Path = suboutputDirnamePath + "/" + convert_sv_file2;
		convertVcf(sv_file2, convert_sv_file2_Path, ref_file, mate_filename, snv_filename);
	}else{
		sv_file_name = sv_file1;
		sv_file_name = PathqueryDot(sv_file_name);
		mate_filename = mateItemFilename;
		//path
		sv_file_name = Pathquerybackslash(sv_file_name);
		sv_file_name_Path = suboutputDirnamePath + "/" + sv_file_name;
		convertVcf(sv_file1, sv_file_name_Path, ref_file, mate_filename, snv_filename);
		//benchmark
		convert_sv_file2 = PathqueryDot(sv_file2);
		mate_filename = "Benckmark_" + mateItemFilename;
		convert_sv_file2 = Pathquerybackslash(convert_sv_file2);
		convert_sv_file2_Path = suboutputDirnamePath + "/" + convert_sv_file2;
		convertVcf(sv_file2, convert_sv_file2_Path, ref_file, mate_filename, snv_filename);
	}
	outConvertScreenFile.close();

	if(convert_sv_files.size()>1){	//multiple data sets are evaluated			sv_files1
		for(string& sv_file1 : convert_sv_files){
			SVStat(ref_file, sv_file1, convert_sv_file2_Path, convert_sv_files, tool_names);
		}
	}else
		SVStat(ref_file, sv_file_name_Path, convert_sv_file2_Path, sv_files1, tool_names);

	ResultPresentation(sv_files1, outputPathname, tool_names, outputBasicMetricschart, MeticsValues, MeticsValues1);

}

// SV stat
void SVStat(string &ref_file, string &user_file, string &benchmark_file, vector<string> &sv_files1, vector<string> &tool_names){

//	mkdir(outputPathname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if(sv_files1.size()>1){
		for(size_t i=1; i<=sv_files1.size();i++){
			if(user_file == sv_files1.at(i-1)){
				if(tool_names.size()>0) outputInsideToolDirname = tool_names.at(i-1);
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
				}
				outputInsideToolDirname = outputPathname + outputInsideToolDirname;
				mkdir(outputInsideToolDirname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
				statScreenFilename = outputInsideToolDirname +'/' + statScreenFilename;
				outStatScreenFile.open(statScreenFilename);
				break;
			}
		}
	}else{
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
	refRegSizeStat(user_file, benchmark_file, maxValidRegThres, sv_files1);

	cout << "\n\n############# Stage 2: Num statistics: #############" << endl;
	outStatScreenFile << "\n\n############# Stage 2: Num statistics: #############" << endl;
	SVNumStat(user_file, benchmark_file, ref_file,  maxValidRegThres, outputPathname, sv_files1);

	// compute SV size difference statistics
	cout << "\n\n############# Stage 3: SV size difference statistics: #############" << endl;
	outStatScreenFile << "\n\n############# Stage 3: SV size difference statistics: #############" << endl;
	SVSizeDifStat(user_file, benchmark_file, maxValidRegThres, sv_files1);

	cout << "\n\n############# Stage 4: SV size and num statistics: #############" << endl;
	outStatScreenFile << "\n\n############# Stage 4: SV size and num statistics: #############" << endl;
	SVSizeNumStat(user_file, benchmark_file, ref_file, maxValidRegThres, sv_files1);

	outStatScreenFile.close();
	//Evaluate multiple data sets and initialize the path for the next data save
	if(sv_files1.size()>1)	FolderInit();
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

	cout << "############# Parameters for 'stat' command: #############" << endl;
	cout << "      Reference file: " << ref_file << endl;
	cout << " User-called SV file: " << user_file << endl;
	cout << "   Benchmark SV file: " << benchmark_file << endl;

	cout << " Maximal region size: " << maxValidRegThres << endl;
	cout << " Overlap extend size: " << extendSize << endl;
	cout << "   Number of threads: " << num_threads << endl;
	cout << "   Type match level: " << typeMatchLevel << endl;
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
	outStatScreenFile << "   Number of threads: " << num_threads << endl;
	outStatScreenFile << "    Type match level: " << typeMatchLevel << endl;
	outStatScreenFile << "    Output directory: " << outputPathname.substr(0, outputPathname.size()-1) << endl;
	outStatScreenFile << "Long SV regions file: " << longSVFilename << endl << endl;
}
