#include "global.h"
#include "dataLoader.h"
#include "sv_stat.h"
#include "util.h"


int main(int argc, char *argv[]){

	string command;

	if (argc < 2) { showUsage(); return 1; }

    if (strcmp(argv[1], "-h") == 0 or strcmp(argv[1], "help") == 0 or strcmp(argv[1], "--help") == 0) {
        if (argc == 2) { showUsage(); return 0; }
        argv++;
        argc = 2;
    }

    if (strcmp(argv[1], "convert")==0){
    	command = "convert";
    	return parseConvert(argc-1, argv+1);
    }else if(strcmp(argv[1], "stat")==0){
    	command = "stat";
    	return parseStat(argc-1, argv+1);
    }else{
    	cerr << "invalid command: " << argv[1] << endl << endl;;
    	showUsage(); return 1;
    }

	return 0;
}

// show usage
void showUsage(){

	cout << "Usage:  sv_stat  <command> [options]" << endl << endl;

	cout << "Commands:" << endl;
	cout << "     convert      convert SV results to 6-column BED (or 9-column BEDPE) file format" << endl;
	cout << "     stat         compute the SV statistics" << endl;
}

// show usage for convert command
void showUsageConvert(){

	cout << "Usage:  sv_stat convert [options] <infile> <outfile>" << endl << endl;

	cout << "Options:" << endl;
	cout << "     -f STR       SV input file format (required):" << endl;
	cout << "                  bed: BED/BEDPE format" << endl;
	cout << "                  vcf: VCF format" << endl;
	cout << "                  csv: CSV format" << endl;
	cout << "     -d INT       remove duplicated variant items [1]: 1 for yes, 0 for no" << endl;
	cout << "     -h           show this help message and exit" << endl;
}

// show usage for stat command
void showUsageStat(){

	cout << "Usage:  sv_stat stat [options] <USER_SV_FILE> <BENCHMARK_SV_FILE>" << endl << endl;

	cout << "Description:" << endl;
	cout << "     USER_SV_FILE         User called SV result file." << endl;
	cout << "     BENCHMARK_SV_FILE    Benchmark SV file." << endl << endl;

	cout << "Options:" << endl;
	cout << "     -m           valid maximal region size for statistics: [0]" << endl;
	cout << "                  0 is for all variant size are valid, and while positive" << endl;
	cout << "                  values are for the valid maximal region size, then longer" << endl;
	cout << "                  regions are omitted and saved to the file specified with -l" << endl;
	cout << "     -o FILE      output path name for SV statistics: [output]" << endl;
	cout << "     -l FILE      file name for long SV regions: [long_sv_reg.bed]" << endl;
	cout << "     -h           show this help message and exit" << endl;
}

// parse the parameters for convert command
int parseConvert(int argc, char **argv)
{
	int opt;
	string sv_format, in_file, out_file, remove_dup_str = "1";
	bool remove_dup_falg;

	while( (opt = getopt(argc, argv, ":f:d:h")) != -1 ){
		switch(opt){
			case 'f': sv_format = optarg; break;
			case 'd': remove_dup_str = optarg; break;
			case 'h': showUsageConvert(); exit(0);
			case '?': cout << "unknown option -" << (char)optopt << endl; exit(1);
			case ':': cout << "the option -" << (char)optopt << " needs a value" << endl; exit(1);
		}
	}

	if(remove_dup_str.compare("1")==0) remove_dup_falg = true;
	else if(remove_dup_str.compare("0")==0) remove_dup_falg = false;
	else{
		cout << "Error: Please specify the correct value for '-d' option" << endl << endl;
		showUsageConvert();
		return 1;
	}

	opt = argc - optind; // the number of SAMs on the command line
	if(opt==2) {
		in_file = argv[optind];
		out_file = argv[optind+1];
	}else { showUsageConvert(); return 1; }

	if(sv_format.compare("bed")==0)
		convertBed(in_file, out_file, remove_dup_falg);
	else if(sv_format.compare("vcf")==0)
		convertVcf(in_file, out_file, remove_dup_falg);
	else if(sv_format.compare("csv")==0)
		convertCsv(in_file, out_file, remove_dup_falg);
	else if(sv_format.compare("nm")==0)  // private usage: nm
		convertNm(in_file, out_file, remove_dup_falg);
	else{
		cout << "Error: Please specify the correct SV file format" << endl << endl;
		showUsageConvert();
		return 1;
	}

	return 0;
}

// parse the parameters for stat command
int parseStat(int argc, char **argv){
	int32_t opt;
	string sv_file1, sv_file2;

	if (argc < 2) { showUsageStat(); return 1; }

	maxValidRegThres = 0;
	while( (opt = getopt(argc, argv, ":m:o:l:h")) != -1 ){
		switch(opt){
			case 'm': maxValidRegThres = stoi(optarg); break;
			case 'o': outputPathname = optarg; break;
			case 'l': longSVFilename = optarg; break;
			case 'h': showUsageStat(); exit(0);
			case '?': cout << "unknown option -" << (char)optopt << endl; exit(1);
			case ':': cout << "the option -" << (char)optopt << " needs a value" << endl; exit(1);
		}
	}

	if(maxValidRegThres<0){
		cout << "Error: Please specify the correct positive value for option: -m" << endl;
		showUsageStat();
		return 1;
	}
	if(outputPathname.at(outputPathname.size()-1)!='/') outputPathname += "/";

	opt = argc - optind; // the number of SAMs on the command line
	if(opt==2) {
		sv_file1 = argv[optind];
		sv_file2 = argv[optind+1];
	}else { showUsageStat(); return 1; }

	SVStat(sv_file1, sv_file2);

	return 0;
}

// SV stat
void SVStat(string &sv_file1, string &sv_file2){

	mkdir(outputPathname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	cout << "############# Phage 1: SV size statistics: #############" << endl;
	refRegSizeStat(sv_file1, sv_file2, maxValidRegThres);

	cout << "\n\n############# Phage 2: Num statistics: #############" << endl;
	SVNumStat(sv_file1, sv_file2, maxValidRegThres, outputPathname);

	// compute SV size difference statistics
	cout << "\n\n############# Phage 3: SV size difference statistics: #############" << endl;
	SVSizeDifStat(sv_file1, sv_file2, maxValidRegThres);

	cout << "\n\n############# Phage 4: SV type and num statistics: #############" << endl;
	SVTypeNumStat(sv_file1, sv_file2, maxValidRegThres);
}

