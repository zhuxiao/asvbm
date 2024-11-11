#include "paras.h"


int main(int argc, char *argv[]){

	program_cmd_str = getProgramCmdStr(argc, argv);
	all_cmd_program = getAllProgramCmdStr(argc, argv);

	if (argc < 2) { showUsageStat(); return 1; }
	string gnuplot_version = getProgramVersion("gnuplot --version");
	if(gnuplot_version.empty()){
		cerr << "Cannot find the 'gnuplot', please make sure it is correctly installed and the executable file 'gnuplot' or its soft link is included in the '$PATH' directory." << endl;
		exit(1);
	}

	if (strcmp(argv[1], "-h") == 0 or strcmp(argv[1], "help") == 0 or strcmp(argv[1], "--help") == 0) {
		if (argc == 2) { showUsageStat(); return 0; }
		argv++;
		argc = 2;
	}else if(strcmp(argv[1], "-v") == 0 or strcmp(argv[1], "--version") == 0){
		show_version();
		return 0;
	}

	return parseStat(argc, argv);
}

