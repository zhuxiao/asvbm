#include "paras.h"


int main(int argc, char *argv[]){
	string command;
	program_cmd_str = getProgramCmdStr(argc, argv);
	all_cmd_program = getAllProgramCmdStr(argc, argv);
	if (argc < 2) { show(); return 1; }
	string gnuplot_version = getProgramVersion("gnuplot --version");
	if(gnuplot_version.empty()){
		cerr << "Cannot find the 'gnuplot', please make sure it is correctly installed and the executable file 'gnuplot' or its soft link is included in the '$PATH' directory." << endl;
		exit(1);
	}

	if (strcmp(argv[1], "-h") == 0 or strcmp(argv[1], "help") == 0 or strcmp(argv[1], "--help") == 0) {
		if (argc == 2) { show(); return 0; }
		argv++;
		argc = 2;
	}else if(strcmp(argv[1], "-v") == 0 or strcmp(argv[1], "--version") == 0){
		show_version();
		return 0;
	}
	if (strcmp(argv[1], "create")==0){
		command = "create";
		create_bench_flag = true;
		return parseCreate(argc-1, argv+1);
//		return parseConvert(argc-1, argv+1);
	}else if(strcmp(argv[1], "stat")==0){
		command = "stat";
		return parseStat(argc-1, argv+1);
	}else{
		cerr << "invalid command: " << argv[1] << endl << endl;;
		show(); return 1;
	}

	return 0;


}

