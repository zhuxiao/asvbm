#include "paras.h"


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

