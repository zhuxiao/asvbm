#include "meminfo.h"

// global variables
int64_t mem_total = 0;		// in kB
//int64_t mem_free;		// in kB
//int64_t mem_available;	// in kB
int64_t swap_total = 0;	// in kB
int64_t extend_total = 0;	// in kB, the minimum of mem_total and swap_total

int64_t mem_seqAln = 0;		// in kB
double mem_use_block_factor = 0.9;
double extend_use_block_factor = 0.1;
int32_t mem_wait_seconds = 2;	// wait time

//int32_t work_num = 0;

pthread_mutex_t mutex_mem = PTHREAD_MUTEX_INITIALIZER;

// get memory information (in kB, 1024 bytes) by awk (or grep) command
int64_t getMemInfo(const char *name, int32_t index){
	int64_t info_int = -1;
	FILE *stream;
	char tmp[256], info[128];

	//sprintf(tmp, "cat /proc/meminfo | grep %s | awk '{print $%d}'", name, index);
	sprintf(tmp, "cat /proc/meminfo | awk '$1 ~/%s/' | awk '{print $%d}'", name, index);
	stream = popen(tmp, "r");
	if(fread(info, 1, sizeof(info), stream)>0){
		info_int = atol(info);
	}
	pclose(stream);

	return info_int;
}

