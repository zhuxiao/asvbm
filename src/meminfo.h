#ifndef SRC_MEMINFO_H_
#define SRC_MEMINFO_H_

#include <iostream>
#include <mutex>
#include <chrono>
#include <condition_variable>

using namespace std;

// global variables
extern int64_t mem_total;		// in kB
//int64_t mem_free;		// in kB
//int64_t mem_available;	// in kB
extern int64_t swap_total;	// in kB
extern int64_t extend_total; // in kB

extern int64_t mem_seqAln;		// in kB
extern double mem_use_block_factor;
extern double extend_use_block_factor;
extern int32_t mem_wait_seconds;	// wait time

//extern int32_t won rk_num;

extern pthread_mutex_t mutex_mem;

int64_t getMemInfo(const char *name, int32_t index);


#endif /* SRC_MEMINFO_H_ */
