/*
 * benchmark_integration.h
 *
 *  Created on: Jan 5, 2025
 *      Author: node
 */

#ifndef SRC_BENCHMARK_INTEGRATION_H_
#define SRC_BENCHMARK_INTEGRATION_H_

#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <unordered_set>
#include <tuple>
#include <thread>
#include <future>
#include "extvab.h"
#include "gnuplotcall.h"
#include "num_stat.h"
#include "vcf_header.h"


void integrationBenchmark(string &outputFile, string &ref_file);
#endif /* SRC_BENCHMARK_INTEGRATION_H_ */
