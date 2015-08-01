/*
 * WigRepeatAnalyzer.h
 *
 *  Created on: Oct 30, 2012
 *      Author: spinelli
 */

#ifndef WIGREPEATANALYZER_H_
#define WIGREPEATANALYZER_H_

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

//#include "regex.h"
#include "RepeatLocus.h"
#include "WigLocus.h"

class WigRepeatAnalyzer {
public:
	WigRepeatAnalyzer();
	virtual ~WigRepeatAnalyzer();
	void analyzeRepeat( std::string wig_path, std::string wig_file, const char* repeat_masker_file, const char* output_folder);
	static void usage();

private:
	RepeatLocus* readRepeatLocus( std::ifstream& repeatmasker_input_file);
	WigLocus* readWigLocus( std::ifstream& wig_input_file, WigLocus* last_locus, bool replace);
	void accumulateValue( std::map< std::string ,long int*> * _map, std::string key, long int value);

};

#endif /* WIGREPEATANALYZER_H_ */
