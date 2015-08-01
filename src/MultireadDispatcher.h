/*
 * MultireadDispatcher.h
 *
 *  Created on: Sep 17, 2012
 *      Author: spinelli
 */

#ifndef MULTIREADDISPATCHER_H_
#define MULTIREADDISPATCHER_H_

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <map>
#include <stdio.h>
#include <sys/stat.h>
#include <string>
#include <vector>
#include "R.h"

#include "AlignedTag.h"

class MultireadDispatcher {
public:
	MultireadDispatcher( std::string);
	virtual ~MultireadDispatcher();
	void analyzeTags( char*,char**);
	static void usage();
	void loadReference( char*);

private:
	std::string outputDir;
	std::string baseOutputName;
	int chromosomNumber;
	std::string* chromosomName;
};

#endif /* MULTIREADDISPATCHER_H_ */
