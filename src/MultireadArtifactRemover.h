/*
 * ArtifactRemover.h
 *
 *  Created on: Oct 15, 2012
 *      Author: spinelli
 */

#ifndef ARTIFACTREMOVER_H_
#define ARTIFACTREMOVER_H_

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
#include <list>
#include <libgen.h>
#include "R.h"

#include "AlignedTag.h"

class ArtifactRemover {
public:
	ArtifactRemover( int, char*, char*);
	ArtifactRemover( int, char*);
	virtual ~ArtifactRemover();
	static void usage();
	void loadReference( char*);
	void loadData( char*);
	void removeArtifact( int);
	void outputResult(char**);
	static bool compareAlignedTags( AlignedTag* , AlignedTag* ) ;

private:
	int verbosity;
	int chromosomNumber;
	int* chromosomLength;
	std::string* chromosomName;
	int totalNumberTags;
	int unireadArtifacts;
	int multireadArtifacts;
	std::string baseOutputName;
	std::map<int, std::list<AlignedTag*>*>* tagToAligns;
	std::map<int, std::list<AlignedTag*>*>* chromToTags;
	void checkTagAccumulation( std::list<AlignedTag*>&, int, FILE*);
};

#endif /* ARTIFACTREMOVER_H_ */
