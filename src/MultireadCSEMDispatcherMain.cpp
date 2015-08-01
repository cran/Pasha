#include<ctime>
#include<cmath>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<string>
#include<iostream>
#include<fstream>
#include<vector>
#include<algorithm>
#include<sys/stat.h>
#include <limits.h>
#include "R.h"

#include "MultireadCSEMPosPair.h"
#include "MultireadCSEMBcTree.h"

using namespace std;

const int STRLEN=1024;

struct Alignment {
	string dir;
	int rid;
	int cid;
	int pos;
	double frac;

	Alignment() {
		cid = 0;
		dir = "";
		pos = 0;
		rid = 0;
		frac = 0.0;
	}

	Alignment( int cid, string dir, int pos, int rid, double frac = 0.0) {
		this->cid = cid;
		this->dir = dir;
		this->pos = pos;
		this->rid = rid;
		this->frac = frac;
	}

};

int ROUND = 0;
int iterationNumber = 200; // 200 by default

int windowSize, halfws;
int m, n, ns; // m chromosomes, n reads, ns tot # of alignments
int *seqLen;
string *chromosomName;

//nUniq unique reads, nMulti multi reads
int nUniq, nMulti;

//# of alignments for multireads
int nAligns;

vector<int> s;
vector<Alignment> uniq, multi;

bcTree **bctrees;

/**
 * Load genome information from genome reference file
 *
 */
void loadRef( char* refF) {
	FILE *fi = NULL;

	//load ref file
	Rprintf( "\nLoading genome information");
	fi = fopen( refF, "r");
	fscanf(fi, "%d", &m);
	Rprintf( "\n|-- Number of chromosomes: %d", m);
	seqLen = new int[m];
	Rprintf( "\n|-- Sequence lengths");
	for( int i = 0; i < m; i++){
		fscanf(fi, "%d", &seqLen[i]);
	}

	Rprintf( "\n|-- Chromosom names");
	char* token = new char[20];
	chromosomName = new string[m];
	for( int i = 0; i < m; i++){
		fscanf(fi, "%s", token);
		chromosomName[i] = token;
	}
	Rprintf( "\n|-- Done");
	fclose( fi);
}

/**
 * Load the data from bowtie aligned file
 *
 *
 */
void loadData( istream &fin) {

	Rprintf( "\nParsing data file");

	string line;
	int crid = -1;
	int rid;
	string dir;
	int pos, cid;
	char *num, *num2, str[STRLEN];
	char *token;
	vector<Alignment> aligns;

	n = nUniq = nMulti = 0;
	ns = nAligns = 0;
	s.clear();
	uniq.clear();
	multi.clear();

	crid = -1;
	aligns.clear();
	// Parse the file line
	while( getline( fin, line)){

		// Get the line tokens
		size_t pt = line.find( ':');
		// Strand
		assert( pt != string::npos);
		dir = line[pt - 1];
		// Read ID
		rid = atoi( line.substr( 0, pt - 1).c_str());
		strcpy( str, line.substr( pt + 2).c_str());
		// chromosome id (starts from 0)
		num = strtok( str, ",");
		assert( num != NULL);
		cid = atoi( num);
		// position
		num2 = strtok(NULL, ",");
		assert( num2 != NULL);
		pos = atoi( num2);

		// Classify the tag positions group
		if( crid != rid != 0){
			if( crid >= 0){
				if( aligns.size() == 1){
					uniq.insert( uniq.end(), aligns.begin(), aligns.end());
					++nUniq;
				}
				else{
					if( nAligns == 0){
						s.push_back( 0);
					}
					multi.insert( multi.end(), aligns.begin(), aligns.end());
					nAligns += aligns.size();
					s.push_back( nAligns);
					++nMulti;
				}
			}
			crid = rid;
			aligns.clear();
		}

		aligns.push_back( Alignment( cid, dir, pos, rid));

		++ns;
		if( ns % 10000000 == 0)
			Rprintf( "\n|-- Number of line already parsed: %d", ns);
	}

	// Classify the last group of tag in the file
	if( crid >=0 ){
		if( aligns.size() == 1){
			uniq.insert( uniq.end(), aligns.begin(), aligns.end());
			++nUniq;
		}
		else{
			if( nAligns == 0)
				s.push_back( 0);
			multi.insert( multi.end(), aligns.begin(), aligns.end());
			nAligns += aligns.size();
			s.push_back( nAligns);
			++nMulti;
		}
	}

	n = nUniq + nMulti;
	Rprintf( "\n|-- Parsing Finished!");
	Rprintf( "\n|-- Number of Uniread tags = %d", nUniq);
	Rprintf( "\n|-- Number of Multiread tags = %d", nMulti);
	Rprintf( "\n|-- Total number of tags = %d", n);
}

int constructBCTrees() {

	Rprintf( "\nComputing BCTrees...");
	int nPos;
	PosPair *arr = new PosPair[ns];

	for( int i = 0; i < nUniq; i++){
		arr[i].cid = uniq[i].cid;
		arr[i].pos = uniq[i].pos;
	}
	for( int i = 0; i < nAligns; i++){
		arr[i + nUniq].cid = multi[i].cid;
		arr[i + nUniq].pos = multi[i].pos;
	}
	sort( arr, arr + ns);
	nPos = unique( arr, arr + ns) - arr;

	if( m <= 0){
		error( "Number of chromosomes is less than 1!");
		return -1;
	}

	bctrees = new bcTree*[m];
	memset( bctrees, 0, sizeof(bcTree*) * m);

	int curpos = 0, curcid, curend;

	do{
		curcid = arr[curpos].cid;
		curend = curpos + 1;
		while( curend < nPos && curcid == arr[curend].cid)
			++curend;
		bctrees[curcid] = new bcTree( curend - curpos, arr + curpos);
		curpos = curend;
	} while( curpos < nPos);

	for( int i = 0; i < m; i++)
		if( bctrees[i] == NULL){
			bctrees[i] = new bcTree();
		}

	delete[] arr;
	Rprintf( "\nFinish Construct BCTrees");

	return 0;
}

/**
 * Dispatch the signal among all
 *
 */
int EM() {
	Rprintf( "\nDispatching signal...");
	double ll, ll_old;
	double genomeL;
	vector<double> rec;

	char cid;
	double tot;
	int seql, pos, lb, ub;

	genomeL = 0.0;
	for( int i = 0; i < m; i++)
		genomeL += seqLen[i];

	// Uniform Prior, initialize
	ll = -n * log( genomeL);

	for( int i = 0; i < nUniq; i++){
		cid = uniq[i].cid;
		pos = uniq[i].pos;
		seql = seqLen[cid];
		lb = max( -1, pos - halfws - 1);
		ub = min( seql - 1, pos + halfws);
		ll += log( double( ub - lb));
		uniq[i].frac = 1.0;
	}

	for( int i = 0; i < nMulti; i++){
		tot = 0.0;
		rec.clear();
		for( int j = s[i]; j < s[i + 1]; j++){
			cid = multi[j].cid;
			pos = multi[j].pos;
			seql = seqLen[cid];
			lb = max( -1, pos - halfws - 1);
			ub = min( seql - 1, pos + halfws);
			rec.push_back( double( ub - lb));
			tot += (ub - lb);
		}
		ll += log( tot);
		for( int j = s[i]; j < s[i + 1]; j++){
			multi[j].frac = rec[j - s[i]] / tot;
		}
	}
	if( ll-1 == ll){
		error( "Error : inconsistent value of initialized parameter : %f", ll);
		return -1;
	}
	Rprintf( "\nll_initial = %f", ll);

	ROUND = 0;

	do{
		ll_old = ll;

		time_t a, b;
		long tcall;

		a = time( NULL);
		tinsert = tcall = 0;

		//clear bctrees
		for( int i = 0; i < m; i++)
			bctrees[i]->clear();

		//M step  ??? Somehow M step @_@|||
		for( int i = 0; i < nUniq; i++){
			tcall++;
			bctrees[uniq[i].cid]->insert( uniq[i].pos, uniq[i].frac);
		}
		for( int i = 0; i < nMulti; i++)
			for( int j = s[i]; j < s[i + 1]; j++){
				tcall++;
				bctrees[multi[j].cid]->insert( multi[j].pos, multi[j].frac);
			}
		b = time( NULL);

		a = time( NULL);
		tcount = tcall = 0;

		//E step
		ll = 0.0;
		for( int i = 0; i < nUniq; i++){
			cid = uniq[i].cid;
			pos = uniq[i].pos;
			seql = seqLen[cid];
			lb = max( -1, pos - halfws - 1);
			ub = min( seql - 1, pos + halfws);
			ll += log( bctrees[cid]->count( ub) - bctrees[cid]->count( lb));
			uniq[i].frac = 1.0;

			tcall += 2;
		}
		for( int i = 0; i < nMulti; i++){
			tot = 0.0;
			rec.clear();
			for( int j = s[i]; j < s[i + 1]; j++){
				cid = multi[j].cid;
				pos = multi[j].pos;
				seql = seqLen[cid];
				lb = max( -1, pos - halfws - 1);
				ub = min( seql - 1, pos + halfws);
				double val = bctrees[cid]->count( ub) - bctrees[cid]->count( lb);
				rec.push_back( val);
				tot += val;

				tcall += 2;
			}
			ll += log( tot);
			for( int j = s[i]; j < s[i + 1]; j++){
				multi[j].frac = rec[j - s[i]] / tot;
			}
		}
		ll -= n * log( 1.0 * n);

		b = time( NULL);

		++ROUND;
		Rprintf( "\nROUND %d : %f; %f", ROUND, ll, ll_old);

	} while( ll > ll_old + 0.1 && ROUND < iterationNumber);

	return 0;
}

/**
 * Output the result to file
 *
 */
void output( string input_file_path, string output_dir, char** returned_output_file) {

	// Retrieve the base name of the file (without extension and suffix)
	string base_output_name = input_file_path;
	size_t extension_index = base_output_name.find( ".bow");
	if( extension_index != string::npos){
		base_output_name = base_output_name.substr( 0, extension_index);
	}
	size_t slash_index = base_output_name.rfind( "/");
	if( slash_index != string::npos){
		base_output_name = base_output_name.substr( slash_index + 1, base_output_name.length());
	}

	// test if the output dir exists
	struct stat sb;
	if( stat( output_dir.c_str(), &sb) != 0){
		string command = "mkdir " + output_dir;
		system( command.c_str());
	}

	// outputing result
	string output_file_path_1 = output_dir + "/" + base_output_name + "_csemDispatch.txt";
	Rprintf( "\nWriting result file : %s", output_file_path_1.c_str());

	FILE *fo = fopen( output_file_path_1.c_str(), "wb");
	for( int i = 0; i < nAligns; i++){
		fprintf( fo, "%s\t%s\t%d\t%.15g\n", chromosomName[multi[i].cid].c_str(), multi[i].dir.c_str(), multi[i].pos, multi[i].frac);
	}
	fclose( fo);

	strcpy( returned_output_file[0],output_file_path_1.c_str());
}

/**
 * Free memory before exiting
 *
 *
 */
void freeAll() {
	delete[] seqLen;
	for( int i = 0; i < m; i++){
		delete bctrees[i];
	}
	delete[] bctrees;
}

/**
 * Help on software usage
 *
 */
void usage() {

	Rprintf( "\nUsage: csem [--experimentFileName|-i] <bowtie_aligned_file> [--analysisDir|-o] <output_folder> [--genomeReference|-r] <gen_ref> [--windowSize|-w] <size> [--interation|-b] <interation>");
	Rprintf( "\n     - <bowtie_aligned_file> : name of the output file of Bowtie (.bow)");
	Rprintf( "\n     - <output_folder> : name of the output folder where result file would be written. This folder path is relative to working dir");
	Rprintf( "\n     - <gen_ref> : The path to the file containing the information on the genome : chromosom number, lengths and names");
	Rprintf( "\n     - <size> : the window size for the multiread allocator");
	Rprintf( "\n     - <iterations> : the number of iterations for the multiread allocator");
}



extern "C" {

using namespace std;

	 void C_CSEMDispatch ( char** r_file_name, char** r_output_dir, char** r_reference, int* r_window_size, int* r_iteration_number, char** returnedOutputFile) {

		// Test the passed arguments
		if( r_file_name == NULL || r_output_dir == NULL || r_reference == NULL || r_window_size == NULL || r_iteration_number == NULL){
			error("Error: NULL argument passed to function.");
			return;
		}

		// Convert the R Vector to suitable types
		char* file_name = r_file_name[0];
		char* output_dir = r_output_dir[0];
		char* reference = r_reference[0];
		windowSize = r_window_size[0];
		iterationNumber = r_iteration_number[0];

		//windowSize = atoi( argv[2]);
		assert( windowSize & 1);
		halfws = windowSize / 2;

		//UPPERBOUND = atoi( argv[3]);

		// Load the genome reference info
		if( reference != NULL){
			loadRef( reference);
		}
		else{
			error( "Error : No genome reference specified : %s", output_dir);
			usage();
			return;
		}

		// Execute the analysis
		if( file_name != NULL){
			if( output_dir != NULL){
				// Load the data
				ifstream fin;
				fin.open( file_name);
				loadData( fin);
				fin.close();

				// Construct the BCTrees (binary counting trees)
				if( constructBCTrees() < 0){
					error( "An issue occurred during BCTrees building");
					return;
				}

				int result = EM();

				if( result < 0){
					error( "An issue occurred during dispatching");
					return;
				}

				// Output the result
				output( file_name, output_dir, returnedOutputFile);
			}
			else{
				error( "Error : No output dir specified : %s", output_dir);
				usage();
			}
		}
		else{
			error( "Error : No input file specified : %s", file_name);
			usage();
		}

		freeAll();
	 }
}
