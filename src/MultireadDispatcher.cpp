/*
 * MultireadDispatcher.cpp
 *
 *  Created on: Sep 17, 2012
 *      Author: spinelli
 *
 */

/**
 * This tool get information from a bowtie aligned file (bowtie output format with --concise option) and retrieve from it
 * the tags that are aligned to several places in the genome. The output is a file containing, for each tag position, the
 * fraction of signal that must be assigned to it. The current version of this tool spread the signal of a tag uniformly on
 * its respective positions.
 *
 * The output file format is : <chromosom>	<strand>	<position>	<fraction>
 */

#include "MultireadDispatcher.h"

using namespace std;

MultireadDispatcher::MultireadDispatcher( string out_dir) {

	outputDir = out_dir;
}

MultireadDispatcher::~MultireadDispatcher() {

}

/**
 * Parse the Bowtie result file to dispatch signal over multireads.
 *
 */
void MultireadDispatcher::analyzeTags( char* input_file_path, char** returned_output_file) {

	// open the file
	Rprintf("\nOpening source file : %s", input_file_path);
	ifstream input_file;
	input_file.open( input_file_path);

	// define map and vector used to classify tags
	map<int, vector<AlignedTag*>*>* tag_to_aligns = new map<int, vector<AlignedTag*>*>;

	// Parse the File and filter data
	Rprintf("\nParsing Source File...");
	string line;
	vector<AlignedTag*>* align_list = new vector<AlignedTag*>;
	int previous_tag_id = -1;
	char *token, str[200];
	int line_number = 0;
	while( getline( input_file, line)){

		// Read the tag strand
		size_t pt = line.find( ':');
		assert( pt != string::npos);
		char strand = line[pt - 1];
		// Read the tag ID
		int tag_id = atoi( line.substr( 0, pt - 1).c_str());
		strcpy( str, line.substr( pt + 2).c_str());
		// Read the tag chrom ID
		char* num = strtok( str, ",");
		assert( num != NULL);
		int chrom_id = atoi( num);
		// Read the tag position
		char* num2 = strtok(NULL, ",");
		assert( num2 != NULL);
		int position = atoi( num2);

		// create the AligedTag from the retrieved info
		AlignedTag* align = new AlignedTag( tag_id, strand, chrom_id, position);

		// Test if the tag is the same as the previous tag
		// If so, simply add the AlignedTag to the list
		if( previous_tag_id == tag_id){
			align_list->push_back( align);
		}
		// If not
		else{
			// save the AlignedTag list (apart for first run) to the map mapping tag name to AlignedTag list
			if( previous_tag_id >= 0){
				tag_to_aligns->insert( pair<int, vector<AlignedTag*>*>( previous_tag_id, align_list));
			}
			// create a new AlignedTag list with the current tag info
			align_list = new vector<AlignedTag*>;
			align_list->push_back( align);
		}

		// Memorize tag name
		previous_tag_id = tag_id;

		line_number++;
		if( line_number % 10000000 == 0){
			Rprintf("\nNumber of lines already parsed : %d", line_number);
		}
	}

	if( line_number == 0){
		throw "No data found from source file";
	}

	// Save the last block of lines
	if( previous_tag_id >= 0){
		tag_to_aligns->insert( pair<int, vector<AlignedTag*>*>( previous_tag_id, align_list));
	}

	Rprintf("\nParsing finished");
	input_file.close();

	// Retrieve the base name of the file (without extension and parent folder)
	baseOutputName = input_file_path;
	size_t extension_index = baseOutputName.find( ".bow");
	if( extension_index != string::npos){
		baseOutputName = baseOutputName.substr( 0, extension_index);
	}
	size_t slash_index = baseOutputName.rfind( "/");
	if( slash_index != string::npos){
		baseOutputName = baseOutputName.substr( slash_index + 1, baseOutputName.length());
	}

	// test if the output dir exists
	struct stat sb;
	if( stat( outputDir.c_str(), &sb) != 0){
		string command = "mkdir " + outputDir;
		system( command.c_str());
	}

	// outputing result
	string output_file_path_1 = outputDir + "/" + baseOutputName + "_uniformDispatch.txt";
	FILE *file_out_1 = fopen( output_file_path_1.c_str(), "wb");

	Rprintf("\nDispatching multiread signal.");
	Rprintf("\n   Output file is %s", output_file_path_1.c_str());

	// print the multireads after
	Rprintf("\n   Saving multireads");
	for( map<int, vector<AlignedTag*>*>::iterator tag_itr = tag_to_aligns->begin(); tag_itr != tag_to_aligns->end(); ++tag_itr){
		int nb_position = tag_itr->second->size();
		if( nb_position > 1){
			double frac = 1.0 / (double) nb_position;
			for( vector<AlignedTag*>::iterator vect_itr = tag_itr->second->begin(); vect_itr != tag_itr->second->end(); ++vect_itr){
				fprintf( file_out_1, "%s\t%c\t%d\t%.15g\n", chromosomName[(*vect_itr)->getChromID()].c_str(), (*vect_itr)->getStrand(), (*vect_itr)->getPosition(), frac);
			}
		}
	}
	fclose( file_out_1);

	strcpy( returned_output_file[0],output_file_path_1.c_str());

	Rprintf("\nFinished");

}

/**
 *  Load the Genome reference file to get the chromosom number and names
 *
 */
void MultireadDispatcher::loadReference( char* reference_path) {

	FILE *fi = NULL;

	//load reference file
	fi = fopen( reference_path, "r");

	// get the number of chromosom
	fscanf( fi, "%d", &chromosomNumber);
	Rprintf("\nChromosom number = %d", chromosomNumber);

	if( chromosomNumber <= 0){
		throw "No Chromosom definition found in reference.";
	}

	// get the length of each chromosom
	int* seqLen = new int[chromosomNumber];
	for( int i = 0; i < chromosomNumber; i++){
		fscanf( fi, "%d", &seqLen[i]);
	}
	// get the name of each chromosom
	chromosomName = new string[chromosomNumber];
	char* token = new char[20];
	for( int i = 0; i < chromosomNumber; i++){
		fscanf( fi, "%s", token);
		chromosomName[i] = token;
	}

	fclose( fi);
}

/**
 * Help on software usage
 *
 */
void MultireadDispatcher::usage() {
	Rprintf( "\n Usage: MultireadDispatcher [--experimentFileName|-e] <bowtie_aligned_file> [--analysisDir|-o] <output_folder> [--genomeReference|-r] <gen_ref>\n");
	Rprintf( "\n     - <bowtie_aligned_file> : name of the output file of Bowtie (.bow)\n");
	Rprintf( "\n     - <output_folder> : name of the output folder where result file would be written. This folder path is relative to working dir\n");
	Rprintf( "\n     - <gen_ref> : The path to the file containing the information on the genome : chromosom number, lengths and names\n");

}
