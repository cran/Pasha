/*
 * ArtifactRemover.cpp
 *
 *  Created on: Oct 15, 2012
 *      Author: spinelli
 */

/**
 * This tools aims to remove artifacts coming from multiread tags accumulation from a Bowtie output file
 * ('.bow' obtained using the --concise option).
 * An accumulation of tag at a given position is considered as an artifact if the number in the pile
 * is greater than the total number of tags in the experiment divided by the provided threshold, the result multiplied
 * by the number of repeats the tag know i.e.
 *
 *   #tags > ( #total_tags / threshold) * (#repeats)
 *
 * When an artifact is detected, only the limit number of tag at the position is kept, the exceeding tags are removed.
 */

#include "MultireadArtifactRemover.h"

using namespace std;

ArtifactRemover::ArtifactRemover( int verbose, char* input_file_path, char* output_path) {

	//#ifdef _WIN32
	//string PATH_SEPARATOR="\\";
	//#else
	string PATH_SEPARATOR="/";
	//#endif

	unireadArtifacts = 0;
	multireadArtifacts = 0;
	verbosity = verbose;

	// Retrieve the base name of the file (without extension and suffix)
	string base_input_file_path = input_file_path;
	string input_file_base_name = basename(input_file_path);
	size_t found = input_file_base_name.find( ".bow");
	if( found != string::npos){
		baseOutputName = input_file_base_name.substr( 0, found);
	}
	else{
		baseOutputName = input_file_base_name;
	}
	baseOutputName = output_path + PATH_SEPARATOR + baseOutputName;
}


ArtifactRemover::ArtifactRemover( int verbose, char* input_file_path) {

	unireadArtifacts = 0;
	multireadArtifacts = 0;
	verbosity = verbose;

	// Retrieve the base name of the file (without extension and suffix)
	string base_input_file_path = input_file_path;
	size_t found = base_input_file_path.find( ".bow");
	if( found != string::npos){
		baseOutputName = base_input_file_path.substr( 0, found);
	}
	else{
		baseOutputName = base_input_file_path;
	}
}

ArtifactRemover::~ArtifactRemover() {

}

/**
 *  Load the Genome reference file to get the chromosom number and names
 *
 */
void ArtifactRemover::loadReference( char* reference_path) {


	Rprintf( "\nReading Genome Reference file...");
	FILE *fi = NULL;

	//load reference file
	fi = fopen( reference_path, "r");

	// get the number of chromosom
	fscanf( fi, "%d", &chromosomNumber);
	Rprintf("\nGot chromosome number = %d", chromosomNumber);

	if( chromosomNumber <= 0){
		throw "No Chromosome definition found in reference.";
	}

	// get the length of each chromosome
	chromosomLength = new int[chromosomNumber];
	for( int i = 0; i < chromosomNumber; i++){
		fscanf( fi, "%d", &chromosomLength[i]);
	}
	Rprintf("\nGot chromosome lengths.");
	// get the name of each chromosome
	chromosomName = new string[chromosomNumber];
	char* token = new char[20];
	for( int i = 0; i < chromosomNumber; i++){
		fscanf( fi, "%s", token);
		chromosomName[i] = token;
	}
	Rprintf("\nGot chromosome names.");

	fclose( fi);

	Rprintf("\nDone.");
}

/**
 * Parse the Bowtie result file to load data
 *
 */
void ArtifactRemover::loadData( char* input_file_path) {

	// open the file
	Rprintf("\nOpening source file : %s", input_file_path);
	ifstream input_file;
	input_file.open( input_file_path);

	// Define a map used to classify tags by tag ID
	tagToAligns = new map<int, list<AlignedTag*>*>;

	// Define a map used to classify tags by chromosom
	chromToTags = new map<int, list<AlignedTag*>*>;

	// Parse the File and filter data
	Rprintf("\n|-- Parsing Source File...");
	list<AlignedTag*>* align_list = new list<AlignedTag*>;
	int previous_tag_id = -1;
	char *token, str[200];
	string line;
	int line_number = 0;
	int total_number_tag_position = 0;
	totalNumberTags = 0;
	while( getline( input_file, line)){

		// Read the tag strand
		size_t pt = line.find( ':');
		assert( pt != string::npos);
		char strand = (line.substr( pt-1, pt))[0];
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

		// create the AlignedTag from the retrieved info
		AlignedTag* align = new AlignedTag( tag_id, strand, chrom_id, position);
		total_number_tag_position++;

		// Add the tag position to its chromosom list
		map<int, list<AlignedTag*>*>::iterator chrom_tags = chromToTags->find( chrom_id);
		if( chrom_tags == chromToTags->end()){
			list<AlignedTag*>* new_list = new list<AlignedTag*>;
			align->setRepeat( 1);
			new_list->push_back( align);
			chromToTags->insert( pair<int, list<AlignedTag*>*>( chrom_id, new_list));
		}
		else{
			align->setRepeat( 1);
			chrom_tags->second->push_back( align);
		}

		// Test if the tag is the same as the previous tag
		// If so, simply add the AlignedTag to the current list
		if( previous_tag_id == tag_id){
			align_list->push_back( align);
		}
		// If tag changed
		else{
			// save the AlignedTag list (only for multiread and apart for first run) to the map
			if( previous_tag_id >= 0){
				int align_list_size = align_list->size();
				if( align_list_size > 1){
					for( list<AlignedTag*>::iterator tag_itr = align_list->begin(); tag_itr != align_list->end(); ++tag_itr){
						(*tag_itr)->setRepeat( align_list_size);
					}
					tagToAligns->insert( pair<int, list<AlignedTag*>*>( previous_tag_id, align_list));
				}
				totalNumberTags++;
			}
			// create a new AlignedTag list with the current tag info
			align_list = new list<AlignedTag*>;
			align_list->push_back( align);
		}

		// Memorize tag name
		previous_tag_id = tag_id;

		line_number++;
		if( line_number % 10000000 == 0){
			Rprintf("\n|---|--Number of lines already parsed : %d", line_number);
		}
	}

	// Save the last block of lines if it concerns a multiread
	int align_list_size = align_list->size();
	if( previous_tag_id >= 0 && align_list_size > 1){
		for( list<AlignedTag*>::iterator tag_itr = align_list->begin(); tag_itr != align_list->end(); ++tag_itr){
			(*tag_itr)->setRepeat( align_list_size);
		}
		tagToAligns->insert( pair<int, list<AlignedTag*>*>( previous_tag_id, align_list));
	}
	totalNumberTags++;

	Rprintf("\n|--Number of tags retrieved : %d", totalNumberTags);
	Rprintf("\n|--Number of tag positions retrieved : %d", total_number_tag_position);
	Rprintf("\n|--Number of chromosom lists : %d", chromToTags->size());
	Rprintf("\n|--Parsing finished");

	if( totalNumberTags == 0|| total_number_tag_position == 0){
		throw "The number of tag or tag position retrieved is null.";
	}
	input_file.close();

	// Sorting list of AlignedTag by chromosom
	Rprintf("\nSorting AlignedTag chromosom lists");
	for( map<int, list<AlignedTag*>*>::iterator list_itr = chromToTags->begin(); list_itr != chromToTags->end(); ++list_itr){
		Rprintf("\n|--Sorting Chrom ID %d:%d tag positions", list_itr->first, (list_itr->second)->size());
		list_itr->second->sort( compareAlignedTags);
	}
	Rprintf("\n|--Sorting finished");
}

/**
 * Remove artifacts from loaded data
 *
 */
void ArtifactRemover::removeArtifact( int artifact_ratio) {

	// open a write file to save statistics on artifacts
	string artifact_stat_path = baseOutputName + "_ArtifactStats.txt";
	FILE *artifact_stat_file = fopen( artifact_stat_path.c_str(), "wb");

	// Parse the list of chromosom
	Rprintf("\nRemoving artefact");
	double uniread_limit_of_artifact = totalNumberTags / (double) artifact_ratio;
	for( map<int, list<AlignedTag*>*>::iterator chrom_itr = chromToTags->begin(); chrom_itr != chromToTags->end(); ++chrom_itr){
		Rprintf("\n|-- Analyzing Chrom ID %d : %d tag positions", chrom_itr->first, (chrom_itr->second)->size());
		int old_position = -1;
		char old_strand = ' ';
		list<AlignedTag*> current_list;
		// Parse the list of tags of a chromosom
		int tags_number = 0;
		for( list<AlignedTag*>::iterator list_itr = (chrom_itr->second)->begin(); list_itr != (chrom_itr->second)->end(); ++list_itr){
			int current_position = (*list_itr)->getPosition();
			char current_strand = (*list_itr)->getStrand();

			// If the position and strand are the same as previous, accumulate the AlignedTags
			if( current_position == old_position && current_strand == old_strand){
				current_list.push_back( *list_itr);
			}
			// If the position or strand changed, analyze the previous AlignedTag accumulation and remove the exceeding AlignedTag if required
			else{
				if( old_position >= 0){
					// Check if the accumulation contains an artifact
					checkTagAccumulation( current_list, uniread_limit_of_artifact, artifact_stat_file);
				}
				// Clear the old list and create a new one
				current_list.clear();
				current_list.push_back( *list_itr);
				old_position = current_position;
				old_strand = current_strand;
			}

			tags_number++;
			if( tags_number % 10000000 == 0){
				Rprintf("\n|---|--Number of tags already parsed : %d", tags_number);
			}
		}
		// Check if the last accumulation contains an artifact
		checkTagAccumulation( current_list, uniread_limit_of_artifact, artifact_stat_file);
	}

	fclose( artifact_stat_file);

	Rprintf("\nFinal results : ");
	Rprintf("\nUniread artifacts found : %d", unireadArtifacts);
	Rprintf("\nMultiread artifacts found : %d", multireadArtifacts);
	Rprintf("\nRemoving Artifacts finished.");
}

/**
 * Check if the accumulation contains an artifact
 *
 */
void ArtifactRemover::checkTagAccumulation( list<AlignedTag*>& accumulation_list, int uniread_limit_of_artifact, FILE* artifact_stat_file) {

	int accumulation_list_size = accumulation_list.size();

	// Test if the list exceed the uniread artifact threshold
	if( accumulation_list_size > uniread_limit_of_artifact){
		// Retrieve the number of repeats concerned by the tags at this position
		list<AlignedTag*>::iterator tag_itr = accumulation_list.begin();
		int repeat_number = (*tag_itr)->getRepeat();
		if( repeat_number <= 0){
			Rprintf("\n|---|--Issue with tag in accumulation : No repeat number information. TagID = %d", (*tag_itr)->getTagID());
			return;
		}
		// If repeat number is 1, it means the accumulation is composed of unireads
		if( repeat_number == 1){
			unireadArtifacts++;
			return;
		}

		// Test if the accumulation exceed the limit number of tags for multiread artifact, if not bypass it
		int multiread_limit_of_artifact = uniread_limit_of_artifact * repeat_number;
		if( accumulation_list_size > multiread_limit_of_artifact){
			multireadArtifacts++;
			fprintf( artifact_stat_file, "%s\t%d\t%d\t%d\t%d\n", chromosomName[(*tag_itr)->getChromID()].c_str(), (*tag_itr)->getPosition(), accumulation_list_size, repeat_number,
					chromosomLength[(*tag_itr)->getChromID()]);
			if( verbosity > 0){
				Rprintf("\n|---|-- Artifact retrieved at position chrID %d : %d : %c", (*tag_itr)->getChromID(), (*tag_itr)->getPosition(), (*tag_itr)->getStrand());
				Rprintf("\n|---|---|--Accumulation size = %d", accumulation_list_size);
				Rprintf("\n|---|---|--Repeat number for that accumulation = %d", repeat_number);
				Rprintf("\n|---|---|--Uniread artifact limit = %d", uniread_limit_of_artifact);
				Rprintf("\n|---|---|--Multiread artifact limit = %d", multiread_limit_of_artifact);
			}

			// Parse the list of tag in the accumulation and remove the exceeding AlignedTags from the list of the tag positions
			// so that they will not be saved in the final output. By default, it keeps in the accumulation a number of tag equals
			// to the number of repeat.
			int count = 0;
			for( list<AlignedTag*>::iterator artifact_itr = accumulation_list.begin(); artifact_itr != accumulation_list.end(); ++artifact_itr){
				if( count >= repeat_number){
					int artifact_tag_id = (*artifact_itr)->getTagID();
					map<int, list<AlignedTag*>*>::iterator artifact_align_list_itr = tagToAligns->find( (*artifact_itr)->getTagID());
					if( artifact_align_list_itr != tagToAligns->end()){
						artifact_align_list_itr->second->remove( *artifact_itr);
					}
				}
				count++;
			}
		}
	}
}

/**
 * Output the list of remaining AlignedTags to a bowtie format file
 *
 */
void ArtifactRemover::outputResult( char** returned_output_file) {

	string output_file_path_1 = baseOutputName + "_NoA.bow";
	FILE *file_out_1 = fopen( output_file_path_1.c_str(), "wb");

	Rprintf("\nSaving result");
	for( map<int, list<AlignedTag*>*>::iterator tag_itr = tagToAligns->begin(); tag_itr != tagToAligns->end(); ++tag_itr){
		for( list<AlignedTag*>::iterator vect_itr = tag_itr->second->begin(); vect_itr != tag_itr->second->end(); ++vect_itr){
			fprintf( file_out_1, "%d%c:<%d,%d,0>\n", (*vect_itr)->getTagID(), (*vect_itr)->getStrand(), (*vect_itr)->getChromID(), (*vect_itr)->getPosition());
		}
	}
	fclose( file_out_1);
	strcpy( returned_output_file[0],output_file_path_1.c_str());
	Rprintf("\nResult saved.");
}

/**
 * Compare two AlignedTags to classify them by ascending order of chromosom and position
 *
 */
bool ArtifactRemover::compareAlignedTags( AlignedTag* first, AlignedTag* second) {

	int chrom1 = first->getChromID();
	int chrom2 = second->getChromID();

	if( chrom1 != chrom2){
		return chrom1 < chrom2;
	}
	else{
		int first_pos = first->getPosition();
		int second_pos = second->getPosition();
		if( first_pos != second_pos){
			return first_pos < second_pos;
		}
		else{
			return first->getStrand() == second->getStrand();
		}
	}
}

/**
 * Help on software usage
 *
 */
void ArtifactRemover::usage() {
	Rprintf( " Usage: ArtifactRemover [--experimentFileName|-e] <bowtie_aligned_file> [--artifactRatio|-a] <tag_ratio> [--genomeReference|-r] <gen_ref> [--verbosity|-v] <verbose>\n");
	Rprintf( "     - <bowtie_aligned_file> : name of the output file of Bowtie (.bow)\n");
	Rprintf( "     - <tag_ratio> : ratio of tag indicating the limit for artifact detection\n");
	Rprintf( "     - <gen_ref> : The path to the file containing the information on the genome : chromosom number, lengths and names\n");
	Rprintf( "     - <verbose> : the verbosity level. 0=info; 1=trace\n");
}
