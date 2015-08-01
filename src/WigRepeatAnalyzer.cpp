/*
 * WigRepeatAnalyzer.cpp
 *
 *  Created on: Oct 30, 2012
 *      Author: spinelli
 */

/** This tool considers the signal provided in a WIG file and compare it to the locations provided in a RepeatMasker file.
 * The coverage of the WIG on the repeat locii and the part of the signal weight in them are computed. Statistics are produced
 * respect to the repeat classes (SINE, LINE...) and their sub-families. Statistics for positive signal and negative signal are separated.
 *
 * Output is provided in 4 files. All output files have the WIG filename as base name extended with a suffix:
 *
 * * _repeatClassCoverage.txt : contains the information about the coverage of the repeat classes on the WIG. Format:
 *   <class_name>		<signal_sign>		<class_coverage>		<total_coverage_sign>			<total_coverage_both_sign>
 *
 * * _repeatFamilyCoverage.txt : contains the information about the coverage of the repeat families on the WIG. Format:
 *   <family_name>		<signal_sign>		<family_coverage>		<total_coverage_sign>			<total_coverage_both_sign>
 *
 * * _repeatClassWeight.txt : contains the information about the weight of the repeat classes on the WIG. Format:
 *   <class_name>		<signal_sign>		<class_coverage>		<total_coverage_sign>			<total_coverage_both_sign>
 *
 * * _repeatFamilyWeight.txt : contains the information about the weight of the repeat families on the WIG. Format:
 *   <family_name>		<signal_sign>		<family_weight>		<total_weight_sign>			<total_weight_both_sign>
 */


#include "WigRepeatAnalyzer.h"
#include "R.h"

using namespace std;

WigRepeatAnalyzer::WigRepeatAnalyzer() {

}

WigRepeatAnalyzer::~WigRepeatAnalyzer() {

}

void WigRepeatAnalyzer::analyzeRepeat( string wig_path, string wig_file, const char* repeat_masker_file, const char* base_output_folder) {

	map<string, long int*>* class_size = new map<string, long int*>;
	map<string, long int*>* class_weight = new map<string, long int*>;
	map<string, long int*>* family_size = new map<string, long int*>;
	map<string, long int*>* family_weight = new map<string, long int*>;
	long int total_wig_size_pos = 0;
	long int total_wig_size_neg = 0;
	long int total_wig_weight_pos = 0;
	long int total_wig_weight_neg = 0;

	// Open repeat masker file and initialize repeat variables
	ifstream repeatmasker_input_file;
	repeatmasker_input_file.open( repeat_masker_file);
	RepeatLocus* current_repeat_locus = readRepeatLocus( repeatmasker_input_file);
	string old_repeat_chr = current_repeat_locus->getChromosom();
	double total_weight;

	// Open wig file and initialize wig variables
	ifstream wig_input_file;
	string wig_full_path = (string) wig_path + "/" + (string) wig_file;
	Rprintf( "Analyzing file : %s\n", wig_full_path.c_str());
	wig_input_file.open( wig_full_path.c_str());
	WigLocus* last_wig_locus = NULL;
	WigLocus* current_wig_locus = readWigLocus( wig_input_file, last_wig_locus, false);
	// Accumulate the size and weight of the wig locus
	if( current_wig_locus->getValue() > 0){
		total_wig_size_pos += current_wig_locus->getStep();
		total_wig_weight_pos += current_wig_locus->getStep() * current_wig_locus->getValue();
	}
	else if( current_wig_locus->getValue() < 0){
		total_wig_size_neg += current_wig_locus->getStep();
		total_wig_weight_neg -= current_wig_locus->getStep() * current_wig_locus->getValue();
	}

	// Compare Repeat Masker file entries and wig file entries
	bool change_chrom = false;
	do{
		// If there is a repeat locus
		if( current_repeat_locus != NULL){
			if( current_wig_locus != NULL){
				string current_repeat_chr = current_repeat_locus->getChromosom();
				string current_wig_chr = current_wig_locus->getChromosom();
				//Rprintf( "current_repeat_locus= %s:%d", current_repeat_locus->getChromosom().c_str(), current_repeat_locus->getStart());
				//Rprintf( "current_wig_locus= %s:%d", current_wig_locus->getChromosom().c_str(), current_wig_locus->getStart());

				// If chromosom are not the same, go ahead in the one of the file depending on which have changed chromsom
				if( current_repeat_chr.compare( current_wig_chr) != 0){
					if( !change_chrom){
						Rprintf( "|-- Need to change chromosom : WIG chromosom = %s; Repeat Chromosom = %s\n", current_wig_chr.c_str(), current_repeat_chr.c_str());
						change_chrom = true;
					}
					// If chromosom of repeat locus is constant, it means the wig file has changed chromosom and we
					// have to go ahead in the repeat file
					if( old_repeat_chr.compare( current_repeat_chr) == 0){
						delete current_repeat_locus;
						current_repeat_locus = readRepeatLocus( repeatmasker_input_file);
						continue;
					}
					// If the chromosom of repeat locus changed, it means we have to go ahead in the wig file
					else{
						// Remember old locus
						last_wig_locus = current_wig_locus;
						// Create new locus from the old one
						current_wig_locus = readWigLocus( wig_input_file, last_wig_locus, false);
						// Accumulate the size and weight of the wig locus
						if( current_wig_locus != NULL){
							if( current_wig_locus->getValue() > 0){
								total_wig_size_pos += current_wig_locus->getStep();
								total_wig_weight_pos += current_wig_locus->getStep() * current_wig_locus->getValue();
							}
							else if( current_wig_locus->getValue() < 0){
								total_wig_size_neg += current_wig_locus->getStep();
								total_wig_weight_neg -= current_wig_locus->getStep() * current_wig_locus->getValue();
							}
						}
						// Remove old locus
						delete last_wig_locus;
						continue;
					}
				}
				// If the chromosoms are the same, test the overlap
				else{
					if( change_chrom){
						Rprintf( "|---|-- Chromosom changed: WIG chromosom = %s; Repeat Chromosom = %s\n", current_wig_chr.c_str(), current_repeat_chr.c_str());
						change_chrom = false;
					}
					int r_start = current_repeat_locus->getStart();
					int r_end = current_repeat_locus->getEnd();
					int w_start = current_wig_locus->getStart();
					int w_end = current_wig_locus->getEnd();

					// The test below look if the repeat locus and the wig locus have an intersection
					// If so, accumulate corresponding scores
					if( !(r_end < w_start || r_start > w_end)){
						long int size = min( w_end, r_end) - max( w_start, r_start);
						long int weight = size * current_wig_locus->getValue();
						if( current_wig_locus->getValue() > 0){
							accumulateValue( class_size, current_repeat_locus->getClasse() + "+", size);
							accumulateValue( family_size, current_repeat_locus->getFamily() + "+", size);
							accumulateValue( class_weight, current_repeat_locus->getClasse() + "+", weight);
							accumulateValue( family_weight, current_repeat_locus->getFamily() + "+", weight);
						}
						else if( current_wig_locus->getValue() < 0){
							accumulateValue( class_size, current_repeat_locus->getClasse() + "-", size);
							accumulateValue( family_size, current_repeat_locus->getFamily() + "-", size);
							accumulateValue( class_weight, current_repeat_locus->getClasse() + "-", -weight);
							accumulateValue( family_weight, current_repeat_locus->getFamily() + "-", -weight);
						}
					}

					// Advance in the file of the locus that is at the most left position
					if( r_end <= w_end){
						old_repeat_chr = current_repeat_chr;
						delete current_repeat_locus;
						current_repeat_locus = readRepeatLocus( repeatmasker_input_file);
					}
					else{
						// Remember the old locus
						last_wig_locus = current_wig_locus;
						// Build the new locus from the old one
						current_wig_locus = readWigLocus( wig_input_file, last_wig_locus, false);
						// Accumulate the size and weight of the wig locus
						if( current_wig_locus != NULL){
							if( current_wig_locus->getValue() > 0){
								total_wig_size_pos += current_wig_locus->getStep();
								total_wig_weight_pos += current_wig_locus->getStep() * current_wig_locus->getValue();
							}
							else if( current_wig_locus->getValue() < 0){
								total_wig_size_neg += current_wig_locus->getStep();
								total_wig_weight_neg -= current_wig_locus->getStep() * current_wig_locus->getValue();
							}
						}
						// Remove the old locus
						delete last_wig_locus;
					}

				}
			}
		}
		// If there is no more repeat locus, simply advance on wig to count total sizes and total weight
		else{
			// Remember old locus
			last_wig_locus = current_wig_locus;
			// Create new locus from the old one
			current_wig_locus = readWigLocus( wig_input_file, last_wig_locus, false);
			// Accumulate the size and weight of the wig locus
			if( current_wig_locus != NULL){
				if( current_wig_locus->getValue() > 0){
					total_wig_size_pos += current_wig_locus->getStep();
					total_wig_weight_pos += current_wig_locus->getStep() * current_wig_locus->getValue();
				}
				else if( current_wig_locus->getValue() < 0){
					total_wig_size_neg += current_wig_locus->getStep();
					total_wig_weight_neg -= current_wig_locus->getStep() * current_wig_locus->getValue();
				}
			}
			// Remove old locus
			delete last_wig_locus;
		}

	} while( current_wig_locus != NULL);

	wig_input_file.close();
	repeatmasker_input_file.close();

	Rprintf("Analysis finished.\n");

	Rprintf("Building output file...\n");
	// Retrieve the base name of the file (without extension)
	string base_ouput_file = wig_file;
	size_t extension_index = base_ouput_file.find( ".wig");
	if( extension_index != string::npos){
		base_ouput_file = base_ouput_file.substr( 0, extension_index);
	}

	Rprintf("Building output folder...\n");
	// Test if the output dir exists
	struct stat sb;
	string output_folder = (string) base_output_folder;
	if( stat( output_folder.c_str(), &sb) != 0){
		string command = "mkdir -p " + (string) output_folder;
		system( command.c_str());
	}

	Rprintf( "\nSaving results...\n");
	// outputting results on sizes
	string output_file_path_class_size = output_folder + "/" + base_ouput_file + "_repeatClassCoverage.txt";
	Rprintf( "|-- Writing repeats class size stats to : %s\n", output_file_path_class_size.c_str());
	FILE *file_out_class_size = fopen( output_file_path_class_size.c_str(), "wb");
	for( map<string, long int*>::iterator itr_class = class_size->begin(); itr_class != class_size->end(); ++itr_class){
		string class_name_sign = itr_class->first;
		string class_name = class_name_sign.substr(0, class_name_sign.length()-1);
		string sign = class_name_sign.substr(class_name_sign.length()-1, 1);
		long int total;
		if( sign.compare( "+") == 0){
			total = total_wig_size_pos;
		}
		else{
			total = total_wig_size_neg;
		}
		fprintf( file_out_class_size, "%s\t%s\t%ld\t%ld\t%ld\n", class_name.c_str(), sign.c_str(), *(itr_class->second), total, total_wig_size_pos+total_wig_size_neg);
	}
	fclose( file_out_class_size);

	string output_file_path_family_size = output_folder + "/" + base_ouput_file + "_repeatFamilyCoverage.txt";
	Rprintf( "|-- Writing repeats family size stats to : %s\n", output_file_path_family_size.c_str());
	FILE *file_out_family_size = fopen( output_file_path_family_size.c_str(), "wb");
	for( map<string, long int*>::iterator itr_family = family_size->begin(); itr_family != family_size->end(); ++itr_family){
		string family_name_sign = itr_family->first;
		string family_name = family_name_sign.substr(0, family_name_sign.length()-1);
		string sign = family_name_sign.substr(family_name_sign.length()-1, 1);
		long int total;
		if( sign.compare( "+") == 0){
			total = total_wig_size_pos;
		}
		else{
			total = total_wig_size_neg;
		}
		fprintf( file_out_family_size, "%s\t%s\t%ld\t%ld\t%ld\n", family_name.c_str(), sign.c_str(), *(itr_family->second), total, total_wig_size_pos+total_wig_size_neg);
	}
	fclose( file_out_family_size);

	// outputting results on weights
	string output_file_path_class_weight = output_folder + "/" + base_ouput_file + "_repeatClassWeight.txt";
	Rprintf( "|-- Writing repeats class weight stats to : %s\n", output_file_path_class_weight.c_str());
	FILE *file_out_class_weight = fopen( output_file_path_class_weight.c_str(), "wb");
	for( map<string, long int*>::iterator itr_class = class_weight->begin(); itr_class != class_weight->end(); ++itr_class){
		string class_name_sign = itr_class->first;
		string class_name = class_name_sign.substr(0, class_name_sign.length()-1);
		string sign = class_name_sign.substr(class_name_sign.length()-1, 1);
		long int total;
		if( sign.compare( "+") == 0){
			total = total_wig_weight_pos;
		}
		else{
			total = total_wig_weight_neg;
		}
		fprintf( file_out_class_weight, "%s\t%s\t%ld\t%ld\t%ld\n", class_name.c_str(), sign.c_str(), *(itr_class->second), total, total_wig_weight_pos+total_wig_weight_neg);
	}
	fclose(file_out_class_weight);

	string output_file_path_family_weight = output_folder + "/" + base_ouput_file + "_repeatFamilyWeight.txt";
	Rprintf( "|-- Writing repeats family weight stats to : %s\n", output_file_path_family_weight.c_str());
	FILE *file_out_family_weight = fopen( output_file_path_family_weight.c_str(), "wb");
	for( map<string, long int*>::iterator itr_family = family_weight->begin(); itr_family != family_weight->end(); ++itr_family){
		string family_name_sign = itr_family->first;
		string family_name = family_name_sign.substr(0, family_name_sign.length()-1);
		string sign = family_name_sign.substr(family_name_sign.length()-1, 1);
		long int total;
		if( sign.compare( "+") == 0){
			total = total_wig_weight_pos;
		}
		else{

			total = total_wig_weight_neg;
		}
		fprintf( file_out_family_weight, "%s\t%s\t%ld\t%ld\t%ld\n", family_name.c_str(), sign.c_str(), *(itr_family->second), total, total_wig_weight_pos+total_wig_weight_neg);
	}
	fclose( file_out_family_weight);

}

/**
 *
 *
 */
void WigRepeatAnalyzer::accumulateValue( map<string, long int*>* _map, string key, long int value) {

	map<string, long int*>::iterator key_itr = _map->find( key);
	if( key_itr == _map->end()){
		_map->insert( pair<string, long int*>( key, new long int( value)));
	}
	else{
		*(key_itr->second) += value;
	}

}

/**
 * Read a line in the Fixed step WIG file and build a WigLocus instance that describe the locus
 *
 */
WigLocus* WigRepeatAnalyzer::readWigLocus( ifstream& wig_input_file, WigLocus* last_locus, bool replace) {

	if( wig_input_file.eof()){
		return NULL;
	}

	string line;
	try{

		if( last_locus == NULL){
			bool found = false;
			do{
				getline( wig_input_file, line);
				//Rprintf("1.line=%s\n",line.c_str());
				if( !line.empty()){
					size_t space_index = line.find( " ");
					if( space_index != string::npos){
						if( line.substr( 0, 9).compare( "fixedStep") == 0){
							found = true;
						}
					}
					else{
						found = true;
					}
				}
			} while( !wig_input_file.eof() && !found && !line.empty());
		}
		else{
			getline( wig_input_file, line);
			//Rprintf("2.line=%s\n",line.c_str());
		}

		if( line.empty()){
			return NULL;
		}
		char *token, str[1000];
		strcpy( str, line.c_str());

		string first_token = std::strtok( str, " ");

		//if( !replace){
		//	Rprintf( "\nChange WIG Locus" << endl;
		//}

		// If the line is a decription line, retrieve information on chromosom, start and step
		// Then recall readWiglocus to get the value from the next line
		if( first_token.compare( "fixedStep") == 0){
			// Read the chromosom name
			string chromosom_token = std::strtok( NULL, " ");
			string chromosom = "";
			size_t equal_index = chromosom_token.find( "=");
			if( equal_index != string::npos){
				chromosom = chromosom_token.substr( equal_index + 1);
			}

			// Read the start index
			string start_token = std::strtok( NULL, " ");
			int start = -1;
			equal_index = start_token.find( "=");
			if( equal_index != string::npos){
				start = atoi( start_token.substr( equal_index + 1).c_str());
			}

			// Read the step value
			string step_token = std::strtok( NULL, " ");
			int step = -1;
			equal_index = step_token.find( "=");
			if( equal_index != string::npos){
				step = atoi( step_token.substr( equal_index + 1).c_str());
			}
			//Rprintf( "|-- Wig locus definition : " << chromosom << " / " << start << " / " << step << endl;
			// Build the WigLocus from retrieved information
			WigLocus* locus = new WigLocus( chromosom, start, step, 0);

			// Set the value into the WingLocus by reading the next line of the file
			locus = readWigLocus( wig_input_file, locus, true);

			return locus;
		}
		// If the line contains a value, get it and build the WigLocus
		else{
			int value = atoi( first_token.c_str());
			// If 'replace' is true, it means the last WigLocus must be used and modified
			if( replace){
				last_locus->setValue( value);
				//Rprintf( "|-- Wig locus value : " << last_locus->getChromosom() << ":" << last_locus->getStart() << "-" << last_locus->getEnd() << " (" << last_locus->getStep() << ") " << last_locus->getValue() << endl;
				return last_locus;
			}
			// If 'replace' is false, it means a new WigLocus must be created from the information
			// stored into the last locus
			else{
				int step = last_locus->getStep();
				int start = last_locus->getStart() + step;
				WigLocus* locus = new WigLocus( last_locus->getChromosom(), start, step, value);
				//Rprintf( "|-- Wig locus value : " << locus->getChromosom() << ":" << locus->getStart() << "-" << locus->getEnd() << " (" << locus->getStep() << ") " << locus->getValue() << endl;
				return locus;
			}
		}
	}
	catch( exception& e){
		error( "ERROR : error reading wig line : %s", line.c_str());
		return NULL;
	}

}

/**
 * Read a locus in the RepeatMasker file
 *
 */
RepeatLocus* WigRepeatAnalyzer::readRepeatLocus( ifstream& repeatmasker_input_file) {

	if( repeatmasker_input_file.eof()){
		return NULL;
	}

	string line;
	getline( repeatmasker_input_file, line);
	//Rprintf("RepeatLocus Line = %s\n", line.c_str());

	if( line.empty()){
		return NULL;
	}

	char *token, str[200];
	strcpy( str, line.c_str());
	token = std::strtok( str, "\t");
	token = std::strtok( NULL, "\t");
	token = std::strtok( NULL, "\t");
	token = std::strtok( NULL, "\t");
	token = std::strtok( NULL, "\t");

	// Read the name of the chromosom
	token = std::strtok( NULL, "\t");
	string current_repeat_chrom;
	if( token != NULL){
		current_repeat_chrom = token;
	}
	else{
		error( "A line in the Repeat Masker file is incorrect : %s", line.c_str());
	}

	// Read the repeat locus start position
	token = std::strtok( NULL, "\t");
	int current_repeat_start;
	if( token != NULL){
		current_repeat_start = atoi( token);
	}
	else{
		error( "A line in the Repeat Masker file is incorrect : %s", line.c_str());
	}

	// Read the repeat locus end position
	token = std::strtok( NULL, "\t");
	int current_repeat_end;
	if( token != NULL){
		current_repeat_end = atoi( token);
	}
	else{
		error( "A line in the Repeat Masker file is incorrect : %s", line.c_str());
	}
	token = std::strtok( NULL, "\t");
	token = std::strtok( NULL, "\t");
	token = std::strtok( NULL, "\t");

	// Read the repeat locus class
	token = std::strtok( NULL, "\t");
	string current_repeat_class;
	if( token != NULL){
		current_repeat_class = token;
	}
	else{
		error( "A line in the Repeat Masker file is incorrect : %s", line.c_str());
	}

	// Read the repeat locus family
	token = std::strtok( NULL, "\t");
	string current_repeat_family;
	if( token != NULL){
		current_repeat_family = token;
	}
	else{
		error( "A line in the Repeat Masker file is incorrect : %s", line.c_str());
	}

	//Rprintf( "\nChange Repeat Locus" << endl;
	//Rprintf( "|-- Repeat definition: " << current_repeat_chrom << " / " << current_repeat_start << " / " << current_repeat_end << " / " << current_repeat_class << endl;
	// Build the RepeatLocus
	RepeatLocus* locus = new RepeatLocus( current_repeat_chrom, current_repeat_start, current_repeat_end, current_repeat_class, current_repeat_family);

	return locus;
}

/**
 * Help on software usage
 *
 */
void WigRepeatAnalyzer::usage() {

	Rprintf(
			"This tool analyzes the portion of signal provided in a WIG file that is located in repeat locii. Analysis reports signal coverage (the length of the loci) and signal weight (the length pondered by the signal intensity).\n");
	Rprintf( " The output is composed of 4 files, 2 for signal coverage (classified by repeat class and repeat family) and 2 for signal weight (classified by repeat class and repeat family).\n");
	Rprintf( " VERY IMPORTANT : The WIG file and the repeat masker file MUST have the same chromosom order.\n\n");
	Rprintf( " Usage: WigRepeatAnalyzer [--wigFolder|-d] <wig_folder> [--wigRegex|-i] <wig_regex> [--repeatMasker|-m] <repeat_masker> -- [outputFolder|--o] <out_folder>\n");
	Rprintf( "     - <wig_folder> : path to the WIG file to analyze.\n");
	Rprintf( "     - <wig_regex> : regular expression describing the WIG file to analyze. Can be its complete name.\n");
	Rprintf( "     - <repeat_masker> : path to the file containing repeats definitions.\n");
	Rprintf( "     - <out_folder> : path of the output folder where result files will be stored. If not specified, the wig_folder will be used.\n");
}

