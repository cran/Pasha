#include "WigRepeatAnalyzer.h"
#include "R.h"

using namespace std;

 extern "C" {

	 void C_AnalyzeRepeat( char** r_file_name, char** r_input_folder, char** r_output_folder, char** r_repeat_masker_file_path
			 //, int* is_regex
			 ) {

		Rprintf("Checking variables...\n");
		if( r_file_name == NULL){
			error("Error: input file name is NULLn.");
			return;
		}
		if( r_input_folder == NULL){
					error("Error: input folder is NULL.");
					return;
		}
		if( r_output_folder == NULL){
					error("Error: output folder is NULL.");
					return;
		}
		if( r_repeat_masker_file_path == NULL){
					error("Error: repeat masker file is NULL.");
					return;
		}
//		if( is_regex == NULL){
//					error("Error: regex is NULL.");
//					return;
//		}

		char* file_name = NULL;
		char* regex_file_name = r_file_name[0];
		char* input_folder = r_input_folder[0];
		char* output_folder = r_output_folder[0];
		char* repeat_masker_file_path = r_repeat_masker_file_path[0];

		// test if an input folder has been specified
		if( input_folder == NULL){
			error( "Error : No WIG folder specified.");
			WigRepeatAnalyzer::usage();
			return;
		}

		// test if an input file name (regex) has been specified
		if( r_file_name == NULL){
			error( "Error : No WIG file name specified.");
			WigRepeatAnalyzer::usage();
			return;
		}

		// If the provided filename is delared as a regex (is_regex = 1),
		// search for the corresponding file
//		if( is_regex[0] == 1){
//			Rprintf("Looking for input file as regex...\n");
//
//			// Compile the file name regular expression (regex)
//			regex_t pattern;
//			int compile_regex = regcomp( &pattern, regex_file_name, REG_NOSUB | REG_EXTENDED);
//			if( compile_regex != 0){
//				error( "Error : unable to compile regular expression : %s", regex_file_name);
//				return;
//			}
//
//			// search in the input folder for a unique WIG file corresponding to the given regex
//			bool found_file = false;
//			struct dirent* ep;
//			DIR *dp;
//			dp = opendir( input_folder);
//
//			if( dp != NULL){
//				while( ep = readdir( dp)){
//					int compare = regexec( &pattern, ep->d_name, 0, NULL, 0);
//					if( compare == 0){
//						if( !found_file){
//							file_name = ep->d_name;
//							found_file = true;
//							Rprintf("  Input file found : %s\n", file_name);
//						}
//						else{
//							error( "Error : Several files found with the given regular expression '%s' : ", regex_file_name);
//							error( "   -> %s", file_name);
//							error( "   -> %s", ep->d_name);
//							error( "Please choose an other regular expression to have a single file identified.");
//							return;
//						}
//					}
//				}
//
//				(void) closedir( dp);
//			}
//			else{
//				error( "Error : Couldn't open the input folder : %s", input_folder);
//				return;
//			}
//			regfree( &pattern);
//		}
//		else{
			file_name = r_file_name[0];
			Rprintf("Considering WIG file:%s", file_name);
//		}

		// Check if the output folder is null. If so, use the input one as output
		if( output_folder == NULL){
			output_folder = input_folder;
		}

		if( file_name != NULL){
			if( repeat_masker_file_path != NULL){
				try{
					Rprintf("Launching analysis...");
					WigRepeatAnalyzer analyzer;
					analyzer.analyzeRepeat( input_folder, file_name, repeat_masker_file_path, output_folder);
					Rprintf("Finished.");
					return;
				}
				catch( exception &e){
					error( "ERROR : %s", e.what());
					return;
				}
			}
			else{
				error( "Error : No Repeat Masker file specified : %s", repeat_masker_file_path);
				WigRepeatAnalyzer::usage();
				return;
			}
		}
		else{
			error( "Error : No input WIG file found with name/regex : %s", r_file_name[0]);
			WigRepeatAnalyzer::usage();
			return;
		}

	 }

 }
