#include "MultireadArtifactRemover.h"
#include "R.h"

 extern "C" {

 using namespace std;

 void C_RemoveArtifact (char** r_file_name, char** r_output_path, char** r_reference, int* r_ratio, int* r_verbosity, char** returnedOutputFile) {

	    // Test the passed arguments
	 	if( r_file_name == NULL || r_reference == NULL || r_ratio == NULL || r_verbosity == NULL || r_output_path == NULL){
	 		error("Error: NULL argument passed to function.");
	 		return;
	 	}

	 	// Convert the R Vector to suitable types
	 	char* file_name = r_file_name[0];
	 	char* reference = r_reference[0];
	 	char* output_path = r_output_path[0];
	 	int ratio = r_ratio[0];
	 	int verbosity = r_verbosity[0];

		// Launch the analysis if the parameters are ok
		if( file_name != NULL){
			if( ratio > 0){
				if( reference != NULL){
					try{
						ArtifactRemover analyzer = ArtifactRemover( verbosity, file_name, output_path);
						analyzer.loadReference( reference);
						analyzer.loadData( file_name);
						analyzer.removeArtifact( ratio);
						analyzer.outputResult( returnedOutputFile);
					}
					catch( exception &e){
						error("Error : %s",e.what());
					}
				}
				else{
					error( "Error : No genome reference directory specified.");
					ArtifactRemover::usage();
				}
			}
			else{
				error( "Error : No artifact detection ratio specified : %d", ratio);
				ArtifactRemover::usage();
			}
		}
		else{
			error( "Error : No input file specified : %s", file_name);
			ArtifactRemover::usage();
		}

	}

 } // extern "C"
