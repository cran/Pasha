#include "MultireadDispatcher.h"
#include "R.h"

 extern "C" {

 using namespace std;

 void C_UniformDispatch ( char** r_file_name, char** r_output_dir, char** r_reference, char** returnedOutputFile) {

	    // Test the passed arguments
	 	if( r_file_name == NULL || r_output_dir == NULL || r_reference == NULL ){
	 		error("Error: NULL argument passed to function.");
	 		return;
	 	}

	 	// Convert the R Vector to suitable types
	 	char* file_name = r_file_name[0];
	 	char* output_dir = r_output_dir[0];
	 	char* reference = r_reference[0];


		// Launch the analysis if the parameters are ok
		if( file_name != NULL){
			if( output_dir != NULL){
				if( reference != NULL){
					try{
						MultireadDispatcher analyzer = MultireadDispatcher( output_dir);
						analyzer.loadReference( reference);
						analyzer.analyzeTags( file_name, returnedOutputFile);
					}
					catch( exception &e){
						error( "Error : %s", e.what());
					}
				}
				else{
					error( "Error : No genome reference specified : %s", output_dir);
					MultireadDispatcher::usage();
				}
			}
			else{
				error( "Error : No output dir specified : %s", output_dir);
				MultireadDispatcher::usage();
			}
		}
		else{
			error( "Error : No input file specified : %s", file_name);
			MultireadDispatcher::usage();
		}

	}

 } // extern "C"
