#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


/* Declare foreign language functions and register them in R */
void C_binVector(int *piledValues, double *result, int *binSize, int *piledValues_Size);
void C_elongationEstimation(int *piledValuesPlus, int *piledValuesMinus, double *stepShift, int *maxShift, int *piledValues_Size, int *stepShift_Size);
void C_pileupDouble(int * start, int * fraglength, int * dir,  int * readlength, double * weight, int * nbTags, int * maxCoord, double * res);
void C_RemoveArtifact (char** r_file_name, char** r_output_path, char** r_reference, int* r_ratio, int* r_verbosity, char** returnedOutputFile);
void C_CSEMDispatch ( char** r_file_name, char** r_output_dir, char** r_reference, int* r_window_size, int* r_iteration_number, char** returnedOutputFile);
void C_UniformDispatch ( char** r_file_name, char** r_output_dir, char** r_reference, char** returnedOutputFile);
void C_AnalyzeRepeat( char** r_file_name, char** r_input_folder, char** r_output_folder, char** r_repeat_masker_file_path);

void R_init_Pasha(DllInfo *info)
{
    /* Create the R_CMethodDef array */
    R_CMethodDef cMethods[] = {
        {"C_binVector", (DL_FUNC) &C_binVector, 4, (R_NativePrimitiveArgType[4]) {INTSXP, REALSXP, INTSXP, INTSXP}},
        {"C_elongationEstimation", (DL_FUNC) &C_elongationEstimation, 6, (R_NativePrimitiveArgType[6]) {INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, INTSXP}},
        {"C_pileupDouble", (DL_FUNC) &C_pileupDouble, 8, (R_NativePrimitiveArgType[8]) {INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP}},
        {"C_RemoveArtifact", (DL_FUNC) &C_RemoveArtifact, 6, (R_NativePrimitiveArgType[6]) {STRSXP, STRSXP, STRSXP, INTSXP, INTSXP, STRSXP}},
        {"C_CSEMDispatch", (DL_FUNC) &C_CSEMDispatch, 6, (R_NativePrimitiveArgType[6]) {STRSXP, STRSXP, STRSXP, INTSXP, INTSXP, STRSXP}},
        {"C_UniformDispatch", (DL_FUNC) &C_UniformDispatch, 4, (R_NativePrimitiveArgType[4]) {STRSXP, STRSXP, STRSXP, STRSXP}},
        {"C_AnalyzeRepeat", (DL_FUNC) &C_AnalyzeRepeat, 4, (R_NativePrimitiveArgType[4]) {STRSXP, STRSXP, STRSXP, STRSXP}},
        {NULL, NULL, 0}
    };

    /* Register the routine */
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}
