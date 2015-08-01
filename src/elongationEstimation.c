#include <R.h>

// The function assumes that the piledValues vectors for the positive and negative strands have the same size
void C_elongationEstimation(int *piledValuesPlus, int *piledValuesMinus, double *stepShift, int *maxShift, int *piledValues_Size, int *stepShift_Size) 
{
	// Pointers shifted relatively from the begin of the negative strand vector
	int* sumPointers[(*stepShift_Size)];

	for(int j=0; j<(*stepShift_Size); j++)
	{
		// Positionning the pointers to the positions shifted as specified in "stepShift" (relative to the begin of the negative strand vector)
		sumPointers[j]=piledValuesMinus+((int)stepShift[j]);
		stepShift[j]=0;
	}

	// Now, stepShift will be used to store the resulting sums computed along the main loop for each corresponding pointer

	// We'll go over the vector, trying to avoid putting "negative strand shifted pointers" outside of the allocated space
	for(int i=0; i<((*piledValues_Size)-((*maxShift)+1)); i++)
	{
		for(int j=0; j<(*stepShift_Size); j++)
		{
			// incrementing all the "shifted" pointers for the negative strand and summing the multiplication result 
			// with the current "corresponding" position  at each point of the negative strand vector
			stepShift[j]+=((*(sumPointers[j]++))*((piledValuesPlus[i])));
		}
	}
}

