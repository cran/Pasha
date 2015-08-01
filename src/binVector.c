#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// result has to be initialized before the call to a size of trunc(piledValues_Size/(double)binSize)+1
void C_binVector(int *piledValues, double *result, int *binSize, int *piledValues_Size) 
{
	// Saving the original adress of the vectors (to restore it at the end because it will be incremented along the loop)
	int * originalPointer_piledValues=piledValues;
	double * originalPointer_result=result;

	int currentSum=0;
	int i;

	// Looping on the values, incrementing the "i" from 1 to binSize in the loop
    for(i=1; piledValues < (originalPointer_piledValues+(*piledValues_Size)); i=(i<(*binSize))?(i+1):1)
	{
		// Summing the values and incrementing the pointer of values
		currentSum+=(*(piledValues++));

		//----Rprintf("Value  -> %d\n", *piledValues);
		//----Rprintf("CurSum -> %d\n\n", currentSum);

		// The end of the current bin is reached
		if(i==(*binSize))
		{
			//----Rprintf("End of bin\n");

			// Computing the average for the current bin and put it in the result vector
			(*result)=currentSum/(double)(*binSize);
			// Going to the next element of result
			result++;
			// Reinitialize the counting variable
			currentSum=0;
		}
	}

	//----Rprintf("AfterLoop i -> %d\n\n", i);

	// Piling the last bin (which is probably necessary because it's unlikely to have exactly the good amount of values to fit the binSize)
	// The test is testing equality to one because when we are after the last element piledValues in the loop (getting out of the loop)
	// the incrementation of i correspondiong to the last element has been done already. So if synchronous, the i would be reinitialized to 1.
	if(i!=1)
	{
		//----Rprintf("Treating last bin apart...\n");
		(*result)=currentSum/(double)(i-1);
	}

	// Putting back the pointers to the beginning of the vectors
	piledValues=originalPointer_piledValues;
	result=originalPointer_result;

}
