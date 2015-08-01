#include <R.h>

// All arguments (except maxCoord, nbTags and res) must have a length==nbTags, res must have a length==maxCoord
void C_pileupDouble(int * start, int * fraglength, int * dir,  int * readlength, double * weight, int * nbTags, int * maxCoord, double * res)
{

	int i, j, st, end;

	for(i = (*nbTags) ; i-- ; )
	{

		//Rprintf( "%d", i);
		if( dir[i] == 1 )
		{
			// forward direction ASSUMING STRAND IS REPRESENTED LIKE IN IRANGES AND LAST VERSIONS OF SHORTREAD AS '+ -' and eventually other ones after (which will be ignored)
			end = start[i] + fraglength[i]-1;
			//if( end > (*maxCoord) ) Rprintf( "Fragment is going over the end of the chromosome (fw)\n" );

			for( j = start[i]; j <= end; j++ )
			{
				// Don't write outside of the array
				if(j>0 && j<=(*maxCoord))
				res[j-1] += weight[i]; // because 1-based coordinates in a 0-based array
			}
		}
		else if(dir[i] == 2)
		{
			// reverse direction
			st = start[i] + readlength[i]-1;
			//if( st > (*maxCoord) ) Rprintf( "Fragment is going over the end of the chromosome (bw)\n" );

			end = st - fraglength[i];
			//if( end < 0 ) Rprintf( "Tag extension starting before the chromosome (bw)\n" );

			for( j = st; j > end; j-- )
			{
				// Don't write outside of the array
				if(j>0 && j<=(*maxCoord))
				res[j-1] += weight[i];
			}
		}

	}
}

