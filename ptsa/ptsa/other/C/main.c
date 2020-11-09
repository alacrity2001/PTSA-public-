#include "VALMOD.h"
#include <stdio.h>
#include <stdlib.h>
mpElements * selfJoin(char * file, int lengthSeries, int lengthSubsequence , _Bool bPrintFileMP)
{
    //printf("Launch the selfjoin with %s %d %d\n",file,lengthSeries,lengthSubsequence );
    mpElements * mpEl = selfJoinSTOMP(file, lengthSeries, lengthSubsequence,bPrintFileMP);
    return mpEl;
}

mpElements * join(char * fileA, char * fileB, int lengthSeriesA, int lengthSeriesB , int lengthSubsequence, _Bool bPrintFileMP )
{
    //printf("Launch the join %s on %s,  lengths: %d, %d - length subsequence: %d \n", fileA, fileB, lengthSeriesA, lengthSeriesB, lengthSubsequence);
    mpElements * mpEl = joinSeriesSTOMP(fileA, fileB, lengthSeriesA, lengthSeriesB, lengthSubsequence,bPrintFileMP);
    return mpEl;
}

mpElements * freempElements( mpElements * mpEl )
{
	free(mpEl->matrixProfile);
	free(mpEl->indxProfile);
        free(mpEl->standardProfile);
	free(mpEl);
	
}