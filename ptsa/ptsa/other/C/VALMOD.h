#ifndef VALMOD_H
#define VALMOD_H

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
}
#endif

typedef struct MatrixProfileStuff
{
	double * matrixProfile;
	int * indxProfile;
        double * standardProfile;
	
} mpElements;

mpElements * joinSeriesSTOMP(char * fileTSA, char * fileTSB, int lengthSeriesA, int lengthSeriesB, int lengthSubsequence,_Bool bprintMP);
mpElements * selfJoinSTOMP(char * fileTSA, int lengthSeries, int lengthSubsequence,_Bool bprintMP);
#endif /* VALMOD_H */