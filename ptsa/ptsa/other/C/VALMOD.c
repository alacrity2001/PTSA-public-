
#include <fftw3.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#define INF 1e20  
#include <string.h>
#include <stdbool.h>
#include "FiboKLAMD.h"
#include "VALMOD.h"
//###############################--Functions and structures definition--###############################
//#define SANITY_CHECK 1


typedef struct IndexDistanceProfile{
	double bestDist;
	int bestIdx;
	double * distProfile;
	int * indxProfile;
	int length;
        double standardDev;
} idProf;



typedef struct RangeMatrixProfileStuff
{
	double ** matrixProfiles;
	int ** indexProfiles;
        double **  standardProfiles;
} rangeMpElements;

typedef struct VariableLengthMatrixProfile
{
	double * matrixProfile;
	int * indexProfile;
	int * lengthProfile;
	double * distancesEvolutionMotif;
        
        double * matrixProfileNonLengthNormalized;
        int * indexProfileNonLengthNormalized;
        int * lengthProfileNonLengthNormalized;
        
} VALMAP;




idProf orderedidP; // global for passing the output of Distance profile among functions

double randPersonal (double mu, double sigma);
int reverseSignal(double *in, int n);
double * computeMaxForLB(double * in, int sizeIn, int minSize, int maxSize);
double * slidingDotProduct(double * longSeries, double * query, int sizeSeries, int sizeQuery);
double * computeMeanStd(double * series, int sizeSeries, int sizeQuery);
double * computeSumSumSquared(double * series, int sizeSeries, int sizeQuery);
idProf returnBestDistanceOfARow(double * QT, double * meandStdQuery, double * meanStdSeries, int indxQ, int sizeSeries,int sizeQuery, _Bool bSelfJoin);
idProf returnDistanceProfile(double * QT, double * meandStdQuery, double * meanStdSeries, int indxQ, int sizeSeries,int sizeQuery, _Bool bSelfJoin);
mpElements STAMP(double * Ta, double * Tb, int sizeQuery, int sizeTa, int sizeTb, _Bool selfJoin, const char * fileNameMP);
heapEntry STOMP(rangeMpElements rMP, int offset, double * Ta, double * Tb, int sizeQuery, int sizeTa, int sizeTb, _Bool selfJoin, const char * fileNameMP);
int selfJoinWith_STAMP_and_STOMP();
int basicExperiment();
//###################################################################################################

//Functions, Variables, Data structures realtive to the KLAMD_MOENIZED algorithm
//***************************************************************************************************
void returnDistanceProfileWithMotifsProfileComputation(double * QT, double * meandStdQuery, double * meanStdSeries, int indxQ, int sizeSeries,int sizeQuery, _Bool bSelfJoin, mMxHeap * motifProfile, double * sumSumSQData, double * sumSumSQQuery);
double computeTrueDistanceEntryHeapEntry(double * Ta, double * Tb, int sizeTa, int sizeTb, heapEntry * entry, int actualQuerySize);
//***************************************************************************************************


//VALMOD is the algorithm which lower bounds the Distance Profiles!! It may use the Yan bound or the MOEN bound
void returnDistanceProfileWithMotifsProfileComputation_VALMOD(double * QT, double * meandStdQuery, double * meanStdSeries, int indxQ, int sizeSeries,int sizeQuery, _Bool bSelfJoin, mMxHeap ** distancesProfilesTop, double * sumSumSQData, double * sumSumSQQuery);
heapEntry STOMP_MotifsProfile_VALMOD(VALMAP * multiMP, double * QT, int offset, double * Ta, double * Tb, int sizeQuery, int sizeTa, int sizeTb, _Bool selfJoin, const char * fileNameMP, mMxHeap ** listMotifsProfile);
void VALMOD(VALMAP * multiMP, double * Ta, double * Tb, int sizeQuery, int sizeQueryMAX, int sizeTa, int sizeTb, _Bool selfJoin, const char * fileNameMP);
void VALMOD_RUN(int argc, char ** argv);
int selfJoinWith_STOMPmultilength_VALMOD(int sizeSeries,int min,int max,char * fileTS);
///////////////////////////////////////////////////////////////////////////////////////


_Bool sanityCheckLengths(int lengthSeriesA, int lengthSeriesB,int lengthSubsequence);
_Bool sanityCheckFile(char * nfileTS);
_Bool sanityCheckLengthsSelfJoin(int lengthSeriesA,int lengthSubsequence);


void JOIN_RUN(int argc, char ** argv);
void SELFJOIN_RUN(int argc, char ** argv);
mpElements * _join_STOMP(int sizeSeries, int sizeSeriesB, int length, char * fileTSA, char * fileTSB, _Bool bPrintMP);
mpElements * _selfjoin_STOMP(int sizeSeries, int length, char * fileTSA, _Bool bPrintMP);

int threshold ; // VALMOD: elements of the distance profile kept in memory
heapEntry discord;
///////////
_Bool bComparisonStomp;
_Bool bOnlySTOMP;
char * fileNameValmap;
 int kDiscords;
 
 
double randPersonal (double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
 
  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }
 
  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
 
  call = !call;
 
  return (mu + sigma * (double) X1);
}


int reverseSignal(double *in, int n)
{
	int switchL = ceil(n/2);
	int i;
	for(i=0;i<switchL;i++)
	{
		double sw = in[i];
		in[i] = in[(n-1)-i];
		in[(n-1)-i] = sw;
	}
}

double * slidingDotProduct(double * longSeries, double * query, int sizeSeries, int sizeQuery)
{
	// long series and query should have allocated the same number of bytes, which is the bytes for longSeries * 2
	
	
	fftw_complex *out; // long series in the frequency domain 
	fftw_complex *out2; // query in the frequency domain
	fftw_complex *out3; // dot product in frequency domain
	fftw_plan plan_backward;
	fftw_plan plan_forward;
	fftw_plan plan_forward2;
	reverseSignal(query,sizeQuery); 
	int i;
	
	for (i = sizeSeries; i < (sizeSeries*2); i++ )
	{
		longSeries[i] = 0;
	}
	for (i = sizeQuery; i < (sizeSeries*2); i++ )
	{
		query[i] = 0;
	}
	
	int nc = sizeSeries + 1; // size of the long series and query in frequency domain

	out = fftw_malloc ( sizeof ( fftw_complex ) * nc );
	out2 = fftw_malloc ( sizeof ( fftw_complex ) * nc );
	out3 = fftw_malloc ( sizeof ( fftw_complex ) * nc );
  
	plan_forward = fftw_plan_dft_r2c_1d ( (sizeSeries*2), longSeries, out, FFTW_ESTIMATE );
	plan_forward2 = fftw_plan_dft_r2c_1d ( (sizeSeries*2), query, out2, FFTW_ESTIMATE );
  
	fftw_execute ( plan_forward);
	fftw_execute ( plan_forward2);
  
	// compute the dot product in frequency domain

	for (i = 0; i < nc; i++ )
	{
		out3[i][0] = ((out[i][0] * out2[i][0]) - (out[i][1] * out2[i][1]));
		out3[i][1] = ((out[i][0] * out2[i][1]) + (out[i][1] * out2[i][0])); 
	}
	
	double * finalDotProductSum = fftw_malloc ( sizeof ( double ) * (2*sizeSeries) );
	plan_backward = fftw_plan_dft_c2r_1d ( (2*sizeSeries), out3, finalDotProductSum, FFTW_ESTIMATE );
	fftw_execute ( plan_backward );

	
	double * outputDP =  fftw_malloc (sizeof(double)*(sizeSeries-sizeQuery+1));
	
	int posI=0;
	for (i=(sizeSeries-sizeQuery)+1;i>0;i--)
	{
		double dp = (finalDotProductSum[sizeSeries-i] / (double) (sizeSeries*2)) ;
		//printf ( "%12f \n", dp);
		outputDP[posI] = dp;
		posI++;
	}

	fftw_destroy_plan ( plan_forward );
	fftw_destroy_plan ( plan_forward2);
	fftw_destroy_plan ( plan_backward );
	fftw_free ( out );
	fftw_free ( out2);
	fftw_free ( out3 );
	fftw_free ( finalDotProductSum );

	return outputDP;
}

double * computeMeanStd(double * series, int sizeSeries, int sizeQuery)
{
	double * MeanStdComplete = fftw_malloc(sizeof(double) * ((sizeSeries-sizeQuery+1)*2));
	int j;
	int cont=0;
	double sum=0;
	double sumSq=0;
	for(j=0;j<sizeSeries;j++)
	{
		sum=sum+series[j];
		sumSq=sumSq+series[j]*series[j];
		if((j+1)>=sizeQuery)
		{ 
			double mean = sum/sizeQuery;
			MeanStdComplete[cont] = mean;
			double sigmaSq = ((sumSq/sizeQuery)-mean*mean);
			if (sigmaSq>0)
			{				
				MeanStdComplete[cont+(sizeSeries-sizeQuery+1)] = sqrt(sigmaSq);
			}
			else
			{
				MeanStdComplete[cont+(sizeSeries-sizeQuery+1)] = 0;
			}
			sum=sum-series[j-(sizeQuery-1)];
			sumSq=sumSq-series[j-(sizeQuery-1)]*series[j-(sizeQuery-1)];
			cont++;
		}
	}
	return MeanStdComplete;
}

double * computeSumSumSquared(double * series, int sizeSeries, int sizeQuery)
{
	double * SumSumSquared = fftw_malloc(sizeof(double) * ((sizeSeries-sizeQuery+1)*2));
	int j;
	int cont=0;
	double sum=0;
	double sumSq=0;
	for(j=0;j<sizeSeries;j++)
	{
		sum=sum+series[j];
		sumSq=sumSq+series[j]*series[j];
		if((j+1)>=sizeQuery)
		{ 
			SumSumSquared[cont] = sum;
			SumSumSquared[cont+(sizeSeries-sizeQuery+1)] = sumSq;
			sum=sum-series[j-(sizeQuery-1)];
			sumSq=sumSq-series[j-(sizeQuery-1)]*series[j-(sizeQuery-1)];
			cont++;
		}
	}
	return SumSumSquared;
}

void updateValmap(VALMAP * matrixProfileVL, double realDistance, double actualSizeQueryNext, int indexUpdateVM, int index2)
{
	//UPDATE VALMAP!
	double ln1 = 1.0/(double)actualSizeQueryNext;
	double ln1SQ = sqrt(ln1);
	double normalizedDistance = realDistance * ln1SQ;
	
	if(normalizedDistance < matrixProfileVL->matrixProfile[indexUpdateVM] )
	{
		matrixProfileVL->matrixProfile[indexUpdateVM] = normalizedDistance;
		matrixProfileVL->indexProfile[indexUpdateVM] = index2;
		matrixProfileVL->lengthProfile[indexUpdateVM] = actualSizeQueryNext;
	}
	if(realDistance < matrixProfileVL->matrixProfileNonLengthNormalized[indexUpdateVM] )
	{
		matrixProfileVL->matrixProfileNonLengthNormalized[indexUpdateVM] = realDistance;
		matrixProfileVL->indexProfileNonLengthNormalized[indexUpdateVM] = index2;
		matrixProfileVL->lengthProfileNonLengthNormalized[indexUpdateVM] = actualSizeQueryNext;
	}
}



idProf returnBestDistanceOfARow(double * QT, double * meandStdQuery, double * meanStdSeries, int indxQ, int sizeSeries,int sizeQuery, _Bool bSelfJoin)
{
	// TO DO return the complete Distance and index profile
	
	idProf idP;
	
	double * distanceProfile = (double *) fftw_malloc((sizeSeries-sizeQuery+1)*sizeof(double));
	int * indexProfile = (int *) fftw_malloc((sizeSeries-sizeQuery+1)*sizeof(int));

	double BestDist= INF;
	int besIndxData=-1;
	int i;
	for(i=0;i<(sizeSeries-sizeQuery+1);i++)
	{
		if(indxQ<(i-(sizeQuery/2)) || indxQ>(i+(sizeQuery/2)) || !bSelfJoin) // avoid the trivial matches
		{
			double dist= (2*sizeQuery) * (1- ( (QT[i] - (sizeQuery*meandStdQuery[0]*meanStdSeries[i])) / (sizeQuery*meandStdQuery[1]*meanStdSeries[i+(sizeSeries-sizeQuery+1)]) )) ;
			distanceProfile[i] = dist;
			indexProfile [i] = i;
			if(dist < BestDist)
			{
				besIndxData = i;
				BestDist = dist;
			}
		}
		else
		{
			indexProfile [i] = -1;
		}
	}

	//printf("entry matrix profile: query id: %d , sub-series id: %d , distance: %lf\n",indxQ,besIndxData,finalDist);
	
	idP.bestDist = BestDist;
	idP.bestIdx = besIndxData;
	idP.distProfile = distanceProfile;
	idP.indxProfile = indexProfile;
	
	return idP;
}

idProf returnDistanceProfile(double * QT, double * meandStdQuery, double * meanStdSeries, int indxQ, int sizeSeries,int sizeQuery, _Bool bSelfJoin)
{
	//return a void vector with the best distance and the best index at the first two positions and then the Distance profile
	idProf idP;
	double * distanceProfile = (double *) fftw_malloc((sizeSeries-sizeQuery+1)*sizeof(double));
	int * indexProfile = (int *) fftw_malloc((sizeSeries-sizeQuery+1)*sizeof(int));
	double BestDist= INF;
	int besIndxData=-1;
	int i;
	int contD =0;
        double sumDist = 0;
        double sumDistSq=  0;
        
	for(i=0;i<(sizeSeries-sizeQuery+1);i++)
	{
		if(indxQ<(i-(sizeQuery/2)) || indxQ>(i+(sizeQuery/2)) || !bSelfJoin) // avoid the trivial matches
		{
			double dist= (2*sizeQuery) * (1- ( (QT[i] - (sizeQuery*meandStdQuery[0]*meanStdSeries[i])) / (sizeQuery*meandStdQuery[1]*meanStdSeries[i+(sizeSeries-sizeQuery+1)]) )) ;
			if (dist<0)
			{
                            dist =0;
			}
    
                        sumDist = sumDist+dist;
                        sumDistSq = sumDistSq+(dist*dist);
                        
			distanceProfile[contD]=dist;
			indexProfile[contD]=i;
			contD++;
			if(dist < BestDist)
			{
				besIndxData = i;
				BestDist = dist;
			}
		}
	}
        
        idP.standardDev = 0;
        if(contD>0)
        {
            idP.standardDev = sqrt((sumDistSq/contD) - (sumDist/contD)*(sumDist/contD));
	}
        idP.bestDist = BestDist;
	idP.bestIdx = besIndxData;
	idP.distProfile = distanceProfile;
	idP.indxProfile = indexProfile;
	idP.length = contD;
	
	return idP;
}

void returnDistanceProfileWithMotifsProfileComputation(double * QT, double * meandStdQuery, double * meanStdSeries, int indxQ, int sizeSeries,int sizeQuery, _Bool bSelfJoin, mMxHeap * motifProfile, double * sumSumSQData, double * sumSumSQQuery)
{
	double BestDist= INF;
	int besIndxData=-1;
	int i;
	int contD=0;
	int limit = (sizeSeries-sizeQuery+1);
	for(i=0;i<limit;i++)
	{
		if(indxQ<(i-(sizeQuery/2)) || indxQ > (i+(sizeQuery/2)) || !bSelfJoin) // avoid the trivial matches
		{
			double stdDevData = meanStdSeries[i+(sizeSeries-sizeQuery+1)];
			double stdDevQuery = meandStdQuery[1];
			double dist= (2*sizeQuery) * (1- ( (QT[i] - (sizeQuery*meandStdQuery[0]*meanStdSeries[i])) / (sizeQuery*stdDevQuery*stdDevData) )) ;
			
			if(dist < BestDist)
			{
				besIndxData = i;
				BestDist = dist;
			}
			contD++;
			if(motifProfile!=NULL)
			{
				push(dist, indxQ, i, QT[i], sumSumSQQuery[0], sumSumSQData[i], sumSumSQQuery[1],  sumSumSQData[i+(sizeSeries-sizeQuery+1)], stdDevQuery, stdDevData, sizeQuery, motifProfile);
			}
		}
	}

	orderedidP.bestDist = BestDist;
	orderedidP.bestIdx = besIndxData;
	orderedidP.length = contD;
}


void returnDistanceProfileWithMotifsProfileComputation_VALMOD(double * QT, double * meandStdQuery, double * meanStdSeries, int indxQ, int sizeSeries,int sizeQuery, _Bool bSelfJoin, mMxHeap ** distancesProfilesTop, double * sumSumSQData, double * sumSumSQQuery)
{
	double BestDist= INF;
	int besIndxData=-1;
	int i;
	int contD=0;
	int limit = (sizeSeries-sizeQuery+1);
	for(i=0;i<limit;i++)
	{
		if(indxQ<(i-(sizeQuery/2)) || indxQ > (i+(sizeQuery/2)) || !bSelfJoin) // avoid the trivial matches
		{
			double stdDevData = meanStdSeries[i+(sizeSeries-sizeQuery+1)];
			double stdDevQuery = meandStdQuery[1];
			double dist= (2*sizeQuery) * (1- ( (QT[i] - (sizeQuery*meandStdQuery[0]*meanStdSeries[i])) / (sizeQuery*stdDevQuery*stdDevData) )) ;
			
			if(dist < BestDist)
			{
				besIndxData = i;
				BestDist = dist;
			}
			contD++;
			if(distancesProfilesTop!=NULL)
			{
				push(dist, indxQ, i, QT[i], sumSumSQQuery[0], sumSumSQData[i], sumSumSQQuery[1],  sumSumSQData[i+(sizeSeries-sizeQuery+1)], stdDevQuery, stdDevData, sizeQuery, distancesProfilesTop[indxQ]);
			}
		}
	}
        
	orderedidP.bestDist = BestDist;
	orderedidP.bestIdx = besIndxData;
	orderedidP.length = contD;
}



mpElements STAMP(double * Ta, double * Tb, int sizeQuery, int sizeTa, int sizeTb, _Bool selfJoin, const char * fileNameMP)
{
	
	int limit = (sizeTb-sizeQuery)+1;
	double * matrixProfile = fftw_malloc(sizeof(double)*limit);
	int * indexProfile = fftw_malloc(sizeof(int)*limit);
	double * copyTb = fftw_malloc(sizeof(double)*(sizeTb*2));
	double * TaMeanStdDev = computeMeanStd(Ta, sizeTa, sizeQuery);
	int i;
	FILE *  fMPfile = NULL;
	if(fileNameMP!=NULL){fMPfile = fopen(fileNameMP, "a");}
	
	for(i=0;i<limit;i++)
	{
		memcpy(copyTb, Tb, sizeQuery*(sizeof(double)));
		// copy Tb is the correct query 
		//##################################################--MASS--###############################################################
		double * QT = slidingDotProduct(Ta, copyTb, sizeTa, sizeQuery);
		double * queryMeanStd = computeMeanStd(Tb, sizeQuery, sizeQuery);
		
		idProf dip = returnBestDistanceOfARow( QT, queryMeanStd, TaMeanStdDev, i, sizeTa,sizeQuery,selfJoin); // store the best in the file
		//#########################################################################################################################
		
		matrixProfile[i] = dip.bestDist;
		indexProfile[i] = dip.bestIdx;
		
		if(fMPfile!=NULL){fprintf(fMPfile,"%lf\n",matrixProfile[i]);}
		
		fftw_free(dip.distProfile); // free distance profile since it is not needed
		fftw_free(dip.indxProfile); // free index profile since it is not needed
		fftw_free(queryMeanStd);
		fftw_free(QT);
		Tb++;
	}
	
	fftw_free(copyTb);
	fftw_free(TaMeanStdDev);
	Tb= Tb-(limit*sizeof(double));
	if(fMPfile!=NULL){fclose(fMPfile);}
	
	mpElements output;
	
	output.matrixProfile = matrixProfile;
	output.indxProfile = indexProfile;

	return output;
}

heapEntry STOMP(rangeMpElements rMP, int offset, double * Ta, double * Tb, int sizeQuery, int sizeTa, int sizeTb, _Bool selfJoin, const char * fileNameMP)
{
	FILE *  fMPfile = NULL;
	if(fileNameMP!=NULL){fMPfile = fopen(fileNameMP, "w");}
	heapEntry bestMP;

        int limitA = (sizeTa-sizeQuery)+1;
	int limit = (sizeTb-sizeQuery)+1;

	double * copyTb = fftw_malloc(sizeof(double)*2*sizeTa);
	memcpy(copyTb, Tb, sizeQuery*sizeof(double)); // get the first query of B
        // copyTb is the first query of B
	//##################################################--MASS--###############################################################
	double * QT = slidingDotProduct(Ta, copyTb, sizeTa, sizeQuery);
	double * queryMeanStd = computeMeanStd(Tb, sizeQuery, sizeQuery);
	double * TaMeanStdDev = computeMeanStd(Ta, sizeTa, sizeQuery);
	idProf dip = returnDistanceProfile(QT, queryMeanStd, TaMeanStdDev,0,sizeTa,sizeQuery,selfJoin);
	rMP.matrixProfiles[offset][0] = sqrt(dip.bestDist);
	rMP.indexProfiles[offset][0] = dip.bestIdx;
        rMP.standardProfiles[offset][0] = dip.standardDev;
	bestMP.distance = dip.bestDist;
	bestMP.index1 = 0;
	bestMP.index2 = dip.bestIdx;
	discord.distance = dip.bestDist;
	discord.index1 = 0;
	discord.index2 = dip.bestIdx;
        
	if(fMPfile!=NULL){fprintf(fMPfile,"%lf\n",rMP.matrixProfiles[offset][0]);}
	fftw_free(dip.distProfile); // free distance profile since it is not needed
	fftw_free(dip.indxProfile); // free index profile since it is not needed
	
	//#########################################################################################################################	
	double * copyTa = fftw_malloc(sizeof(double)*(sizeTb*2));
	memcpy(copyTa, Ta, sizeQuery*(sizeof(double))); // get the first query of A
	// copyTa is the first query  of a
	//##################################################--MASS--###############################################################
	double * QT_b = slidingDotProduct(Tb, copyTa, sizeTb, sizeQuery);
	//#########################################################################################################################
	
	int i,j;

	for(i=1;i<limit;i++)
	{
		for(j=(limitA-1);j>=1;j--)
		{
			QT[j] = QT[j-1] - (Ta[j-1] * Tb[i-1]) + (Ta[j+sizeQuery-1] * Tb[i+sizeQuery-1]);
 		}
		
		QT[0] = QT_b[i];
		double * queryMeanStdNext = computeMeanStd(&Tb[i], sizeQuery, sizeQuery);
		idProf dip = returnDistanceProfile(QT, queryMeanStdNext, TaMeanStdDev, i, sizeTa, sizeQuery, selfJoin);
		
		if(bestMP.distance>dip.bestDist)
		{
			bestMP.distance = dip.bestDist;
			bestMP.index1 = i;
			bestMP.index2 = dip.bestIdx;
		}
		else if(discord.distance<dip.bestDist)
		{
			discord.distance = dip.bestDist;
			discord.index1 = i;
			discord.index2 = dip.bestIdx;
		}
		
		rMP.matrixProfiles[offset][i] = sqrt(dip.bestDist);
		rMP.indexProfiles[offset][i] = dip.bestIdx;
                rMP.standardProfiles[offset][i] = dip.standardDev;
		fftw_free(dip.distProfile); // free distance profile since it is not needed
		fftw_free(dip.indxProfile); // free index profile since it is not needed
		fftw_free(queryMeanStdNext);
		if(fMPfile!=NULL){fprintf(fMPfile,"%lf\n",rMP.matrixProfiles[offset][i]);}
	}

	fftw_free(QT);
	fftw_free(QT_b);
	fftw_free(copyTa);
	fftw_free(copyTb);
	fftw_free(queryMeanStd);
	fftw_free(TaMeanStdDev);
	if(fMPfile!=NULL){fclose(fMPfile);}
	return bestMP;
}

double computeTrueDistanceEntryHeapEntry(double * Ta, double * Tb, int sizeTa, int sizeTb, heapEntry * entry, int actualQuerySize)
{

	int newIndexQuery = entry->index1+(actualQuerySize-1);
	int newIndexData = entry->index2+(actualQuerySize-1);
	
	if(newIndexQuery>=sizeTb ||  newIndexData >= sizeTa )
	{
		return -1.0;
	}
	
	entry->dP = entry->dP + (Tb[newIndexQuery]*Ta[newIndexData]);
	entry->s1 = entry->s1 + Tb[newIndexQuery];
	entry->s2 = entry->s2 + Ta[newIndexData];
	entry->ss1 = entry->ss1+(Tb[newIndexQuery]*Tb[newIndexQuery]);
	entry->ss2 = entry->ss2+(Ta[newIndexData]*Ta[newIndexData]);;
	
	double dotProduct = entry->dP;
	double avgQuery,stdQuery;
	avgQuery = entry->s1/actualQuerySize;
	stdQuery = (entry->ss1/actualQuerySize) - (avgQuery*avgQuery);
	
	double avgData,stdDevData,dist;
	avgData =  entry->s2/actualQuerySize;
	stdDevData = ( entry->ss2/actualQuerySize) - (avgData*avgData);
	
	if(stdDevData>0)
		stdDevData = sqrt(stdDevData);
	else
		stdDevData = 0;
	
	if(stdQuery>0)
		stdQuery = sqrt(stdQuery);
	else
		stdQuery = 0;
	
	dist = ((2*actualQuerySize) * (1- ( (dotProduct - (actualQuerySize*avgQuery*avgData)) / (actualQuerySize*stdQuery*stdDevData) ))) ;
	entry->distance = dist;

	return dist;				
}






heapEntry STOMP_MotifsProfile_VALMOD(VALMAP * matrixProfileVL, double * QT, int offset, double * Ta, double * Tb, int sizeQuery, int sizeTa, int sizeTb, _Bool selfJoin, const char * fileNameMP, mMxHeap ** listMotifsProfile)
{
	heapEntry bestMP;
	double * TaMeanStdDev = computeMeanStd(Ta, sizeTa, sizeQuery);	
	int limit = (sizeTb-sizeQuery)+1;
	double * QTChanging = malloc(sizeof(double) * limit);
	//##################################################--MASS--###############################################################
	double * queryMeanStd = computeMeanStd(Tb, sizeQuery, sizeQuery);
	double * minSizeSumSumSq = computeSumSumSquared(Ta, sizeTa, sizeQuery);
	double * querySumSumSq = computeSumSumSquared(Tb, sizeQuery, sizeQuery);
	returnDistanceProfileWithMotifsProfileComputation_VALMOD(QT, queryMeanStd, TaMeanStdDev, 0,sizeTa,sizeQuery, selfJoin, listMotifsProfile, minSizeSumSumSq, querySumSumSq);
	fftw_free(querySumSumSq);
	
	updateValmap(matrixProfileVL, sqrt(orderedidP.bestDist), sizeQuery, 0, orderedidP.bestIdx);
		   
	//multiMP->matrixProfiles[offset][0] = orderedidP.bestDist; // orderedidP is global 
	//multiMP->indexProfiles[offset][0] = orderedidP.bestIdx;

	bestMP.distance = orderedidP.bestDist;
	bestMP.index1 = 0;
	bestMP.index2 = orderedidP.bestIdx;
	//#########################################################################################################################	
	double * copyTa = fftw_malloc(sizeof(double)*(sizeTa*2));
	memcpy(copyTa, Ta, sizeQuery*(sizeof(double))); // get the first query of A
	// copyTa is the first query  of A
	//##################################################--MASS--###############################################################
	double * QT_b = slidingDotProduct(Tb, copyTa, sizeTb, sizeQuery);
	//#########################################################################################################################
	int i,j;
	for(i=1;i<limit;i++)
	{
		if(i==1)
		{
			for(j=(limit-1);j>=1;j--)
			{
				QTChanging[j] = QT[j-1] - Ta[j-1] * Tb[i-1] + Ta[j+sizeQuery-1] * Tb[i+sizeQuery-1];
			}
		}
		else
		{
			for(j=(limit-1);j>=1;j--)
			{
				QTChanging[j] = QTChanging[j-1] - Ta[j-1] * Tb[i-1] + Ta[j+sizeQuery-1] * Tb[i+sizeQuery-1];
			}
		}
		QTChanging[0] = QT_b[i];
			
		double * queryMeanStdNext = computeMeanStd(&Tb[i], sizeQuery, sizeQuery);
		querySumSumSq = computeSumSumSquared(&Tb[i], sizeQuery, sizeQuery);
		returnDistanceProfileWithMotifsProfileComputation_VALMOD(QTChanging, queryMeanStdNext, TaMeanStdDev, i,sizeTa,sizeQuery, selfJoin, listMotifsProfile,minSizeSumSumSq, querySumSumSq);
		fftw_free(querySumSumSq);
		
		updateValmap(matrixProfileVL, sqrt(orderedidP.bestDist), sizeQuery, i, orderedidP.bestIdx);
		
		
		//multiMP->matrixProfiles[offset][i] = orderedidP.bestDist;
		//multiMP->indexProfiles[offset][i] = orderedidP.bestIdx;
		if(bestMP.distance>orderedidP.bestDist)
		{
			bestMP.distance = orderedidP.bestDist;
			bestMP.index1 = i;
			bestMP.index2 = orderedidP.bestIdx;
		}
		fftw_free(queryMeanStdNext);
	}

	fftw_free(copyTa);
	fftw_free(queryMeanStd);
	fftw_free(TaMeanStdDev);	
	fftw_free(QT_b);
	fftw_free(QTChanging);
	fftw_free(minSizeSumSumSq);
	
	return bestMP;
}

void updateTopkDiscordsOld(heapEntry * topKDiscords, mMxHeap * minmaxEntry, int kDiscords, _Bool bFirstTime )
{

    int i=0;
    for(i=0;i<kDiscords;i++)
    {
         heapEntry * lMax = getHead(minmaxEntry);
        
         if(bFirstTime )
         {
             topKDiscords[(kDiscords-i)-1].distance = lMax->distance;
             topKDiscords[(kDiscords-i)-1].index1 = lMax->index1;
             topKDiscords[(kDiscords-i)-1].index2 = lMax->index2;
         }
         else if( topKDiscords[(kDiscords-i)-1].distance < lMax->distance )
         {
             topKDiscords[(kDiscords-i)-1].distance = lMax->distance;
             topKDiscords[(kDiscords-i)-1].index1 = lMax->index1;
             topKDiscords[(kDiscords-i)-1].index2 = lMax->index2;
         }
      
         removeHead(minmaxEntry);

    }
}


void updateTopkDiscords(heapEntry * topKDiscords, mMxHeap * minmaxEntry, int kDiscords, _Bool bFirstTime, int sizeSubSeq)
{

    heapEntry arrayHeaPentry[kDiscords];
    int i,j;
    for(i=0;i<kDiscords;i++)
    {
        heapEntry * lMax = getHead(minmaxEntry);
        arrayHeaPentry[(kDiscords-i)-1].distance = lMax->distance;
        arrayHeaPentry[(kDiscords-i)-1].index1 = lMax->index1;
        arrayHeaPentry[(kDiscords-i)-1].index2 = lMax->index2;
        removeHead(minmaxEntry);
    }

    for(i=0;i<kDiscords;i++)
    {
        _Bool bTrivial = false;
        for(j=0;j<kDiscords;j++)
        {
            if(j!=i && topKDiscords[j].distance >=0)
            {
                int diff1 = abs(topKDiscords[j].index1 - arrayHeaPentry[i].index1);
                //int diff2 = abs(topKDiscords[j].index2 - arrayHeaPentry[i].index2);

                bTrivial = bTrivial || !(diff1>=(sizeSubSeq/2));// || diff2>=(sizeSubSeq/2)); 
                if (bTrivial){
                    break;
                }
            }
        }
        
        if(!bTrivial)
        {
            if(bFirstTime || topKDiscords[i].distance<0)
            {
                topKDiscords[i].distance = arrayHeaPentry[i].distance;
                topKDiscords[i].index1 = arrayHeaPentry[i].index1;
                topKDiscords[i].index2 = arrayHeaPentry[i].index2;
            }
            else if( topKDiscords[i].distance < arrayHeaPentry[i].distance )
            {
                topKDiscords[i].distance = arrayHeaPentry[i].distance;
                topKDiscords[i].index1 = arrayHeaPentry[i].index1;
                topKDiscords[i].index2 = arrayHeaPentry[i].index2;
            }
        }
        
    } 
    
}




void saveKDiscords(heapEntry * topKDiscords, int kDiscords, int size, FILE * fLogDiscords)
{
    int i=0;
    fprintf(fLogDiscords, "size discord: %d\n",size);
    for(i=0;i<kDiscords;i++)
    {
        fprintf(fLogDiscords, "top-%d discords: distance: %lf, index1: %d, index2: %d\n",(i+1), sqrt(topKDiscords[i].distance),topKDiscords[i].index1 ,topKDiscords[i].index2 );
    }
}

// works only for the self join 
//TODO makes for JOIN
void VALMOD(VALMAP * matrixProfileVL , double * Ta, double * Tb, int sizeQuery, int sizeQueryMAX, int sizeTa, int sizeTb, _Bool selfJoin, const char * fileNameMP)
{
       
	FILE * fLog = fopen("Logrun.log","a");
        FILE * fLogDiscords = fopen("Discord.log","a");
        fprintf(fLogDiscords, "size discord for the output: %s\n\n",fileNameValmap);
	int matrixProfileSize = sizeTa-sizeQuery+1;
	int lastMatrixProfileSize = matrixProfileSize;
	mMxHeap ** listMotifProfile = malloc(sizeof(mMxHeap*)*matrixProfileSize);
	int i,j;
	for(i=0;i<matrixProfileSize;i++)
	{
		listMotifProfile[i] = newHeap(threshold,false);
	}
	int rangeSizes=((sizeQueryMAX-sizeQuery)+1);
	double * copyTb = fftw_malloc(sizeof(double)*2*sizeTa);
	int offset;
	int updatePerSize=0;
	int matrixProfilesElementsPerSize=0;
	int motifsPerSize=0;
 	int maxNumberMotifsFound = -1;
	int minNumberMotifsFound = -1;
	double * listLBMinNonValid = malloc(sizeof(double)*(matrixProfileSize));
 
        double * listTrueDistMinNonValid = malloc(sizeof(double)*(matrixProfileSize));
	int * indexLBMinNonValid = malloc(sizeof(int)*(matrixProfileSize));
        int * indexNonValidMax  = malloc(sizeof(int)*(matrixProfileSize));
        int avgDPUPMotifs = 0;
        int avgDPUPDiscords = 0;
        
	for(offset=0;offset<rangeSizes || motifsPerSize==0;offset++)
	{	
            int actualSizeQueryNext = (sizeQuery+offset);
            if (offset==0 || motifsPerSize==0)
            {
                    //####STOMP CALL OF VALMOD#####
                    motifsPerSize = (matrixProfileSize-offset);
                    // if offset > - means that I have to re-compute STOMP since there is no valid motif
                    if(offset > 0)
                    {
                        for(i=0;i<lastMatrixProfileSize;i++)
                        {
                                freeHeap(listMotifProfile[i]);
                        }
                        fftw_free(listMotifProfile);
                        offset--;
                        motifsPerSize = (matrixProfileSize-offset);
                        lastMatrixProfileSize = motifsPerSize;
                        listMotifProfile = malloc(sizeof(mMxHeap *)*lastMatrixProfileSize);
                        for(i=0;i<motifsPerSize;i++)
                        {
                                listMotifProfile[i] = newHeap(threshold,false);
                        }
                        actualSizeQueryNext = (sizeQuery+offset);
                    }

                    memcpy(copyTb, Tb, actualSizeQueryNext*sizeof(double)); // get the first query of B
                    //list * listMotifProfile = newList(); // this is the motifs profile list
                    double * QT = slidingDotProduct(Ta, copyTb, sizeTa, actualSizeQueryNext);
                    // compute the sum and squared sum of everything 
                    heapEntry bestMP = STOMP_MotifsProfile_VALMOD(matrixProfileVL, QT, offset, Ta, Tb, actualSizeQueryNext, sizeTa, sizeTb, selfJoin, fileNameMP, listMotifProfile);

                    double ln1 = 1.0/(double)actualSizeQueryNext;
                    double ln1SQ = sqrt(ln1);		
                    matrixProfileVL->distancesEvolutionMotif[offset] =  bestMP.distance * ln1SQ;


                    printf("Elements # in MP with length %d: %d STOMP (Best dist: %lf Q:%d , D:%d)\n",actualSizeQueryNext,motifsPerSize,sqrt(bestMP.distance),bestMP.index1,bestMP.index2);
                    fprintf(fLog,"Elements # in MP with length %d: %d STOMP (Best dist: %lf Q:%d , D:%d)\n",actualSizeQueryNext,motifsPerSize,sqrt(bestMP.distance),bestMP.index1,bestMP.index2);
                    fftw_free(QT);
		}
		else
		{
                     
                    heapEntry * topKDiscords= malloc(sizeof(heapEntry)*kDiscords); 
                    int kd;
                    for(kd=0;kd<kDiscords;kd++)
                    {
                        topKDiscords[kd].distance = -1 ;
                    }
                    
                    _Bool bInitKDiscords = true;
                    _Bool bInitKDiscordsNotValid = true;
                    //####USING THE LB FOR PRUNING COMPUTATION#####
                    updatePerSize=0; // it counts how many Distance profiles are not valid, (meaning : how many DP should be updated)
                    matrixProfilesElementsPerSize=0;
                    motifsPerSize=0;

                    heapEntry ** validEntries = malloc(sizeof(heapEntry *)*(matrixProfileSize-offset));
                    double  nonValidSmaller = -1;


                    double maxEntryNonValid = -1;
                    double minEntryNonValid = -1;
                    double minnABStrueDist = -1;
                    
                    int nonValidKMatch = 0;
                    int nonValidLB=0;
                    int indNVS=0;

                    double maxGeneral=-1;

                    int offsetMax1=0;
                    int offsetMax2=0;
                    mMxHeap ** nonValidDiscordMinmaxEntry=  malloc(sizeof(mMxHeap *)*(matrixProfileSize-offset)); // non valid dp first k matches
                    for(j=0;j<(matrixProfileSize-offset);j++)
                    {
                        // distance profile with offset j
                        heapEntry * lMax = getHead(listMotifProfile[j]);
                        double * mstq = computeMeanStd(&Tb[j], actualSizeQueryNext, actualSizeQueryNext);
                        double lowerbound = (lMax->lbDistance*(lMax->stdDevQuery*lMax->stdDevQuery))/(mstq[1]*mstq[1]);
                        heapEntry * minEntry = NULL;
                        double minEntryLB = -1;
                        heapEntry * newListSorted = listMotifProfile[j]->heap_;
                        mMxHeap * minmaxEntry=  newHeap(kDiscords,false); // k =10 for discords
                       
                        for(i=1;i<listMotifProfile[j]->size;i++)
                        {
                            // please check the TRIVIAL MATCHES
                            if(newListSorted[i].index1<(newListSorted[i].index2-(actualSizeQueryNext/2)) || 
                               newListSorted[i].index1>(newListSorted[i].index2+(actualSizeQueryNext/2)) || !selfJoin) // avoid the trivial matches
                            {
                                double dist = computeTrueDistanceEntryHeapEntry(Ta, Tb, sizeTa, sizeTb, &newListSorted[i], actualSizeQueryNext);
                                
                                if(dist==-1)
                                        continue; // not valid offset

                                 // to find the kth discords (store the k smallest nearest neighbor)
                                pushTopClassic(dist,j,newListSorted[i].index2,minmaxEntry);
                                
                                if(minEntry==NULL )
                                {
                                        minEntry = &newListSorted[i]; 
                                        minEntryLB = (newListSorted[i].lbDistance*(newListSorted[i].stdDevQuery*newListSorted[i].stdDevQuery))/(mstq[1]*mstq[1]);
                                        //newListSorted[i].lbDistance = minEntryLB;
                                }
                                else
                                {
                                        if(minEntry->distance>dist)
                                        {
                                                minEntry = &newListSorted[i]; 
                                        }
                                        double lb = (newListSorted[i].lbDistance*(newListSorted[i].stdDevQuery*newListSorted[i].stdDevQuery))/(mstq[1]*mstq[1]);
                                        if(minEntryLB < lb)
                                        {
                                                minEntryLB = lb; 
                                        }
                                }		
                            }
                        }

                        // pruning is effective for the best match
                        if(minEntry != NULL && minEntry->distance<lowerbound)
                        {
                            // check if the kth is smaller than the LB

                            if(	minnABStrueDist ==-1)
                            {
                                    minnABStrueDist = minEntry->distance;
                            }
                            else if(minEntry->distance < minnABStrueDist)
                            {
                                    minnABStrueDist = minEntry->distance;
                            }

                            validEntries[matrixProfilesElementsPerSize] = minEntry;
                            // new element of the matrix profile 
                            matrixProfilesElementsPerSize++;

                        }
                        // the distance profile has elements but nothing to say for the best match
                        else if(minEntryLB >=0)
                        {
                           
                            listTrueDistMinNonValid[nonValidLB] = minEntry->distance; 
                            listLBMinNonValid[nonValidLB] = minEntryLB;
                         
                            indexLBMinNonValid[nonValidLB] = j;
                            // store the non valid entry with the smallest LB. (I am sure that this guy correctly lowerbounds the MP  ) 
                            updatePerSize++;
                            if(nonValidSmaller == -1)
                            {
                                    indNVS=nonValidLB;
                                    nonValidSmaller = minEntryLB;
                            }
                            else 
                            {
                                    if(nonValidSmaller > minEntryLB)
                                    {
                                            indNVS=nonValidLB;
                                            nonValidSmaller = minEntryLB;	
                                    }
                            }
                            nonValidLB++;
                        }
                        else
                        {
                                // trival match
                                listLBMinNonValid[nonValidLB] = -1;
                                indexLBMinNonValid[nonValidLB] = j;
                        }
                        

                        
                        // pruning is effective for the kth match contained in minmax 
                        heapEntry * minmax = getHead(minmaxEntry);
                        if(minmax != NULL && minmax->distance<lowerbound)
                        {
 
                            updateTopkDiscords(topKDiscords, minmaxEntry, kDiscords,  bInitKDiscords ,actualSizeQueryNext);
                            if(bInitKDiscords)
                            {
                                bInitKDiscords = false;
                            }
                            freeHeap(minmaxEntry);
                        }
                        else 
                        {
                          
                            nonValidDiscordMinmaxEntry[nonValidKMatch] = minmaxEntry;
                            indexNonValidMax[nonValidKMatch] = j;
                            nonValidKMatch++;
                        }
                        fftw_free(mstq);                     
                    }

                    //check to see how many valid motifs (smallest value of the Matrix Profile) I have among the valid entries.
                    int m;
                    heapEntry * bestMotif = NULL;
                    for(m=0;m<matrixProfilesElementsPerSize;m++)
                    {
                        if(validEntries[m]->distance < nonValidSmaller)
                        {
                                updateValmap(matrixProfileVL, sqrt(validEntries[m]->distance), actualSizeQueryNext,validEntries[m]->index1, validEntries[m]->index2);

                                if(bestMotif == NULL)
                                {
                                        bestMotif = validEntries[m];
                                }
                                else
                                {
                                        if(bestMotif->distance > validEntries[m]->distance )
                                        {
                                                bestMotif = validEntries[m];
                                        }
                                }
                                motifsPerSize++;
                        }
                    }
                    //  update motifs 
                   
                    int totUpdateReal=0; // how many Ditance profiles I nned to update for finding the best motifs
                    if(bestMotif == NULL ) // matrix profile may not be containing the minimum
                    {   
                        //Here I do not want to call STOMP just update the distance profiles which have the max LB smaller than the minimum motif pair
                        // or I could not have the discords
                        int z;
                        double minLB = -1;
                        for(z=0;z<nonValidLB;z++)
                        {
                            if ((listLBMinNonValid[z] <=  minnABStrueDist))
                            {
                                int indexToUpdate = indexLBMinNonValid[z];
                                freeHeap(listMotifProfile[indexToUpdate]); 
                                listMotifProfile[indexToUpdate] =  newHeap(threshold,false);
                                memcpy(copyTb, &Tb[indexToUpdate], actualSizeQueryNext*sizeof(double)); 
                                double * QT = slidingDotProduct(Ta, copyTb, sizeTa, actualSizeQueryNext);
                                double * queryMeanStd = computeMeanStd(&Tb[indexToUpdate], actualSizeQueryNext, actualSizeQueryNext);
                                double * minSizeSumSumSq = computeSumSumSquared(Ta, sizeTa, actualSizeQueryNext);
                                double * querySumSumSq = computeSumSumSquared(&Tb[indexToUpdate], actualSizeQueryNext, actualSizeQueryNext);
                                double * TaMeanStdDev = computeMeanStd(Ta, sizeTa, actualSizeQueryNext);
                                returnDistanceProfileWithMotifsProfileComputation_VALMOD(QT, queryMeanStd, TaMeanStdDev, indexToUpdate,sizeTa,actualSizeQueryNext, selfJoin, listMotifProfile, minSizeSumSumSq, querySumSumSq);
                                int lm;

                                mMxHeap * minmaxEntry=  newHeap(kDiscords,false); // k =10 for discords
                                for(lm=1;lm<listMotifProfile[indexToUpdate]->size;lm++)
                                {
                                    // to find the kth discords (store the k smallest nearest neighbor)
                                    pushTopClassic(listMotifProfile[indexToUpdate]->heap_[lm].distance,indexToUpdate,listMotifProfile[indexToUpdate]->heap_[lm].index2,minmaxEntry);
                                }

                                validEntries[matrixProfilesElementsPerSize] = malloc(sizeof(heapEntry));
                                validEntries[matrixProfilesElementsPerSize]->distance = orderedidP.bestDist;
                                validEntries[matrixProfilesElementsPerSize]->index1 = indexToUpdate;
                                validEntries[matrixProfilesElementsPerSize]->index2 = orderedidP.bestIdx;
                                matrixProfilesElementsPerSize++;

                               
                                updateTopkDiscords(topKDiscords, minmaxEntry, kDiscords,  bInitKDiscords, actualSizeQueryNext );
                                if(bInitKDiscords)
                                {
                                    bInitKDiscords = false;
                                }
                                freeHeap(minmaxEntry);
                            
                                totUpdateReal++;
                                
                              

                                fftw_free(TaMeanStdDev);
                                fftw_free(QT);
                                fftw_free(queryMeanStd);
                                fftw_free(minSizeSumSumSq);
                                fftw_free(querySumSumSq);					
                            }
                            else if (listLBMinNonValid[z] >  minnABStrueDist)
                            {
                                    if(minLB<0)
                                    {
                                            minLB=listLBMinNonValid[z];
                                    }
                                    else if (minLB>listLBMinNonValid[z])
                                    {
                                            minLB = listLBMinNonValid[z];
                                    }
                            }
                        }



                        // update the new min max LB
                        nonValidSmaller = minLB;
                        //re - check 
                        for(m=0;m<matrixProfilesElementsPerSize;m++)
                        {
                            if(validEntries[m]->distance < nonValidSmaller)
                            {
                                updateValmap(matrixProfileVL, sqrt(validEntries[m]->distance), actualSizeQueryNext,  validEntries[m]->index1, validEntries[m]->index2);

                                if(bestMotif == NULL)
                                {
                                        bestMotif = validEntries[m];
                                }
                                else
                                {
                                        if(bestMotif->distance > validEntries[m]->distance )
                                        {
                                                bestMotif = validEntries[m];
                                        }
                                }
                                motifsPerSize++;
                            }
                        } 
                    }

                    if(bestMotif!=NULL)
                    {

                            double ln1 = 1.0/(double)actualSizeQueryNext;
                            double ln1SQ = sqrt(ln1);		
                            matrixProfileVL->distancesEvolutionMotif[offset] = bestMotif->distance * ln1SQ;
                    }

                    if (maxNumberMotifsFound < motifsPerSize)
                    { 
                            maxNumberMotifsFound = motifsPerSize;
                    }
                    if ( (minNumberMotifsFound > motifsPerSize || minNumberMotifsFound <0) )
                    { 
                            minNumberMotifsFound = motifsPerSize;
                    }		


                    if(matrixProfilesElementsPerSize > 0)
                    {
                            printf(" matrix profile elements: %d  !!!!  ", matrixProfilesElementsPerSize);
                    }

                   
                    if(motifsPerSize > 0)
                    {
                        printf("Size %d is over with %d non valid elements of MP and with %d motifs pair ! \nTOP MOTIF PAIR: distance: %lf q:%d, d:%d\n",actualSizeQueryNext,updatePerSize,motifsPerSize,sqrt(bestMotif->distance),bestMotif->index1,bestMotif->index2);
                        fprintf(fLog,"Size %d is over with %d non valid elements of MP and with %d motifs pair ! \nTOP MOTIF PAIR: distance: %lf q:%d, d:%d\n",actualSizeQueryNext,updatePerSize,motifsPerSize,sqrt(bestMotif->distance),bestMotif->index1,bestMotif->index2);
                    }
                    
                    
                    //############ DISCORDS PART
                    //#######################################
                    
                    int totUpdateRealDiscord=0;
                    int d=0;
                    for(d=0;d<nonValidKMatch;d++)
                    {
                        //printf("non valid iterate: %d\tot:%d , totUpdated:%d \n",d,nonValidKMatch,totUpdateRealDiscord);
                        int contDisc= (kDiscords-1);
                        // check if all the k elements non valid are less than the valid max
                        while(contDisc>=0)
                        {
                            heapEntry * h = getHead(nonValidDiscordMinmaxEntry[d]);
                            
                            if(h->distance > topKDiscords[contDisc].distance)
                            {
                                totUpdateRealDiscord++;
                                int indexToUpdate = indexNonValidMax[d];
                                freeHeap(listMotifProfile[indexToUpdate]); 
                                listMotifProfile[indexToUpdate] =  newHeap(threshold,false);
                                memcpy(copyTb, &Tb[indexToUpdate], actualSizeQueryNext*sizeof(double)); 
                                double * QT = slidingDotProduct(Ta, copyTb, sizeTa, actualSizeQueryNext);
                                double * queryMeanStd = computeMeanStd(&Tb[indexToUpdate], actualSizeQueryNext, actualSizeQueryNext);
                                double * minSizeSumSumSq = computeSumSumSquared(Ta, sizeTa, actualSizeQueryNext);
                                double * querySumSumSq = computeSumSumSquared(&Tb[indexToUpdate], actualSizeQueryNext, actualSizeQueryNext);
                                double * TaMeanStdDev = computeMeanStd(Ta, sizeTa, actualSizeQueryNext);
                                returnDistanceProfileWithMotifsProfileComputation_VALMOD(QT, queryMeanStd, TaMeanStdDev, indexToUpdate,sizeTa,actualSizeQueryNext, selfJoin, listMotifProfile, minSizeSumSumSq, querySumSumSq);
                                int lm;

                                mMxHeap * minmaxEntry=  newHeap(kDiscords,false); // k =10 for discords
                                for(lm=1;lm<listMotifProfile[indexToUpdate]->size;lm++)
                                {
                                    // to find the kth discords (store the k smallest nearest neighbor)
                                    pushTopClassic(listMotifProfile[indexToUpdate]->heap_[lm].distance,indexToUpdate,listMotifProfile[indexToUpdate]->heap_[lm].index2,minmaxEntry);
                                }

                                  // store the potential discovery (among the valid DP)
                                updateTopkDiscords(topKDiscords, minmaxEntry, kDiscords,  bInitKDiscords, actualSizeQueryNext );
                                if(bInitKDiscords)
                                {
                                    bInitKDiscords = false;
                                }
                                
                                fftw_free(TaMeanStdDev);
                                fftw_free(QT);
                                fftw_free(queryMeanStd);
                                fftw_free(minSizeSumSumSq);
                                fftw_free(querySumSumSq);
                                
                                freeHeap(minmaxEntry);
                                break;
                            }
                            removeHead(nonValidDiscordMinmaxEntry[d]);
                            contDisc--;
                        }
                    }

                    for(d=0;d<nonValidKMatch;d++)
                    {
                        freeHeap(nonValidDiscordMinmaxEntry[d]);
                    }
                    
                    avgDPUPMotifs = avgDPUPMotifs + totUpdateReal;
                    avgDPUPDiscords = avgDPUPDiscords + totUpdateRealDiscord;
                    
                    fftw_free(nonValidDiscordMinmaxEntry);
                    saveKDiscords(topKDiscords,kDiscords,actualSizeQueryNext,fLogDiscords);
                    
                    fprintf(fLog,"TOP K DISCORD PAIR  : (distance: %lf q:%d, d:%d)\n",sqrt(topKDiscords[kDiscords-1].distance),topKDiscords[kDiscords-1].index1,topKDiscords[kDiscords-1].index2);
                    printf("TOP K DISCORD PAIR : (distance: %lf q:%d, d:%d)\n",sqrt(topKDiscords[kDiscords-1].distance),topKDiscords[kDiscords-1].index1,topKDiscords[kDiscords-1].index2);
                    printf("DP updates discords : %d \n",totUpdateRealDiscord);
                    printf("DP updates motifs : %d\n",totUpdateReal);
                    fprintf(fLog,"DP updates discords : %d \n",totUpdateRealDiscord);
                    fprintf(fLog,"DP updates motifs : %d\n",totUpdateReal);
                     //############ END DISCORDS PART
                    //#######################################

                    int f;
                    // free the extra created entries
                    for(f=0;f<totUpdateReal;f++)
                    {
                            free(validEntries[matrixProfilesElementsPerSize-f-1]);
                    }
                    // free validEntries
                    fftw_free(validEntries);
                    fftw_free(topKDiscords);
          
                    
            }	
	}
	free(listTrueDistMinNonValid);
	free(listLBMinNonValid);
	free(indexLBMinNonValid);
        free(indexNonValidMax);
	for(i=0;i<lastMatrixProfileSize;i++)
	{
		 freeHeap(listMotifProfile[i]);
	}
	
        avgDPUPDiscords = avgDPUPDiscords/(rangeSizes-1);
        avgDPUPMotifs = avgDPUPMotifs/(rangeSizes-1);
        
        
        printf("Statistic: total points series:%d, l_min: %d, l_max: %d,\n \
                distance profile points memory:%d , number discords to find: %d \n ",sizeTa,sizeQuery,sizeQueryMAX,threshold, kDiscords);
	
         
        fprintf(fLog,"Statistic: total points series:%d, l_min: %d, l_max: %d,\n \
                distance profile points memory:%d , number discords to find: %d \n ",sizeTa,sizeQuery,sizeQueryMAX,threshold, kDiscords);
        
        
	printf("MAX MOTIFS FOUND LINEARLY: %d, MIN MOTIFS FOUND LINEARLY: %d\n ",maxNumberMotifsFound,minNumberMotifsFound);
	fprintf(fLog,"MAX MOTIFS FOUND LINEARLY: %d, MIN MOTIFS FOUND LINEARLY: %d\n ",maxNumberMotifsFound,minNumberMotifsFound);
            
        fprintf(fLog,"average DP update for motifs %d\n ",avgDPUPDiscords);
        fprintf(fLog,"average DP update for discords %d\n ",avgDPUPMotifs);
        
        printf("average DP update for discords %d\n ",avgDPUPDiscords);
        printf("average DP update for motifs %d\n ",avgDPUPMotifs);
        
        
	fftw_free(listMotifProfile);
	fftw_free(copyTb);
	fclose(fLog);
        fclose(fLogDiscords);
      
        
}


double * computeMaxForLB(double * in, int sizeIn, int minSize, int maxSize)
{
	double * maxSizes = calloc((maxSize-minSize),sizeof(double));
	int i,j;
	int cont=0;
	for(j=minSize;j<maxSize;j++)
	{
		double * inMeanStd = computeMeanStd(in, sizeIn, j);
		for(i=0;i<(sizeIn-j);i++)
		{
			double y = (in[i+j] - inMeanStd[i] ) / inMeanStd[(i+((sizeIn-j)+1))];
			if(y>maxSizes[cont])
			{
				maxSizes[cont] = y ;
			}
		}
		fftw_free(inMeanStd);
		cont++;
	}
	
	return maxSizes;
}





// toy example to test the implementation of STOMP and STAMP
int selfJoinWith_STAMP_and_STOMP()
{
	int i;
	double *in; // the query 
	double *in2; // the data
	int n1 = 1000; // length of in
	int n2 = 1000; // length of in2
	int qLength = 50;
	
	in = fftw_malloc ( sizeof ( double ) * (n1*2) );
	in2 = fftw_malloc ( sizeof ( double ) * (n2*2) );

	in[0] = randPersonal(0,1);
	
	for ( i = 1; i < n1; i++ )
	{
		in[i] = in[i-1] + randPersonal(0,1);
	}
	
	/*
	in2[0] = randPersonal(0,1);
	for ( i = 1; i < n2; i++ )
	{
		in2[i] = in2[i-1] + randPersonal(0,1);
	}
	
	const char * fileName1 = "JoinAB.mp";
	
	// Join A != B 
	mpElements e3 = STAMP(in,in2, qLength, n1,n2,false,fileName1);
	
	fftw_free(e3.matrixProfile);
	fftw_free(e3.indxProfile);
	*/
	
	memcpy(in2, in, n1*(sizeof(double))); // make the 2 series the same for the selfjoin
	const char * fileName2 = "SelfJoin50.mp";
	//Self join A size 50 
	double t1,t2;
	t1=clock();
	mpElements e1 = STAMP(in,in2, qLength, n1,n1,true,fileName2);
	t2=clock();
	
	//printf("Total Execution Time of STAMP: %lf sec \n ",(t2-t1)/CLOCKS_PER_SEC);
	
	rangeMpElements rMP;
	rMP.matrixProfiles = malloc(sizeof(double*));
	rMP.indexProfiles = malloc(sizeof(int*));
	rMP.matrixProfiles[0] = malloc(sizeof(double)*((n1-qLength)+1));
	rMP.indexProfiles[0] = malloc(sizeof(int)*((n1-qLength)+1));
        rMP.standardProfiles[0] = malloc(sizeof(double)*((n1-qLength)+1));
	const char * fileName2b = "SelfJoin50STOMP.mp";
	t1=clock();
		STOMP(rMP,0,in,in2, qLength, n1,n1,true,fileName2b);
	t2=clock();
	
	//printf("Total Execution Time of STOMP: %lf sec \n ",(t2-t1)/CLOCKS_PER_SEC);
	
	fftw_free(e1.matrixProfile);
	fftw_free(e1.indxProfile);
	fftw_free(rMP.matrixProfiles[0]);
	fftw_free(rMP.indexProfiles[0]);
        fftw_free(rMP.standardProfiles[0]);
	fftw_free(rMP.matrixProfiles);
	fftw_free(rMP.indexProfiles);
        fftw_free(rMP.standardProfiles);
	fftw_free ( in );
	fftw_free ( in2 );

	return 1;
}

int basicExperiment()
{
	//###########################################################################################
	//Compute the dot product with fft
	
	int i;
	double *in; // the query 
	double *in2; // the data
	double *originalSeries; // the query 
	double *originalquery; // the data
	int n = 5000000; // length of the long series
	int qLength = 5000; // length of the query
	
	in = fftw_malloc ( sizeof ( double ) * (n*2) );
	in2 = fftw_malloc ( sizeof ( double ) * (n*2) );

	originalSeries = fftw_malloc ( sizeof ( double ) * n );
	originalquery = fftw_malloc ( sizeof ( double ) * qLength );
	
	in[0] = randPersonal(0,1);
	in2[0] = randPersonal(0,1);
	originalSeries[0] = in[0];	
	originalquery[0] = in2[0];
	
	for ( i = 1; i < n; i++ )
	{
		in[i] = in[i-1] + randPersonal(0,1);
		originalSeries[i] = in[i];
	}
	for ( i = 1; i < qLength; i++ )
	{
		in2[i] = in2[i-1] + randPersonal(0,1);
		originalquery[i] = in2[i];
	}

	double t1,t2;
	t1= clock();
	double * dProducts = slidingDotProduct(in,in2, n, qLength);
	t2= clock();
	//printf("Total Execution Time in frequency domain: %lf sec \n ",(t2-t1)/CLOCKS_PER_SEC);
	
	//###########################################################################################
	//Compute naively the dot product
	//printf ( "dot product computed as usual:\n");
	double t3,t4;
	t3= clock();
	int j;
	for ( i = 0; i < ((n-qLength)+1); i++ )
	{
		double dp=0;
		for ( j = 0; j < qLength; j++ )
		{
			dp = dp + (originalquery[j] * originalSeries[i+j]);
		}
		//printf ( "%12f \n", dp);
	}
	t4= clock();
	//printf("Total Execution Time in time domain: %lf sec \n ",(t4-t3)/CLOCKS_PER_SEC);

	//Release the memory associated with the plans.
	
	fftw_free ( in );
	fftw_free ( in2 );
	fftw_free ( dProducts );
	fftw_free ( originalSeries );
	fftw_free ( originalquery );
}

// join A,B


mpElements * _join_STOMP(int sizeSeries, int sizeSeriesB, int length, char * fileTSA, char * fileTSB, _Bool bPrintMP)
{
	int i;
	double *in; // the quey 
	double *in2; // the data
	int n1 = 10; // length of in
	int n2 = 10; // length of in2
	int qLength = 2;
	
	if(sizeSeries>0 && length>0)
	{
		n1 = sizeSeries; // length of in
		n2 = sizeSeriesB; // length of in2
		qLength = length;
	}
	
	in = fftw_malloc ( sizeof ( double ) * (n1*2) );
	in2 = fftw_malloc ( sizeof ( double ) * (n2*2) );

	if(fileTSA==NULL)
	{
            in[0] = randPersonal(0,1) + 5;
            FILE * fLastSeriesUsed = fopen("../DATA/lastSeriesUsed.ts","w");
            fprintf(fLastSeriesUsed,"%lf\n",in[0]);
            for ( i = 1; i < n1; i++ )
            {
                    in[i] = in[i-1] + (randPersonal(0,1)/20);
                    //in[i] =  randPersonal(0,1);
                    fprintf(fLastSeriesUsed,"%lf\n",in[i]);
            }
            fclose(fLastSeriesUsed);

            FILE * fLastSeriesUsedBIN = fopen("../DATA/lastSeriesUsed.bin","wb");
            fwrite(&in[0],sizeof(double),1,fLastSeriesUsedBIN);
            for ( i = 1; i < n1; i++ )
            {
                    //in[i] =  randPersonal(0,1);
                    fwrite(&in[i],sizeof(double),1,fLastSeriesUsedBIN);
            }
            fclose(fLastSeriesUsedBIN);
	}
	else
	{
            if (strstr(fileTSA, ".bin") != NULL) 
            {
                    FILE * ts = fopen(fileTSA,"rb");
                    float point;
                    int cont=0;
                    while (cont<n1) 
                    {
                            fread(&in[cont], sizeof(double), 1, ts);
                            cont++;
                    }
                    //printf("first point of the series: %lf\n",in[0]);
                    fclose(ts);	
            }
            else
            {
                    FILE * ts = fopen(fileTSA,"r");
                    char line[256];
                    int cont=0;
                    while (fgets(line, sizeof(line), ts) && cont<n1) 
                    {
                            /* note that fgets don't strip the terminating \n, checking its
                               presence would allow to handle lines longer that sizeof(line) */
                            in[cont] = atof(line);
                            cont++;
                    }
                    //printf("first point of the series: %lf\n",in[0]);
                    fclose(ts);	
            }
	}

        if(fileTSB==NULL)
	{
            in2[0] = randPersonal(0,1) + 5;
            FILE * fLastSeriesUsed = fopen("../DATA/lastSeriesUsed2.ts","w");
            fprintf(fLastSeriesUsed,"%lf\n",in2[0]);
            for ( i = 1; i < n2; i++ )
            {
                in2[i] = in2[i-1] + (randPersonal(0,1)/20);
                fprintf(fLastSeriesUsed,"%lf\n",in2[i]);
            }
            fclose(fLastSeriesUsed);

            FILE * fLastSeriesUsedBIN = fopen("../DATA/lastSeriesUsed2.bin","wb");
            fwrite(&in2[0],sizeof(double),1,fLastSeriesUsedBIN);
            for ( i = 1; i < n2; i++ )
            {
                fwrite(&in2[i],sizeof(double),1,fLastSeriesUsedBIN);
            }
            fclose(fLastSeriesUsedBIN);
	}
	else
	{
            if (strstr(fileTSB, ".bin") != NULL) 
            {
                FILE * ts = fopen(fileTSB,"rb");
                float point;
                int cont=0;
                while (cont<n2) 
                {
                    fread(&in2[cont], sizeof(double), 1, ts);
                    cont++;
                }
                //printf("first point of the series: %lf\n",in2[0]);
                fclose(ts);	
            }
            else
            {
                FILE * ts = fopen(fileTSB,"r");
                char line[256];
                int cont=0;
                while (fgets(line, sizeof(line), ts) && cont<n2) 
                {
                    /* note that fgets don't strip the terminating \n, checking its
                       presence would allow to handle lines longer that sizeof(line) */
                    in2[cont] = atof(line);
                    cont++;
                }
                //printf("first point of the series: %lf\n",in2[0]);
                fclose(ts);	
            }
	}
        
	
	//printf("size seriesA : %d, size seriesB : %d, query:%d \n",n1,n2,qLength);

	double t1,t2;

        rangeMpElements rMP;
        
        rMP.matrixProfiles = malloc(sizeof(double*));
        rMP.indexProfiles = malloc(sizeof(int*));
        rMP.standardProfiles = malloc(sizeof(double*));
        rMP.matrixProfiles[0] = malloc(sizeof(double)*(n2-(qLength)+1));
        rMP.indexProfiles[0] = malloc(sizeof(int)*(n2-(qLength))+1);
        rMP.standardProfiles[0] = malloc(sizeof(double)*(n2-(qLength))+1);

        double t3,t4;

        FILE * fLog = fopen("Logrun.log","a");
        
        char * nameFileMP = NULL;
        if(bPrintMP)
        {
            nameFileMP = malloc(500);
            snprintf(nameFileMP, 500,"JOIN_%d_%d_%d_last.ts",n1,n2,qLength);
        }

        t3=clock();

        heapEntry bestMP = STOMP(rMP,0,in,in2,qLength,n1,n2,false,nameFileMP);
        
        //printf("STOMP with subsequence length  %d over with best distance: %lf, Q:%d, D:%d \n",qLength,sqrt(bestMP.distance),bestMP.index1,bestMP.index2);
        //printf("STOMP with subsequence length  %d over with discord distance: %lf, Q:%d, D:%d \n",qLength,sqrt(discord.distance),discord.index1,discord.index2);
        if(fileTSA==NULL)
        {
            fprintf(fLog,"STOMP JOIN on 2 RANDOM WALKS, length SERIES_A: %d, length SERIES_B: %d, with subsequence length: %d OVER with best distance: %lf, Q:%d, D:%d \n",sizeSeries,sizeSeriesB,qLength,bestMP.distance,bestMP.index1,bestMP.index2);
        }
        else
        {
            fprintf(fLog,"STOMP JOIN: %s on %s, length SERIES_A: %d, length SERIES_B: %d, with subsequence length: %d OVER with best distance: %lf, Q:%d, D:%d \n",fileTSA,fileTSB,sizeSeries,sizeSeriesB,qLength,bestMP.distance,bestMP.index1,bestMP.index2);
        }
        t4=clock();
        
        if(!(nameFileMP==NULL))
        {
            fftw_free(nameFileMP);
        }
        
        fclose(fLog);

        //printf("Total Execution Time of STOMP: %lf sec \n ",(t4-t3)/CLOCKS_PER_SEC);
        fLog = fopen("Logrun.log","a");
        if(fileTSA==NULL)
        {
            fprintf(fLog,"Total Execution Time of STOMP JOIN  on 2 Random Walks DATA: %lf sec \n ",(t4-t3)/CLOCKS_PER_SEC);
        }
        else
        {
            fprintf(fLog,"Total Execution Time of STOMP JOIN %s on %s: %lf sec \n ",fileTSA,fileTSB,(t4-t3)/CLOCKS_PER_SEC);
        }

        fclose(fLog);


        mpElements * mPiP = malloc(sizeof(mpElements));
        mPiP->matrixProfile= rMP.matrixProfiles[0];
        mPiP->indxProfile=  rMP.indexProfiles[0];
	mPiP->standardProfile = rMP.standardProfiles[0];
        fftw_free(rMP.matrixProfiles);
        fftw_free(rMP.indexProfiles);	
        fftw_free(rMP.standardProfiles);	
	fftw_free(in);
	fftw_free(in2);
	return mPiP;
}



mpElements * _selfjoin_STOMP(int sizeSeries, int length, char * fileTSA, _Bool bPrintMP)
{
	int i;
	double *in; // the quey 
	double *in2; // the data
	int n1 = 10; // length of in
	int n2 = 10; // length of in2
	int qLength = 2;
	
	if(sizeSeries>0 && length>0)
	{
		n1 = sizeSeries; // length of in
		n2 = sizeSeries; // length of in2
		qLength = length;
	}
	
	in = fftw_malloc ( sizeof ( double ) * (n1*2) );
	in2 = fftw_malloc ( sizeof ( double ) * (n2*2) );

	if(fileTSA==NULL)
	{
            in[0] = randPersonal(0,1) + 5;
            FILE * fLastSeriesUsed = fopen("../DATA/lastSeriesUsed.ts","w");
            fprintf(fLastSeriesUsed,"%lf\n",in[0]);
            for ( i = 1; i < n1; i++ )
            {
                    in[i] = in[i-1] + (randPersonal(0,1)/20);
                    //in[i] =  randPersonal(0,1);
                    fprintf(fLastSeriesUsed,"%lf\n",in[i]);
            }
            fclose(fLastSeriesUsed);

            FILE * fLastSeriesUsedBIN = fopen("../DATA/lastSeriesUsed.bin","wb");
            fwrite(&in[0],sizeof(double),1,fLastSeriesUsedBIN);
            for ( i = 1; i < n1; i++ )
            {
                    //in[i] =  randPersonal(0,1);
                    fwrite(&in[i],sizeof(double),1,fLastSeriesUsedBIN);
            }
            fclose(fLastSeriesUsedBIN);
	}
	else
	{
            if (strstr(fileTSA, ".bin") != NULL) 
            {
                    FILE * ts = fopen(fileTSA,"rb");
                    float point;
                    int cont=0;
                    while (cont<n1) 
                    {
                            fread(&in[cont], sizeof(double), 1, ts);
                            cont++;
                    }
                    //printf("first point of the series: %lf\n",in[0]);
                    fclose(ts);	
            }
            else
            {
                    FILE * ts = fopen(fileTSA,"r");
                    char line[256];
                    int cont=0;
                    while (fgets(line, sizeof(line), ts) && cont<n1) 
                    {
                            /* note that fgets don't strip the terminating \n, checking its
                               presence would allow to handle lines longer that sizeof(line) */
                            in[cont] = atof(line);
                            cont++;
                    }
                    //printf("first point of the series: %lf\n",in[0]);
                    fclose(ts);	
            }
	}

        memcpy(in2,in,sizeof ( double ) * (n2*2));
        
	//printf("size seriesA : %d, query:%d \n",n1,qLength);

	double t1,t2;

        rangeMpElements rMP;
        rMP.matrixProfiles = malloc(sizeof(double*));
        rMP.indexProfiles = malloc(sizeof(int*));
        rMP.standardProfiles = malloc(sizeof(double*));
        rMP.matrixProfiles[0] = malloc(sizeof(double)*(n2-(qLength)+1));
        rMP.indexProfiles[0] = malloc(sizeof(int)*(n2-(qLength))+1);
        rMP.standardProfiles[0] = malloc(sizeof(double)*(n2-(qLength))+1);

        double t3,t4;


        FILE * fLog = fopen("Logrun.log","a");
        char * nameFileMP = NULL;
        if(bPrintMP)
        {
            nameFileMP = malloc(500);
            snprintf(nameFileMP, 500,"SELFJOIN_%d_%d_last.ts",n1,qLength);
        }
        
        t3=clock();
        heapEntry bestMP = STOMP(rMP,0,in,in2,qLength,n1,n2,true,nameFileMP);

        //printf("STOMP with subsequence length  %d over with best distance: %lf, Q:%d, D:%d \n",qLength,sqrt(bestMP.distance),bestMP.index1,bestMP.index2);
        //printf("STOMP with subsequence length  %d over with discord distance: %lf, Q:%d, D:%d \n",qLength,sqrt(discord.distance),discord.index1,discord.index2);
        if(fileTSA==NULL)
        {
            fprintf(fLog,"STOMP SELF JOIN on RANDOM WALK, length SERIES_A: %d, with subsequence length: %d OVER with best distance: %lf, Q:%d, D:%d \n",sizeSeries,qLength,bestMP.distance,bestMP.index1,bestMP.index2);
        }
        else
        {
            fprintf(fLog,"STOMP SELF JOIN: %s length SERIES_A: %d, with subsequence length: %d OVER with best distance: %lf, Q:%d, D:%d \n",fileTSA,sizeSeries,qLength,bestMP.distance,bestMP.index1,bestMP.index2);
        }
        t4=clock();
        
        if(!(nameFileMP==NULL))
        {
            fftw_free(nameFileMP);
        }
        fclose(fLog);

        //printf("Total Execution Time of STOMP: %lf sec \n ",(t4-t3)/CLOCKS_PER_SEC);
        fLog = fopen("Logrun.log","a");
        if(fileTSA==NULL)
        {
            fprintf(fLog,"Total Execution Time of STOMP SELF JOIN Random Walk DATA: %lf sec \n ",(t4-t3)/CLOCKS_PER_SEC);
        }
        else
        {
            fprintf(fLog,"Total Execution Time of STOMP SELF JOIN %s : %lf sec \n ",fileTSA,(t4-t3)/CLOCKS_PER_SEC);
        }

        fclose(fLog);


        mpElements * mPiP = malloc(sizeof(mpElements));
        mPiP->matrixProfile= rMP.matrixProfiles[0];
        mPiP->indxProfile=  rMP.indexProfiles[0];
	mPiP->standardProfile = rMP.standardProfiles[0];
        fftw_free(rMP.matrixProfiles);
        fftw_free(rMP.indexProfiles);	
        fftw_free(rMP.standardProfiles);	
	fftw_free(in);
	fftw_free(in2);
	return mPiP;
}


int selfJoinWith_STOMPmultilength_VALMOD(int sizeSeries,int min,int max, char * fileTS)
{
	int i;
	double *in; // the quey 
	double *in2; // the data
	int n1 = 10; // length of in
	int n2 = 10; // length of in2
	int qLengthMin = 2;
	int qLengthMax = 4;
	
	if(sizeSeries>0 && min>0 && max>0)
	{
		n1 = sizeSeries; // length of in
		n2 = sizeSeries; // length of in2
		qLengthMin = min;
		qLengthMax = max;
	}
	
	in = fftw_malloc ( sizeof ( double ) * (n1*2) );
	in2 = fftw_malloc ( sizeof ( double ) * (n2*2) );

	if(fileTS==NULL)
	{
		in[0] = randPersonal(0,1) + 5;
		FILE * fLastSeriesUsed = fopen("../DATA/lastSeriesUsed.ts","w");
		fprintf(fLastSeriesUsed,"%lf\n",in[0]);
		for ( i = 1; i < n1; i++ )
		{
			in[i] = in[i-1] + (randPersonal(0,1)/20);
			//in[i] =  randPersonal(0,1);
			fprintf(fLastSeriesUsed,"%lf\n",in[i]);
		}
		fclose(fLastSeriesUsed);

		FILE * fLastSeriesUsedBIN = fopen("../DATA/lastSeriesUsed.bin","wb");
		fwrite(&in[0],sizeof(double),1,fLastSeriesUsedBIN);
		for ( i = 1; i < n1; i++ )
		{
			//in[i] =  randPersonal(0,1);
			fwrite(&in[i],sizeof(double),1,fLastSeriesUsedBIN);
		}
		fclose(fLastSeriesUsedBIN);
	}
	else
	{
		if (strstr(fileTS, ".bin") != NULL) 
		{
			FILE * ts = fopen(fileTS,"rb");
			float point;
			int cont=0;
			while (cont<n1) 
			{
				fread(&in[cont], sizeof(double), 1, ts);
				cont++;
			}
			//printf("first point of the series: %lf\n",in[0]);
			fclose(ts);	
		}
		else
		{
			FILE * ts = fopen(fileTS,"r");
			char line[256];
			int cont=0;
			while (fgets(line, sizeof(line), ts) && cont<n1) 
			{
				/* note that fgets don't strip the terminating \n, checking its
				   presence would allow to handle lines longer that sizeof(line) */
				in[cont] = atof(line);
				cont++;
			}
			//printf("first point of the series: %lf\n",in[0]);
			fclose(ts);	
		}
	}

	memcpy(in2, in, n1*(sizeof(double))); // make the 2 series the same for the selfjoin
	printf("size series: %d, min query: %d, max query:%d \n",n1,qLengthMin,qLengthMax);

	double t1,t2;
	int o;
	int rangeSizes=((qLengthMax-qLengthMin)+1);

	int lengthVALMAP = n1-qLengthMin+1;
	VALMAP * vlMP =  malloc(sizeof(VALMAP));
	vlMP->matrixProfile =  malloc(sizeof(double)*lengthVALMAP);
	vlMP->indexProfile =  malloc(sizeof(double)*lengthVALMAP);
	vlMP->lengthProfile =  malloc(sizeof(double)*lengthVALMAP);
	vlMP->distancesEvolutionMotif =  malloc(sizeof(double)*rangeSizes);
	vlMP->matrixProfileNonLengthNormalized =  malloc(sizeof(double)*(n1-qLengthMin+1));
	vlMP->indexProfileNonLengthNormalized =  malloc(sizeof(int)*(n1-qLengthMin+1));
	vlMP->lengthProfileNonLengthNormalized =  malloc(sizeof(int)*(n1-qLengthMin+1));
	
	FILE * fLog = NULL; 
	
	if(!bOnlySTOMP)
	{
		/////////////////////////////////////////////////////////////
		t1=clock();
			VALMOD(vlMP, in, in2, qLengthMin, qLengthMax, n1, n1, true, NULL);
		t2=clock();
	
		//printf("Total Execution Time of VALMOD_DP_VALMOD_BOUND_EXACT with %d lengths : %lf sec \n ",(qLengthMax-qLengthMin)+1,(t2-t1)/CLOCKS_PER_SEC);

		fLog = fopen("Logrun.log","a");
		if(fileTS==NULL)
		{
			fprintf(fLog,"Total Execution Time of VALMOD_DP_VALMOD_BOUND_EXACT  with %d lengths on RANDOM walk length: %d-> %lf sec \n ",(qLengthMax-qLengthMin)+1,sizeSeries,(t2-t1)/CLOCKS_PER_SEC);
		}
		else
		{
			fprintf(fLog,"Total Execution Time of VALMOD_DP_VALMOD_BOUND_EXACT  with %d lengths on: %s , length: %d-> %lf sec \n ",(qLengthMax-qLengthMin)+1,fileTS,sizeSeries,(t2-t1)/CLOCKS_PER_SEC);
		}
		fclose(fLog);
	}
	
	
	if(bComparisonStomp || bOnlySTOMP)
	{
		rangeMpElements rMP;
		rMP.matrixProfiles = malloc(sizeof(double*)*((qLengthMax-qLengthMin)+1));
		rMP.indexProfiles = malloc(sizeof(int*)*((qLengthMax-qLengthMin)+1));
                rMP.standardProfiles = malloc(sizeof(double*)*((qLengthMax-qLengthMin)+1));

		for(o=0;o<rangeSizes;o++)
		{
			rMP.matrixProfiles[o] = malloc(sizeof(double)*((n1-(qLengthMin+o))+1));
			rMP.indexProfiles[o] = malloc(sizeof(int)*((n1-(qLengthMin+o))+1));
                        rMP.standardProfiles[o] = malloc(sizeof(double)*((n1-(qLengthMin+o))+1));
		}
		
		double t3,t4;
		int j;
		int cont=0;
		fLog = fopen("Logrun.log","a");
		t3=clock();

		for(j=qLengthMin;j<=qLengthMax;j++)
		{
			char * nameFileMP = malloc(500);
			snprintf(nameFileMP, 500,"%s_l_%d_MATRIX_PROFILE.ts",fileNameValmap,j);
			heapEntry bestMP = STOMP(rMP,cont,in,in2,j,n1,n1,true,nameFileMP);
			//printf("STOMP with subsequence length  %d over with best distance: %lf, Q:%d, D:%d \n",j,sqrt(bestMP.distance),bestMP.index1,bestMP.index2);
			//printf("STOMP with subsequence length  %d over with discord distance: %lf, Q:%d, D:%d \n",j,sqrt(discord.distance),discord.index1,discord.index2);
			if(fileTS==NULL)
			{
				fprintf(fLog,"STOMP on RANDOM WALK, length SERIES: %d, with subsequence length: %d OVER with best distance: %lf, Q:%d, D:%d \n",sizeSeries,j,bestMP.distance,bestMP.index1,bestMP.index2);
			}
			else
			{
				fprintf(fLog,"STOMP on: %s , length SERIES: %d, with subsequence length: %d OVER with best distance: %lf, Q:%d, D:%d \n",fileTS,sizeSeries,j,bestMP.distance,bestMP.index1,bestMP.index2);
			}
			fftw_free(nameFileMP);
			cont++;
		}
		t4=clock();
		fclose(fLog);

		//printf("Total Execution Time of STOMP with %d lengths: %lf sec \n ",(qLengthMax-qLengthMin)+1,(t4-t3)/CLOCKS_PER_SEC);
		fLog = fopen("Logrun.log","a");
		if(fileTS==NULL)
		{
			fprintf(fLog,"Total Execution Time of STOMP on Random Walk DATA with %d lengths: %lf sec \n ",(qLengthMax-qLengthMin)+1,(t4-t3)/CLOCKS_PER_SEC);
		}
		else
		{
			fprintf(fLog,"Total Execution Time of STOMP on %s with %d lengths: %lf sec \n ",fileTS,(qLengthMax-qLengthMin)+1,(t4-t3)/CLOCKS_PER_SEC);
		}
		
		fclose(fLog);

		for(o=0;o<((qLengthMax-qLengthMin)+1);o++) // cycle the MPs
		{
			fftw_free(rMP.matrixProfiles[o]);
			fftw_free(rMP.indexProfiles[o]);	
                        fftw_free(rMP.standardProfiles[o]);		
		}	

		fftw_free(rMP.matrixProfiles);
		fftw_free(rMP.indexProfiles);	
                fftw_free(rMP.standardProfiles);	

	}
	
	// store the variable length matrix profile, the index profile and the length profile

	char * nameFileVALMAP = malloc(500);
	char * nameFileIndexProfile = malloc(500);
	char * nameFilelengthProfile = malloc(500);
	char * nTopTenMatrixProfile = malloc(500);
        char * nTopTenMatrixProfileNN = malloc(500);
	char * evolutionMotifPair = malloc(500);
	

	snprintf(nameFileVALMAP, 500,"%s_a_%d_b_%d_T%d_VALMAP.ts",fileNameValmap,min,max,sizeSeries);
	snprintf(nameFileIndexProfile, 500, "%s_a_%d_b_%d_T%d_IndexProfile.ts",fileNameValmap,min,max,sizeSeries);
	snprintf(nameFilelengthProfile, 500, "%s_a_%d_b_%d_T%d_lengthProfile.ts",fileNameValmap,min,max,sizeSeries);
	snprintf(nTopTenMatrixProfile, 500, "%s_a_%d_b_%d_T%d_TopTenVALMAP.txt",fileNameValmap,min,max,sizeSeries);
	snprintf(nTopTenMatrixProfileNN, 500, "%s_a_%d_b_%d_T%d_TopTenVALMAP_NN.txt",fileNameValmap,min,max,sizeSeries);
        snprintf(evolutionMotifPair, 500, "%s_a_%d_b_%d_T%d_EvolutionMotif.ts",fileNameValmap,min,max,sizeSeries);
	

	
	FILE *  fVALMAP = fopen(nameFileVALMAP, "a");
	FILE *  fIndexProfile = fopen(nameFileIndexProfile, "a");
	FILE *  flengthProfile = fopen(nameFilelengthProfile, "a");
	FILE *  fTopTen = fopen(nTopTenMatrixProfile, "a");
        FILE *  fTopTenNN = fopen(nTopTenMatrixProfileNN, "a");
	FILE *  fEvolution = fopen(evolutionMotifPair, "a");
	
	int l;
	for(l=0;l<lengthVALMAP;l++)
	{
		if(l==0)
		{
			fprintf(fVALMAP,"%lf",vlMP->matrixProfile[l]);
			fprintf(fIndexProfile,"%d",vlMP->indexProfile[l]);
			fprintf(flengthProfile,"%d",vlMP->lengthProfile[l]);
		}
		else
		{
			fprintf(fVALMAP,"\n%lf",vlMP->matrixProfile[l]);
			fprintf(fIndexProfile,"\n%d",vlMP->indexProfile[l]);
			fprintf(flengthProfile,"\n%d",vlMP->lengthProfile[l]);
		}
	}
	
	int nMotifs = (qLengthMax-qLengthMin);
	for(l=0;l<=nMotifs;l++)
	{
		if(l==0)
		{
			fprintf(fEvolution,"%lf",vlMP->distancesEvolutionMotif[l]);
		}
		else
		{
			fprintf(fEvolution,"\n%lf",vlMP->distancesEvolutionMotif[l]);
		}
	}
		
	fclose(fVALMAP);
	fclose(fIndexProfile);
	fclose(flengthProfile);
	fclose(fEvolution);
	
	// top ten matrix profile
	mMxHeap * orderedVALMAP = newHeap(50,false);
	mMxHeap * orderedVALMAPNN = newHeap(50,false);
        
        for(l=0;l<lengthVALMAP;l++)
        {
                pushTopClassic(vlMP->matrixProfile[l], l,vlMP->indexProfile[l], orderedVALMAP);
                pushTopClassic(vlMP->matrixProfileNonLengthNormalized[l], l,vlMP->indexProfileNonLengthNormalized[l], orderedVALMAPNN);  
        }
	
	int rank = 50;
	heapEntry * eVALMAP =  getHead(orderedVALMAP);
        heapEntry * eVALMAPNN =  getHead(orderedVALMAPNN);
	while(eVALMAP!=NULL)
	{
		fprintf(fTopTen,"%d) \t offset1: %d,offset2: %d, normalizedDistance: %lf, length: %d\n", rank, eVALMAP->index1,eVALMAP->index2,eVALMAP->lbDistance,vlMP->lengthProfile[eVALMAP->index1]);
		fprintf(fTopTenNN,"%d) \t offset1: %d,offset2: %d, normalizedDistance: %lf, length: %d\n", rank, eVALMAPNN->index1,eVALMAPNN->index2,eVALMAPNN->lbDistance,vlMP->lengthProfileNonLengthNormalized[eVALMAPNN->index1]);
                removeHead(orderedVALMAP);
                removeHead(orderedVALMAPNN);
		eVALMAP = getHead(orderedVALMAP);
                eVALMAPNN = getHead(orderedVALMAPNN);
		rank--;
	}
	
	freeHeap(orderedVALMAP);
        freeHeap(orderedVALMAPNN);
        
	fclose(fTopTen);
	fclose(fTopTenNN);
        fftw_free(evolutionMotifPair);
	fftw_free(nameFileVALMAP);
	fftw_free(nameFileIndexProfile);
	fftw_free(nameFilelengthProfile);
	fftw_free(nTopTenMatrixProfile);
	fftw_free(nTopTenMatrixProfileNN);
        
	fftw_free(vlMP->matrixProfile);
	fftw_free(vlMP->indexProfile);
	fftw_free(vlMP->lengthProfile);
        
        fftw_free(vlMP->matrixProfileNonLengthNormalized);
	fftw_free(vlMP->indexProfileNonLengthNormalized);
	fftw_free(vlMP->lengthProfileNonLengthNormalized);
	fftw_free(vlMP->distancesEvolutionMotif);
	
	fftw_free(vlMP);

	fftw_free(in);
	fftw_free(in2);
	
	return 1;
}

void VALMOD_RUN(int argc, char ** argv)
{
	time_t t;
    /* Intializes random number generator */
    //srand((unsigned) time(&t));
	int sizeSeries =0;
	int min =0;
	int max =0;
	threshold=1;
	char * fileTS = NULL;
	if(argc>1)
	{
		sizeSeries = atoi(argv[1]);
		min = atoi(argv[2]);
		max = atoi(argv[3]);
		threshold = atoi(argv[5]);
                kDiscords = atoi(argv[6]);
		if(argc>7)
		{
			fileNameValmap = argv[7];
		}
		else
		{
			fileNameValmap = malloc(500);
			snprintf(fileNameValmap,500,"DATA");
		}
		if(argc>8)
		{
			fileTS = argv[8];
		}
		bComparisonStomp = (atoi(argv[4]) == 1);
		bOnlySTOMP = (atoi(argv[4]) == 2);
	}
        
	selfJoinWith_STOMPmultilength_VALMOD(sizeSeries, min, max,fileTS);
}

// fixed length JOIN
void JOIN_RUN(int argc, char ** argv)
{
	time_t t;
    /* Intializes random number generator */
     //srand((unsigned) time(&t));
	int sizeSeriesA =500;
        int sizeSeriesB =200;
	int subLength = 20;
	char * fileTSA = NULL;
        char * fileTSB = NULL;
	
        if(argc> 6)
	{
            sizeSeriesA = atoi(argv[2]);
            sizeSeriesB = atoi(argv[3]);
            subLength = atoi(argv[4]);
            fileTSA = argv[5];
            fileTSB = argv[6];
	}
        else
        {
            sizeSeriesA = atoi(argv[2]);
            sizeSeriesB = atoi(argv[3]);
            subLength = atoi(argv[4]);
        }
        
        if(fileTSA!=NULL)
        {
            fileNameValmap = malloc(500);
            snprintf(fileNameValmap,500,"%s_join_%s",argv[5],argv[6]);
        }
        else
        {
            fileNameValmap = malloc(500);
            snprintf(fileNameValmap,500,"RW_join");
        }
	mpElements * mpEl = _join_STOMP(sizeSeriesA,sizeSeriesB,subLength,fileTSA,fileTSB, true);
        
        free(mpEl->matrixProfile);   
        free(mpEl->indxProfile);  
        free(mpEl->standardProfile);  
        free(mpEl);   
        free(fileNameValmap);   
}

void SELFJOIN_RUN(int argc, char ** argv)
{
	time_t t;
    /* Intializes random number generator */
     //srand((unsigned) time(&t));
	int sizeSeriesA =500;
	int subLength = 20;
	char * fileTSA = NULL;
	
        if(argc> 4)
	{
            sizeSeriesA = atoi(argv[2]);
            subLength = atoi(argv[3]);
            fileTSA = argv[4];
	}
        else
        {
            sizeSeriesA = atoi(argv[2]);
            subLength = atoi(argv[3]);
        }
        
        if(fileTSA!=NULL)
        {
            fileNameValmap = malloc(500);
            snprintf(fileNameValmap,500,"%s_selfjoin",fileTSA);
        }
        else
        {
            fileNameValmap = malloc(500);
            snprintf(fileNameValmap,500,"RW_selfjoin");
        }
	mpElements * mpEl = _selfjoin_STOMP(sizeSeriesA,subLength,fileTSA,true);
        
        free(mpEl->matrixProfile);   
        free(mpEl->indxProfile);   
        free(mpEl->standardProfile);  
        free(mpEl);   
        free(fileNameValmap);   
}







_Bool sanityCheckFile(char * nfileTS)
{
    FILE * fileTS = fopen(nfileTS,"r");
    _Bool bCheckOk = (fileTS!=NULL);
    if(bCheckOk)
    {
        fclose(fileTS);
    }
    return bCheckOk;
}


_Bool sanityCheckLengths(int lengthSeriesA, int lengthSeriesB,int lengthSubsequence)
{
    _Bool bCheckOk = true;
    bCheckOk = bCheckOk || (lengthSubsequence>lengthSeriesB);
    bCheckOk = bCheckOk || (lengthSubsequence>lengthSeriesA);
    bCheckOk = bCheckOk || (lengthSeriesA>0);
    bCheckOk = bCheckOk || (lengthSeriesB>0);
    bCheckOk = bCheckOk || (lengthSubsequence>0);
    
    return bCheckOk;
}

_Bool sanityCheckLengthsSelfJoin(int lengthSeriesA,int lengthSubsequence)
{
    _Bool bCheckOk = true;
    bCheckOk = bCheckOk || (lengthSubsequence>lengthSeriesA);
    bCheckOk = bCheckOk || (lengthSeriesA>0);
    bCheckOk = bCheckOk || (lengthSubsequence>0);
    
    return bCheckOk;
}

mpElements * selfJoinSTOMP(char * fileTSA, int lengthSeriesA, int lengthSubsequence, _Bool bprintMP)
{
   // fileNameValmap = malloc(500);
   // snprintf(fileNameValmap,500,"%s_selfjoin",fileTSA);
    
    mpElements * matrixProfileiProfile = NULL;

    if (sanityCheckLengthsSelfJoin(lengthSeriesA,lengthSubsequence) && sanityCheckFile(fileTSA))
    {
        matrixProfileiProfile = _selfjoin_STOMP(lengthSeriesA,lengthSubsequence,fileTSA,bprintMP);
    }

    return matrixProfileiProfile;
}


mpElements * joinSeriesSTOMP(char * fileTSA, char * fileTSB, int lengthSeriesA, int lengthSeriesB, int lengthSubsequence, _Bool bprintMP)
{
   // fileNameValmap = malloc(500);
   // snprintf(fileNameValmap,500,"%s_join_%s",fileTSA,fileTSB);
    
    mpElements * matrixProfileiProfile = NULL;

    if (sanityCheckLengths(lengthSeriesA, lengthSeriesB,lengthSubsequence) && sanityCheckFile(fileTSA) && sanityCheckFile(fileTSB))
    {
        matrixProfileiProfile = _join_STOMP(lengthSeriesA,lengthSeriesB,lengthSubsequence,fileTSA,fileTSB,bprintMP);
    }
    

    return matrixProfileiProfile;

}

