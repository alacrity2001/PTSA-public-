
#include <stdlib.h>
#include <stdio.h> 
#include <string.h>
#include <stdbool.h>
#include "FiboKLAMD.h"
#include <math.h>

void sinkDown(int k, mMxHeap * heap){
	
	if (heap->bMin)
	{
		int smallest =k;
		if(2*k<heap->position && heap->heap_[smallest].lbDistance > heap->heap_[2*k].lbDistance  && heap->heap_[smallest].distance >=0 && heap->heap_[2*k].distance>=0){
			smallest = 2*k;
		}
		if(2*k+1 <heap->position &&  heap->heap_[smallest].lbDistance > heap->heap_[2*k+1].lbDistance && heap->heap_[smallest].distance >=0 &&  heap->heap_[2*k+1].distance >=0){
			smallest = 2*k+1;
		}
		if(smallest!=k){
			heapEntry temp = heap->heap_[k];
			heap->heap_[k] = heap->heap_[smallest];
			heap->heap_[smallest] = temp;
			sinkDown(smallest,heap);
		}
	}
	else
	{
		int greater =k;
		if(2*k < heap->position && heap->heap_[greater].lbDistance < heap->heap_[2*k].lbDistance && heap->heap_[greater].distance >=0 && heap->heap_[2*k].distance>=0){
			greater = 2*k;
		}
		if(2*k+1 < heap->position &&  heap->heap_[greater].lbDistance < heap->heap_[2*k+1].lbDistance && heap->heap_[greater].distance >=0 && heap->heap_[2*k+1].distance >=0){
			greater = 2*k+1;
		}
		if(greater!=k){
			heapEntry temp = heap->heap_[k];
			heap->heap_[k] = heap->heap_[greater];
			heap->heap_[greater] = temp;
			sinkDown(greater,heap);
		}
		
	}
}


void bubbleUp(mMxHeap * heap)
{
	if(heap->bMin)
	{
		int pos = heap->position-1;
		while(pos>0 && heap->heap_[pos/2].lbDistance > heap->heap_[pos].lbDistance && heap->heap_[pos/2].distance >=0 && heap->heap_[pos].distance >=0 ){
			heapEntry y = heap->heap_[pos];
			heap->heap_[pos]=heap->heap_[pos/2];
			heap->heap_[pos/2] = y;
			pos = pos/2;
		}
	}
	else
	{
		int pos = heap->position-1;
		while(pos>0 && heap->heap_[pos/2].lbDistance < heap->heap_[pos].lbDistance && heap->heap_[pos/2].distance >=0 && heap->heap_[pos].distance >=0){
			heapEntry y = heap->heap_[pos];
			heap->heap_[pos]=heap->heap_[pos/2];
			heap->heap_[pos/2] = y;
			pos = pos/2;
		}
	}
}

int push(double distance, int index1, int index2, double dP, double s1, double s2, double ss1, double ss2,double stD1, double stD2, int length, mMxHeap * heap)
{

	int info = -1; // if we replace the max/min, it will be the pos of the old max/min
	double lbDistance;
	double q = ((dP/length)-((s1/length)*(s2/length)))/(stD1*stD2);
	if(q>0)
	{
		lbDistance = length*(1-(q*q));
	}
	else
	{
		lbDistance = length;
	}
	
	if (heap->size < heap->maxElements)
	{
		info = heap->size;
		heap->size=heap->size+1;
	}
	else
	{
		// VALMOD algorithm specific
		heapEntry * head = getHead(heap);
		if(!heap->bMin)
		{
			if(lbDistance < head->lbDistance){
			
				info = head->info;
				removeHead(heap);
			}
			else
				return -1;
		}
		else
		{
			if(lbDistance > head->lbDistance){
				info = head->info;
				removeHead(heap);
			}
			else
				return -1;
		}
	}
	
	int insertActual; 
	if(heap->position==0)
	{
		insertActual = heap->position+1;
		heap->position = 2;
	}
	else
	{
		insertActual = heap->position;
		heap->position=heap->position+1;
		
	}
	
	
	heap->heap_[insertActual].lbDistance = lbDistance;
	heap->heap_[insertActual].stdDevQuery = stD1;
	heap->heap_[insertActual].distance = distance;
	heap->heap_[insertActual].index1 = index1;
	heap->heap_[insertActual].index2 = index2;
	heap->heap_[insertActual].dP = dP;
	heap->heap_[insertActual].s1 = s1;
	heap->heap_[insertActual].s2 = s2;
	heap->heap_[insertActual].ss1 = ss1;
	heap->heap_[insertActual].ss2 = ss2;
	heap->heap_[insertActual].info = info;
	
	if (heap->position>2)
	{
		bubbleUp(heap);
	}

	return info;
}

int pushTopClassic(double distance, int index1, int index2, mMxHeap * heap)
{
	int info =0;
	if (heap->size < heap->maxElements)
	{
		info = heap->size;
		heap->size=heap->size+1;
	}
	else
	{
		// KLAMD algorithm specific
		heapEntry * head = getHead(heap);
		if(!heap->bMin)
		{
			if(distance < head->lbDistance){
			
				info = head->info;
				removeHead(heap);
			}
			else
				return -1;
		}
		else
		{
			if(distance > head->lbDistance){
				info = head->info;
				removeHead(heap);
			}
			else
				return -1;
		}
	}
	
	int insertActual; 
	if(heap->position==0)
	{
		insertActual = heap->position+1;
		heap->position = 2;
	}
	else
	{
		insertActual = heap->position;
		heap->position=heap->position+1;
		
	}
	heap->heap_[insertActual].distance = distance;
	heap->heap_[insertActual].lbDistance = distance;
	heap->heap_[insertActual].index1 = index1;
	heap->heap_[insertActual].index2 = index2;
	
	if (heap->position>2)
	{
		bubbleUp(heap);
	}

	return info;
}




heapEntry * popHead(mMxHeap * heap)
{
	//TBD
	return NULL;
}

void removeHead(mMxHeap * heap)
{
    
        heap->heap_[1]=heap->heap_[heap->position-1];
        heap->heap_[heap->position-1].distance = -1; // it will be discarded, hopefully
        heap->position= heap->position - 1;		
        sinkDown(1,heap);
        
}


heapEntry * getHead(mMxHeap * heap)
{
	if(heap->position>=2)
	{
		return &heap->heap_[1];
	}
	else
	{
		return NULL;
	}
}


mMxHeap * newHeap(unsigned long long  maxElements,_Bool bMin)
{
	mMxHeap * heap = malloc(sizeof(mMxHeap));
	heap->heap_ =  malloc(sizeof(heapEntry)*(maxElements+1));
	heap->heap_[0].distance =-1;
	heap->size = 0;
	heap->position = 0;
	heap->maxElements = maxElements;
	heap->bMin = bMin;
        
	return heap;
}

void freeHeap(mMxHeap * heap)
{
	free(heap->heap_);
	free(heap);
	heap = NULL;
}




