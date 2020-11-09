#ifndef FiboKLAMD_H
#define FiboKLAMD_H

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
}
#endif

// the heap entry
typedef struct entry{
	double distance;
	double lbSize;
	double lbDistance;
	double stdDevQuery;
	int index1;
	int index2;
	double dP;
	double s1;
	double s2;
	double ss1;
	double ss2;
	int info;
} heapEntry;

// the heap
typedef struct heap{
	heapEntry * heap_;
	unsigned long long size;
	unsigned long long  maxElements;
	_Bool bMin;
	int position;
} mMxHeap;
int pushTopClassic(double distance, int index1, int index2, mMxHeap * heap);
int push(double distance, int index1, int index2, double dP, double s1, double s2, double ss1, double ss2,double stD1, double stD2, int length,mMxHeap * heap);
mMxHeap * newHeap(unsigned long long  maxElements, _Bool bMin); // if bMin true is a min heap otherwise is a max heap
void freeHeap(mMxHeap * heap);
heapEntry * popHead(mMxHeap * heap);
heapEntry * getHead(mMxHeap * heap);
void removeHead(mMxHeap * heap);
#endif /* FiboKLAMD_H */