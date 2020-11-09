#######################################################################
#######################################################################
#####				MATRIX PROFILE API BACKEND					 ######
#######################################################################
#######################################################################


from ctypes import *
import os
import imp
############### WRAP C CODE ###############
class mpElements(Structure):
	pass

mpElements._fields_ = [("matrixProfile", POINTER(c_double)),("indxProfile",  POINTER(c_int)),("standardProfile",  POINTER(c_double))]

_mp = cdll.LoadLibrary(imp.find_module('valmod')[1])

_mp.selfJoin.argtypes = [POINTER(c_char),c_int,c_int,c_bool]
_mp.selfJoin.restype  = POINTER(mpElements)

_mp.join.argtypes = [POINTER(c_char),POINTER(c_char),c_int,c_int,c_int,c_bool]
_mp.join.restype  = POINTER(mpElements)

_mp.freempElements.argtypes = [POINTER(mpElements)]
_mp.freempElements.restype = c_int
###########################################



def selfJoin(nameFileSeries, lengthSeries, subsequenceLength):
	"""
		INPUT
			nameFileSeries: path of the file containing the time series data to analize
							(the time series file should contain a point for each line)
		    lengthSeries: the total number of point contained in the file
			subsequenceLength: the length of the sliding window used to compute the matrix profile

		OUTPUT
		 	tuple(matrix profile, index profile)
			matrix profile: is a list containg the distances between the subsequence at that position and its nearest neigbhor subsequence in the overall sequence
			index profile: is a list containg the index of subsequence with minimum distance for the subsequence at that position
	"""
	return Join(nameFileSeries, nameFileSeries, lengthSeries, lengthSeries, subsequenceLength,True)


def Join(nameFileSeriesA, nameFileSeriesB, lengthSeriesA, lengthSeriesB, subsequenceLength,selfJoin=False):
	"""
		INPUT
			nameFileSeriesA/B: path of the file containing the time series data to analize
							(the time series file should contain a point for each line)
		    lengthSeriesA/B: the total number of point contained in the file
			subsequenceLength: the length of the sliding window used to compute the matrix profile

		OUTPUT
		 	tuple(matrix profile, index profile)
			matrix profile: is a list containg the distances between the subsequence B at that position and its nearest neigbhor subsequence in A
			index profile: is a list containg the index of subsequence with minimum distance for the subsequence at that position
	"""
	global _mp
	matrixProfile, indexProfile = [], []

	bytesNameFileSeriesA= nameFileSeriesA.encode('utf-8')
	bytesNameFileSeriesB= nameFileSeriesB.encode('utf-8')

	if selfJoin:
		result = _mp.selfJoin(bytesNameFileSeriesA,c_int(lengthSeriesA),c_int(subsequenceLength),c_bool(False))
	else:
		result = _mp.join(bytesNameFileSeriesA, bytesNameFileSeriesB, c_int(lengthSeriesA), c_int(lengthSeriesB), c_int(subsequenceLength), c_bool(False))

	if result:
		mp_len = range(((lengthSeriesB-subsequenceLength)+1))
		c = result.contents
		matrixProfile = [c.matrixProfile[i] for i in mp_len]
		indexProfile = [c.indxProfile[i] for i in  mp_len]
		_mp.freempElements(result)
	return matrixProfile,indexProfile
