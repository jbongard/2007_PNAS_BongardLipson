
#include "stdafx.h"

#ifndef _RESULT_H
#define _RESULT_H

#include "matrix.h"

class RESULT {

public:
	MATRIX *outputData;

public:

	RESULT(void);
	~RESULT(void);
	void   AddNoise(void);
	double ComputeDifference(RESULT *otherResult);
	int    GetLengthOfTest(void);
	void   Print(void);
	void   UpdateVarRanges(void);
};

#endif