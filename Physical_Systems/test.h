

#include "stdafx.h"

#ifndef _TEST_H
#define _TEST_H

#include "matrix.h"

class TEST {

public:
	MATRIX *y0;
	double fitness;
	int chosenIndexIntoData;

public:

	TEST(MATRIX *targetData);
	TEST(MATRIX *targetData, int startingIndex);
	TEST(TEST *otherTest);
	~TEST(void);
	double ClosestTestDistance(int numOtherTests,TEST **otherTests);
	void Cross(TEST *otherTest);
	void Mutate(MATRIX *targetData);
	void Print(int wait);
	void Save(ofstream *outFile);
	double Similarity(TEST *other);
	int  TestExists(int numT, TEST **t);
	void ZeroIt(void);

private:
	int  Equal(TEST *other);
};

#endif