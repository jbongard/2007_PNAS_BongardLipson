
#include "stdafx.h"

#ifndef _EXP_PHASE_H
#define _EXP_PHASE_H

#include "model.h"
#include "test.h"

class EXP_PHASE {

public:
	int numTests;
	int numModels;
	TEST **tests;
	TEST **children;

	MODEL **models;

	RESULT **resultMatrix;

public:

	EXP_PHASE(MODEL **m);
	~EXP_PHASE(void);
	TEST *Evolve(int numT1, TEST **t1, int numT2, TEST **t2);

private:
	void ComputeFinalValues(RESULT **results, MATRIX *finalValues);
	void CreateChildPopulation(void);
	void CreateHorizontalResults(RESULT **results, TEST **testsToEvaluate, int i);
	void CreateNewClimbers(void);
	void CreateResultMatrix(int numPreviousTests, TEST **previousTests);
	void DestroyResults(RESULT **results);
	void DestroyResultMatrix(int numPreviousTests);
	int  Done(void);
	void Evaluate(int useChildren, int numT1, TEST **t1, int numT2, TEST **t2);
	int  ParentChildSimilarity(int index);
	void Print(void);
	void ReplaceClimbers(void);
	void ReplaceParentsWithChildren(void);
	void SeekFinalDisagreement(int i, TEST **testsToEvaluate, MATRIX *finalValues, MATRIX *variances, 
		                       int numPreviousTests, TEST **previousTests);
	void SeekPastDisagreement(int i, TEST **testsToEvaluate, RESULT **results, int numPreviousTests); 
	void SeekTotalDisagreement(int i, TEST **testsToEvaluate, RESULT **results,
								int numPreviousTests, TEST **previousTests,
								int numBankedTests,   TEST **bankedTests);
	void SeekExplosions(int i, TEST **testsToEvaluate, MATRIX *finalValues);
	void ShuffleTests(void);
	void SortTests(void);
	TEST **WhichTestsToUse(int useChildren);

};

#endif