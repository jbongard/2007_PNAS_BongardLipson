/* ---------------------------------------------------
   FILE:     bodyPlan.h
	AUTHOR:   Josh Bongard
	DATE:     October 2, 2000
	FUNCTION: This class contains all information for
				 a single physical segment of an organism
				 in the MathEngine environment.
 -------------------------------------------------- */

#include "stdafx.h"

#ifndef _MODEL_H
#define _MODEL_H

#include "diffEq.h"
#include "result.h"
#include "test.h"

class MODEL {

public:
	MATRIX *params;
	DIFF_EQ **diffEqs;
	double fitness;
	double parentFitness;
	int exploding;

public:
	MODEL(MODEL *parent, int copyParams, int copyEqs);
	MODEL(ifstream *inFile);
	~MODEL(void);
	void Cross(MODEL *sister);
	void Evaluate(TEST *test, RESULT *result, int lengthOfTest);
	void Evaluate(TEST *test, RESULT *result, int lengthOfTest, RESULT *targetResult);
	int  GetMaxDepth(void);
	void Mutate(void);
	void Print(void);
	void PrintInfix(ofstream *outStream);
	void Save(ofstream *outFile);
	double Similarity(MODEL *other);
	int  Size(void);
	int  WorseThanParent(void);

private:
	void CreateDiffEqs(void);
	void CreateRandomParams(void);
	void CrossEqs(MODEL *sister);
	void EvaluateDiffEqs(double currT, MATRIX *currY, int *valid);
	void EvaluateDiffEqs(int i, double currT, MATRIX *currY, MATRIX *nextY, int *valid);
	void MutateEqs(void);
	double SimilarityBetEqs(MODEL *other);
};

#endif