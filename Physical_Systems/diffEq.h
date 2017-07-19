/* ---------------------------------------------------
   FILE:     bodyPlan.h
	AUTHOR:   Josh Bongard
	DATE:     October 2, 2000
	FUNCTION: This class contains all information for
				 a single physical segment of an organism
				 in the MathEngine environment.
 -------------------------------------------------- */

#include "stdafx.h"

#ifndef _DIFF_EQ_H
#define _DIFF_EQ_H

#include "eqNode.h"
#include "matrix.h"

class DIFF_EQ {

public:
	EQ_NODE *root;
	MATRIX  *softTree;

public:
	DIFF_EQ(void);
	DIFF_EQ(DIFF_EQ *parent);
	DIFF_EQ(ifstream *inFile);
	~DIFF_EQ(void);
	void Clean(void);
	int  Count(void);
	void Cross(DIFF_EQ *sister);
	double Evaluate(double currT, MATRIX *y, MATRIX *params, int *valid);
	int  GetMaxDepth(void);
	int  IsAHit(int currEq);
	void Mutate(void);
	void Print(void);
	void PrintInfix(void);
	void Save(ofstream *outFile);
	double Similarity(DIFF_EQ *other);

private:
	void CreateRandom(void);
	void CrossRoot(DIFF_EQ *sister);
	void CrossSoftTree(DIFF_EQ *sister);
	void MutateNeutrally(void);
	void MutateSoftTree(void);
	void MutateRoot(void);
};

#endif