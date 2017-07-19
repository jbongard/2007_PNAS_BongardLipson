

#include "stdafx.h"

#ifndef _EQ_NODE_H
#define _EQ_NODE_H

#include "matrix.h"

class EQ_NODE {

public:
	int nodeType;
	double index;
	EQ_NODE *parent;
	EQ_NODE *left;
	EQ_NODE *right;
	double minVal;
	double maxVal;

public:
	EQ_NODE(int depth);
	EQ_NODE(int nt, double i);
	EQ_NODE(EQ_NODE *parentEq);
	EQ_NODE(char *s,int start, int end);
	~EQ_NODE(void);
	void   AllCalculateTermDegrees(int currTerm, MATRIX *termDegrees);
	MATRIX *AllExtractTermsAndDegrees(void);
	void   AllFindTerms(int *numTerms, EQ_NODE **terms);
	int    AllMultsBelow(void);
	void   AllSplit(void);
	void   Clean(void);
	int    Count(int currCount);
	double Evaluate(double currT, MATRIX *y, MATRIX *params, int *valid);
	EQ_NODE *Find(int *counter);
	int    GetDepth(void);
	int    IsEqualTo(EQ_NODE *otherTree);
	void   Mutate(void);
	void   MutateNeutrally(void);
	void   Print(void);
	void   PrintInfix(void);
	void   ResetNodeRanges(void);
	void   Save(ofstream *outFile);

private:
	void   AddANewTerm(void);
	void   CollapseMult(void);
	void   CollapseParameters(void);
	void   CollapsePlus(void);
	void   CollapseSubtraction(void);
	void   ExpandParam(void);
	void   ExtractNodeOp(char *s, int start);
	void   InheritParams(EQ_NODE *donatingNode);
	int	   IsParam(void);
	int	   IsPlusMinusMult(void);
	int    IsTerminal(void);
	int	   IsTwoParity(void);
	int	   IsVar(void);
	int	   IsZeroParity(void);
	void   MutateCreep(void);
	void   MutateCrunch(void);
	void   MutateDumbSnip(void);
	void   MutateGrow(void);
	void   MutateMerge(void);
	void   MutateNeutrallyLeft(void);
	void   MutateNeutrallyRight(void);
	void   MutateNode(void);
	void   MutateOneParity(void);
	void   MutateSnip(void);
	void   MutateSplit(void);
	void   MutateSquish(void);
	void   MutateSwapOperators(void);
	void   MutateSwapSubtrees(void);
	void   MutateTerminal(int allowNudge);
	void   MutateTwoParity(void);
	void   MutateVal(void);
	void   MutateZeroParity(void);
	void   PrintNode(void);
	void   SaveNode(ofstream *outFile);
};

#endif