/* ---------------------------------------------------
   FILE:     bodyPlan.h
	AUTHOR:   Josh Bongard
	DATE:     October 2, 2000
	FUNCTION: This class contains all information for
				 a single physical segment of an organism
				 in the MathEngine environment.
 -------------------------------------------------- */

#include "stdafx.h"
#include "fstream.h"

#ifndef _MATRIX_H
#define _MATRIX_H

class MATRIX {

public:
	int length;
	int width;
	double *vals;

public:
	MATRIX(int ln, int wd);
	MATRIX(int ln, int wd, double val);
	MATRIX(MATRIX *m);
	MATRIX(ifstream *inFile);
	MATRIX(char *fileName);
	~MATRIX(void);
	double   AbsoluteDifference(MATRIX *m);
	double   AbsoluteDifference(MATRIX *m, int startWd, int endWd);
	double   AbsoluteDifference(MATRIX *m, int wd);
	double   AbsoluteDifference(MATRIX *m, int wd, double normFactor);
	double   AbsoluteDifferenceOfTails(MATRIX *m);
	void	 Add(int i, int j, double val);
	void     Add(MATRIX *him);
	void     AddNoise(double amt);
	int		 AllRowsAlsoIn(MATRIX *other);
	int      BinaryToDecimal(int i, int j, int ln);
	void     ChangePairInRow(int i);
	int		 ComputeFSM(MATRIX *str, MATRIX *states, int i, int strIndex, int strLength, int q0, MATRIX *F);
	void	 ComputeVariances(MATRIX *vars);
	void     Copy(MATRIX *him);
	void     CopyRow(int myI, MATRIX *him, int hisI);
	MATRIX  *CountValues(int i1, int i2, int j1, int j2, int maxVal);
	void	 CreateChain(int depth);
	void	 CreateIdentity(void);
	void	 CreateParity(void);
	void	 CreatePermutation(int min, int max);
	void	 CrossRow(int myRow, MATRIX *him, int hisRow);
	void	 CrossSubTrees(int i,MATRIX *him);
	void	 DecreaseFullColumns(void);
	int      Equal(MATRIX *other);
	double   EqualColumnVals(int col1, int col2);
	int      Exploding(void);
	double   Evaluate(int i, double currT, MATRIX *y, MATRIX *params);
	void	 FillColumn(int myCol, int hisCol, MATRIX *m);
	void     FillRow(int myRow, int hisRow, MATRIX *m);
	void     Flip(int i, int j);
	void	 FlipRandomBit(void);
	void	 FlipRandomBitInRow(int i);
	double   Get(int i, int j);
	MATRIX  *GetColumn(int j);
	void     GetMaxValsInColumn(int j, MATRIX *maxVals);
	void     GetMaxValsInColumn(int j, int iFirst, int iLast, MATRIX *maxVals);
	MATRIX  *GetRow(int i);
	int		 In(double val);
	void     IncreaseEmptyColumns(void);
	void	 InitColumns(int colSum);
	void     InsertValSomewhereInEachRow(double val);
	double   MaxDiffCol(int myCol, MATRIX *him, int hisCol);
	double   MaxDiffColChange(int myCol, MATRIX *him, int hisCol);
	double   MaxDiff(MATRIX *other);
	double	 MaxValInColumn(int j);
	double   Mean(void);
	double   MeanOfColumn(int j);
	double   MinValInColumn(int j);
	int		 MostSimilarRow(MATRIX *r, int i1, int i2);
	double	 MSQColumn(int myCol, MATRIX *him, int hisCol);
	double   MSQRow(int myRow, MATRIX *him, int hisRow);
	void	 Mult(double val);
	void     Mutate(double min, double max);
	void	 Nudge(double min, double max);
	void	 Perturb(double min, double max);
	void	 Print(int wait);
	void     Randomize(int maxVal);
	void     RandomizeColumn(int j, int maxVal);
	void     RandomizeRow(int i, int maxVal);
	void	 Replace(MATRIX *m);
	double	 RollingMean(MATRIX *other, int h, int w);
	int		 RowDifference(int myRow, int hisRow, MATRIX *m);
	void     Save(ofstream *outFile);
	void	 SelectUniquelyFrom(int maxVal);
	void	 Set(int i, int j, double val);
	void	 SetAllTo(int val);
	void     SetRow(int i, int val);
	double	 SumOfColumn(int j);
	double	 SumOfColumn(int i1, int i2, int j);
	double   SumOfIndices(MATRIX *indices,int i, int j1, int j2);
	double	 SumOfRow(int i);
	void     SwapRow(int i, MATRIX *him);
	double   VarOfColumn(int j);
	void	 Write(ofstream *outFile);
	void	 WriteAndRename(char *fileName);
	void	 ZeroColumn(int j);
	void	 ZeroIt(void);

private:
	int		BinaryToDecimal(int j);
	int	    Contains(int val);
	MATRIX *DecimalToBinary(int val, int maxValue);
	double  EvaluateNonTerminal(int i, double currT, MATRIX *y, MATRIX *params);
	double  EvaluateTerminal(int i, double currT, MATRIX *y, MATRIX *params);
	int		FindFirstValue(int j);
	double  Mean(int i, int startJ, int endJ);
	void    NudgePairInRow(int i);
	void	InitColumn(int j, int colSum);
	void	InitRandomly(void);
	void    ReadFromFile(ifstream *inFile);
	void    ReplacePairInRow(int i);
	double  RollingMean(MATRIX *other, int j, int h, int w);
};

#endif