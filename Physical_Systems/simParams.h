/* ---------------------------------------------------
   FILE:     simParams.h
	AUTHOR:   Josh Bongard
	DATE:     October 5, 2000
	FUNCTION: This class contains all miscellaneous
				 data and functions for this simulation.
 -------------------------------------------------- */

#include "stdafx.h"
#include "fstream.h"
#include "time.h"

#ifndef _SIM_PARAMS_H
#define _SIM_PARAMS_H

#include "eqNode.h"
#include "matrix.h"
#include "model.h"

class SIM_PARAMS {

public:
	int		 randSeed;
	int		 performingEstimation;
	int		 evolvingTests;
	int		 useTheBank;
	int		 evolveEqs;
	int		 evolveParams;
	MATRIX  *tempY;
	MATRIX  *kn1;
	MATRIX  *kn2;
	MATRIX  *kn3;
	MATRIX  *kn4;
	char	targetFileName[100];
	int		numParams;
	int     numEqs;
	int		lengthOfExperiment;
	int		targetDepth;
	int		writeToScreen;
	int		softGP;
	int		maxCycles;
	int		modelEvals;
	int		targetEvals;
	EQ_NODE **storedEqs;
	MATRIX	*previousFits;
	MATRIX  *storedFits;
	MATRIX  *storedDepths;
	MATRIX  *storedSizes;
	int		currGen;
	int		evaluatingTarget;
	MATRIX *currY;
	MATRIX *nextY;
	double  testMin;
	double  testMax;
	int		snip;
	int		creep;
	int		worstTest;
	time_t	startTime;
	EQ_NODE **previousEqs;
	int		expPhaseMaxGens;
	int		dumbSnip;
	int		prune;
	int		testDist;
	int		mutateCrunch;
	int		mutateSplit;
	int		mutateMerge;
	int		mutateNeutrally;
	int		estPhaseMaxGens;
	int		fullDis;
	int		resultDist;
	int		addNoise;
	double	noiseAmt;
	MATRIX  *minVarValues;
	MATRIX  *maxVarValues;
	int		fullTestDist;
	MATRIX **targetMatrices;
	int		useTimer;
	int		maxMinutes;
	int		dataPoints;

public:	
	SIM_PARAMS(int argc, char **argv);
	~SIM_PARAMS(void);
	void   CreateTargetMatrices(MODEL *target);
	int	   DiffEqUnderScrutiny(void);
	void   FileDelete(char *fileName);
	int    FileExists(char *fileName);
	void   FileRename(char *src, char *dest);
	int    FindStrEnd(char *s, int pos, int count);
	int    FlipCoin(void);
	int    GenOfChange(void);
	int    GenOfChangeNext(void);
	int    GetExperimentRegime(void);
	ofstream *GetOutFile(char *fileName);
	void   ParseParameters(int argc, char **argv);
	double Rand(double min, double max);
	int    RandInt(int min, int max);
	char  *StrCpy(char *s, int start, int end);
	void   TuneNoiseByData(void);

private:
	void   WriteEndFile(void);
};

#endif