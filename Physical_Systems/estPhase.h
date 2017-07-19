
#include "stdafx.h"

#ifndef _EST_PHASE_H
#define _EST_PHASE_H

#include "model.h"

class EST_PHASE {

public:
	MODEL **models;
	MODEL **children;

	int numOfTests;
	TEST **tests;
	RESULT **results;
	int numModels;
	int currCycle;
	int stagnationCounter;
	double lastFitness;

private:
	ofstream *fitFile;
	ofstream *depthFile;
	ofstream *modelFile;
	ofstream *sizeFile;

public:

	EST_PHASE(MODEL *target, int c, int testsSoFar, TEST **t, RESULT **r);
	EST_PHASE(MODEL **m, int c, int testsSoFar, TEST **t, RESULT **r);
	~EST_PHASE(void);
	MODEL **Evolve(void);
	double GetBestFitness(void);
	void Print(int toModel);

private:
	void CloseDataFiles(void);
	void CreateChildPopulation(void);
	void CreateNewClimbers(void);
	int  Done(void);
	void Evaluate(int evaluateParents);
	int  FitnessExists(double fit);
	void InitModels(MODEL *target);
	void OpenDataFiles(int append);
	void OpenDepthDataFile(int append);
	void OpenFitnessDataFile(int append);
	void OpenModelDataFile(int append);
	void OpenSizeDataFile(int append);
	int  ParentChildSimilarity(int index);
	void PrintBestModel(void);
	void PrintProgress(void);
	void RecoverEqs(void);
	void ReplaceClimbers(void);
	void ReplaceParentsWithChildren(void);
	void RestoreEqs(void);
	void ShuffleModels(void);
	void SaveData(void);
	void SaveDepths(void);
	void SaveFits(void);
	void SaveHits(void);
	void SaveModels(void);
	void SaveSizes(void);
	void SortModels(void);
	void StoreAllEqs(void);
	void StoreEqs(void);
};

#endif