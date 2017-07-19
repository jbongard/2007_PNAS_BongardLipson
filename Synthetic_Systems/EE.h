

#include "stdafx.h"

#ifndef _EE_H
#define _EE_H

#include "estPhase.h"
#include "expPhase.h"

#include "model.h"
#include "test.h"
#include "result.h"

class EE {

public:
	MODEL  *target;
	MODEL  **models;
	TEST   **tests;
	TEST   **bankedTests;
	int    *ages;
	RESULT **bankedResults;
	int    numBankedTests;
	int    numTests;
	RESULT **results;
	int	   iterations;
	int    currCycle;
	EST_PHASE *estimationPhase;
	EXP_PHASE *explorationPhase;

private:
	ofstream *testFile;
	ofstream *realFitFile;
	ofstream *bankedFile;

public:

	EE(void);
	~EE(void);
	void PerformBatch(void);
	void PerformInference(void);

private:
	void AgeBankedTests(int append);
	void BackupModels(void);
	void DepositTest(void);
	void DestroyModels(void);
	void DestroyResults(void);
	void DestroyTarget(void);
	void DestroyTests(void);
	void DoInference(void);
	void Estimate(int currIteration);
	void Explore(int currIteration);
	void HalveTestSearchEffort(void);
	void InitResults(void);
	void InitTests(void);
	void LoadTarget(void);
	void Print(void);
	void RestoreOldModels(void);
	void SaveRealFitness(int append);
	void SaveTest(TEST *test, int append);
	int  TestTooHard(void);
	int  WithdrawTest(void);
};

#endif