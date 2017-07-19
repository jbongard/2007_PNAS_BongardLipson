

#include "stdafx.h"
#include "stdlib.h"

#include "time.h"

#ifndef _EE_CPP
#define _EE_CPP

#include "constants.h"
#include "EE.h"
#include "simParams.h"

extern SIM_PARAMS *simParams;

EE::EE(void) {

	LoadTarget();

	InitTests();

	InitResults();

	models = NULL;
}

EE::~EE(void) {

	DestroyModels();

	DestroyResults();

	DestroyTests();

	DestroyTarget();
}

void EE::PerformBatch(void) {

	for (currCycle=0;currCycle<simParams->maxCycles;currCycle++) {

		tests[currCycle] = new TEST(target->targetData);

		if ( currCycle == 0 )
			SaveTest(tests[currCycle],false);
		else
			SaveTest(tests[currCycle],true);

		results[currCycle] = new RESULT;

		simParams->evaluatingTarget = true;
		target->Evaluate(tests[currCycle],results[currCycle],simParams->lengthOfExperiment);
		simParams->evaluatingTarget = false;
		simParams->targetEvals++;
		if ( simParams->addNoise )
			results[currCycle]->AddNoise();
		results[currCycle]->UpdateVarRanges();
		results[currCycle]->Print();
	}

	numTests = simParams->maxCycles;

	Estimate(simParams->maxCycles-1);
}

void EE::PerformInference(void) {

	currCycle = 0;

	if ( simParams->evolvingTests )
		tests[numTests++] = new TEST(target->targetData);
	else
		tests[numTests++] = new TEST(target->targetData,0);

	SaveTest(tests[numTests-1],false);

	results[numTests-1] = new RESULT;

	simParams->evaluatingTarget = true;
	target->Evaluate(tests[numTests-1],results[numTests-1],simParams->lengthOfExperiment);
	simParams->evaluatingTarget = false;
	simParams->targetEvals++;
	results[numTests-1]->UpdateVarRanges();
	if ( simParams->addNoise ) {
		//simParams->TuneNoiseByData();
		results[numTests-1]->AddNoise();
	}
	currCycle++;

	if ( !simParams->useTimer ) {

		while ( currCycle<simParams->maxCycles )

			DoInference();
	}
	else {

		time_t currTime = time(NULL);
		double minutesSinceStart = (currTime - simParams->startTime)/60.0;

		while ( minutesSinceStart < simParams->maxMinutes ) {

			DoInference();

			currTime = time(NULL);
			minutesSinceStart = (currTime - simParams->startTime)/60.0;
		}
	}
}

// ----------------- Private Functions ------------------

void EE::AgeBankedTests(int append) {

	RESULT *tempResult;
	double diff;
	char fileName[50];


	sprintf(fileName,"data/run%d_%d_banked.dat",simParams->randSeed,simParams->GetExperimentRegime());

	if ( append )
		bankedFile = new ofstream(fileName,ios::app);
	else
		bankedFile = new ofstream(fileName);

	(*bankedFile) << currCycle << " " << numBankedTests << " ";

	for (int i=0;i<numBankedTests;i++) {

		ages[i]++;

		tempResult = new RESULT;

		models[0]->Evaluate(bankedTests[i],tempResult,bankedResults[i]->GetLengthOfTest());

		diff = bankedResults[i]->ComputeDifference(tempResult);

		if ( simParams->writeToScreen )
			printf("%d %3.3f\n",ages[i],diff);

		(*bankedFile) << diff << " ";

		delete tempResult;
		tempResult = NULL;
	}

	(*bankedFile) << "\n";

	bankedFile->close();
	delete bankedFile;
	bankedFile = NULL;
}

void EE::BackupModels(void) {

	if ( !simParams->previousEqs ) {
		simParams->previousEqs = new EQ_NODE * [simParams->numEqs*NUM_MODELS];

		for (int i=0;i<simParams->numEqs*NUM_MODELS;i++)

			simParams->previousEqs[i] = NULL;

	}

	int k=0;

	for (int i=0;i<NUM_MODELS;i++)

		for (int j=0;j<simParams->numEqs;j++) {

			if ( !simParams->previousEqs[k] )
				delete simParams->previousEqs[k];

			simParams->previousEqs[k] = new EQ_NODE(models[i]->diffEqs[j]->root);
			simParams->previousEqs[k]->ResetNodeRanges();
			k++;
		}
}

void EE::DepositTest(void) {

	bankedTests[numBankedTests] = tests[numTests-1];
	tests[numTests-1] = NULL;

	ages[numBankedTests] = 0;

	bankedResults[numBankedTests] = results[numTests-1];
	results[numTests-1] = NULL;

	numTests--;
	numBankedTests++;
}

void EE::DestroyModels(void) {

	if ( models ) {

		for (int i=0;i<NUM_MODELS;i++) {

			if ( models[i] ) {
				delete models[i];
				models[i] = NULL;
			}
		}

		delete[] models;
		models = NULL;
	}
}

void EE::DestroyResults(void) {

	for (int i=0;i<numTests;i++) {

		if ( results[i] ) {
			delete results[i];
			results[i] = NULL;
		}
	}

	delete[] results;
	results = NULL;

	if ( simParams->useTheBank ) {

		for (i=numBankedTests;i>0;i--) {

			delete bankedResults[i-1];
			bankedResults[i-1] = NULL;
		}

		delete[] bankedResults;
		bankedResults = NULL;
	}
}

void EE::DestroyTarget(void) {

	delete target;
	target = NULL;
}

void EE::DestroyTests(void) {

	for (int i=0;i<numTests;i++) {

		if ( tests[i] ) {
			delete tests[i];
			tests[i] = NULL;
		}
	}

	delete[] tests;
	tests = NULL;

	if ( simParams->useTheBank ) {

		for (i=numBankedTests;i>0;i--) {

			delete bankedTests[i-1];
			bankedTests[i-1] = NULL;
		}

		delete[] bankedTests;
		bankedTests = NULL;

		delete[] ages;
		ages = NULL;
	}
}

void EE::DoInference(void) {

	Estimate(currCycle);

	if ( currCycle<(simParams->maxCycles-1) ) {

		if ( simParams->useTheBank && WithdrawTest() );

		else {
			Explore(currCycle);

			results[numTests-1] = new RESULT;

			simParams->evaluatingTarget = true;
			target->Evaluate(tests[numTests-1],results[numTests-1],simParams->lengthOfExperiment);
			simParams->evaluatingTarget = false;
			simParams->targetEvals++;
			if ( simParams->addNoise )
				results[numTests-1]->AddNoise();
			results[numTests-1]->UpdateVarRanges();
		}

		if ( simParams->useTheBank ) {

			if ( currCycle == 1 )
				AgeBankedTests(false);
			else
				AgeBankedTests(true);
		}
	}

	currCycle++;
}

void EE::Estimate(int currIteration) {

	if ( models )
		estimationPhase = new EST_PHASE(models,currIteration,numTests,tests,results);
	else
		estimationPhase = new EST_PHASE(target,currIteration,numTests,tests,results);

	models = estimationPhase->Evolve();

	if ( simParams->useTheBank ) {
		
		if ( TestTooHard() ) {
			DepositTest();
			HalveTestSearchEffort();
			RestoreOldModels();
		}
		else
			BackupModels();
	}

	if ( currIteration == 1 )
		SaveRealFitness(false);
	else
		SaveRealFitness(true);

	Print();

	delete estimationPhase;
	estimationPhase = NULL;
}

void EE::Explore(int currIteration) {

	if ( simParams->evolvingTests ) {

		explorationPhase = new EXP_PHASE(models,target->targetData);

		tests[numTests++] = explorationPhase->Evolve(numTests,tests,numBankedTests,bankedTests,target->targetData);

		delete explorationPhase;
		explorationPhase = NULL;
	}
	else {
	
		//tests[numTests++] = new TEST(target->targetData);
		//Random testing
		
		int startIndex = currIteration*STARTING_EXPERIMENT_LENGTH;

		while ( startIndex > (target->targetData->length-STARTING_EXPERIMENT_LENGTH) )
			startIndex = startIndex - target->targetData->length;

		while ( startIndex < 0 )
			startIndex = startIndex + STARTING_EXPERIMENT_LENGTH;

		tests[numTests++] = new TEST(target->targetData,startIndex);
		//Uniform testing
	}

	SaveTest(tests[numTests-1],true);
}

void EE::HalveTestSearchEffort(void) {

	simParams->expPhaseMaxGens = int(double(simParams->expPhaseMaxGens) / 2.0);
}

void EE::InitResults(void) {

	results = new RESULT * [simParams->maxCycles];

	for (int i=0;i<simParams->maxCycles;i++)
		results[i] = NULL;

	bankedResults = new RESULT * [simParams->maxCycles];

	for (i=0;i<simParams->maxCycles;i++)
		bankedResults[i] = NULL;
}

void EE::InitTests(void) {

	tests = new TEST * [simParams->maxCycles];
	for (int i=0;i<simParams->maxCycles;i++)
		tests[i] = NULL;

	numTests = 0;

	bankedTests = new TEST * [simParams->maxCycles];
	ages = new int[simParams->maxCycles];

	for (i=0;i<simParams->maxCycles;i++) {
		bankedTests[i] = NULL;
		ages[i] = 0;
	}

	numBankedTests = 0;
}

void EE::LoadTarget(void) {

	ifstream *inFile;
	inFile = new ifstream(simParams->targetFileName);

	target = new MODEL(inFile);

	inFile->close();
	delete inFile;
	inFile = NULL;

	simParams->numEqs = target->targetData->width;

	simParams->tempY = new MATRIX(1,simParams->numEqs);
	simParams->kn1 = new MATRIX(1,simParams->numEqs);
	simParams->kn2 = new MATRIX(1,simParams->numEqs);
	simParams->kn3 = new MATRIX(1,simParams->numEqs);
	simParams->kn4 = new MATRIX(1,simParams->numEqs);

	//simParams->CreateTargetMatrices(target);
}

void EE::Print(void) {

	if ( simParams->writeToScreen ) {

		printf("[Iteration %d of %d]:\t[tests in play: %d; tests banked: %d]\n",
			currCycle,
			simParams->maxCycles-1,
			numTests,
			numBankedTests);

		estimationPhase->Print(5);
	}
}

void EE::RestoreOldModels(void) {

	int k=0;

	for (int i=0;i<NUM_MODELS;i++)

		for (int j=0;j<simParams->numEqs;j++) {

			delete models[i]->diffEqs[j]->root;

			models[i]->diffEqs[j]->root = simParams->previousEqs[k];

			simParams->previousEqs[k] = NULL;

			k++;
		}
}

void EE::SaveRealFitness(int append) {

	TEST *t;
	RESULT *targetResult;
	RESULT *modelResult;
	double diff = 0.0;

	for (int i=0;i<EVALS_FOR_OBJ_FITNESS;i++) {

		t = new TEST(target->targetData);

		targetResult = new RESULT;
		modelResult = new RESULT;

		simParams->evaluatingTarget = true;
		target->Evaluate(t,targetResult,simParams->lengthOfExperiment);
		simParams->evaluatingTarget = false;

		models[0]->Evaluate(t,modelResult,simParams->lengthOfExperiment);

		diff = diff + targetResult->ComputeDifference(modelResult);

		delete modelResult;
		delete targetResult;
		delete t;
	}

	diff = diff / double(EVALS_FOR_OBJ_FITNESS);

	char fileName[50];

	sprintf(fileName,"data/run%d_%d_real.dat",simParams->randSeed,simParams->GetExperimentRegime());

	if ( append )
		realFitFile = new ofstream(fileName,ios::app);
	else
		realFitFile = new ofstream(fileName);

	time_t currTime = time(NULL);
	double timeSinceStart = currTime - simParams->startTime;

	(*realFitFile) << currCycle << " ";
	(*realFitFile) << diff << " ";
	(*realFitFile) << timeSinceStart << " ";
	(*realFitFile) << simParams->targetEvals << " ";
	(*realFitFile) << simParams->modelEvals << "\n";

	realFitFile->close();
	delete realFitFile;
	realFitFile = NULL;

	if ( diff < TERMINATION_THRESHOLD )
		exit(0);
}

void EE::SaveTest(TEST *test, int append) {

	char fileName[50];

	sprintf(fileName,"data/run%d_%d_tests.dat",simParams->randSeed,simParams->GetExperimentRegime());

	if ( append )
		testFile = new ofstream(fileName,ios::app);
	else
		testFile = new ofstream(fileName);

	(*testFile) << currCycle << " ";

	test->Save(testFile);

	testFile->close();
	delete testFile;
	testFile = NULL;
}

int EE::TestTooHard(void) {

	if ( numTests <= 2 )
		return( false );

	int testTooHard = false;

	int i=0;

	while ( (!testTooHard) && (i<simParams->numEqs) ) {

		if ( (simParams->previousFits->Get(i,0)*BANKED_TEST_CUTOFF) < (simParams->storedFits->Get(i,0)) )
			testTooHard = true;

		i++;
	}

	return( testTooHard );
}

int EE::WithdrawTest(void) {

	int testWithdrawn = false;

	/*
	if ( numBankedTests > 0 ) {

		RESULT *tempResult;
		double diff;
		int testFound = false;
		int i = 0;

		while ( (i<(numBankedTests-1)) && (!testFound) ) {

			tempResult = new RESULT;

			models[0]->Evaluate(bankedTests[i],tempResult,bankedResults[i]->GetLengthOfTest());
			simParams->modelEvals++;

			diff = bankedResults[i]->ComputeDifference(tempResult);

			if ( (diff < (2.0*MODEL_FIT_THRESHOLD)) && 
				 (ages[i]>0) )

				testFound = true;

			else
				i++;

			delete tempResult;
			tempResult = NULL;
		}

		if ( testFound ) {

			tests[numTests] = bankedTests[i];
			bankedTests[i] = NULL;

			results[numTests] = bankedResults[i];
			bankedResults[i] = NULL;

			numTests++;

			while ( i<(numBankedTests-1) ) {
				bankedTests[i] = bankedTests[i+1];
				bankedTests[i+1] = NULL;
				ages[i] = ages[i+1];
				bankedResults[i] = bankedResults[i+1];
				bankedResults[i+1] = NULL;
				i++;
			}
			numBankedTests--;
			testWithdrawn = true;
			
			if ( simParams->writeToScreen )
				printf("Test withdrawn.\n");
		}
	}
	*/

	return( testWithdrawn );
}

#endif