/* ---------------------------------------------------
   FILE:     bodyPlan.cpp
	AUTHOR:   Josh Bongard
	DATE:     October 2, 2000
	FUNCTION: This class contains all information for
				 a single physical segment of an organism
				 in the MathEngine environment.
 -------------------------------------------------- */

#include "stdafx.h"

#ifndef _EXP_PHASE_CPP
#define _EXP_PHASE_CPP

#include "constants.h"
#include "expPhase.h"
#include "math.h"
#include "simParams.h"

extern SIM_PARAMS *simParams;

EXP_PHASE::EXP_PHASE(MODEL **m, MATRIX *targetData) {

	numTests = NUM_TESTS;
	numModels = NUM_MODELS;

	models = m;

	resultMatrix = NULL;

	tests = new TEST * [numTests];

	for (int i=0;i<numTests;i++)
		tests[i] = new TEST(targetData);
}

EXP_PHASE::~EXP_PHASE(void) {

	for (int i=0;i<numTests;i++) {
		delete tests[i];
		tests[i] = NULL;
	}

	delete[] tests;
	tests = NULL;

	models = NULL;
}

TEST *EXP_PHASE::Evolve(int numT1, TEST **t1, int numT2, TEST **t2, MATRIX *targetData) {

	if ( simParams->resultDist )
		CreateResultMatrix(numT1,t1);

	Evaluate(false,numT1,t1,numT2,t2);

	simParams->currGen = 0;

	while ( !Done() ) {

		if ( simParams->writeToScreen )
			Print();

		//CreateChildPopulation();
		CreateNewClimbers(targetData);

		Evaluate(true,numT1,t1,numT2,t2);

		//ReplaceParentsWithChildren();
		ReplaceClimbers();

		SortTests();

		simParams->currGen++;
	}

	TEST *testToReturn = new TEST(tests[0]);

	if ( simParams->writeToScreen )
		testToReturn->Print(false);

	if ( simParams->resultDist )
		DestroyResultMatrix(numT1);

	return( testToReturn );
}

// ----------------- Private Functions ------------------

void EXP_PHASE::ComputeFinalValues(RESULT **results, MATRIX *finalValues) {

	int m;

	for (int j=0;j<NUM_MODELS_FOR_TEST;j++) {

		for (int k=0;k<simParams->numEqs;k++) {

			m=0;
			
			while ( (m<results[j]->outputData->length-1) && 
					( fabs(results[j]->outputData->Get(m,k)) <= pow(10.0,6.0) ) ) {

				m++;
			}

			finalValues->Set(j,k, results[j]->outputData->Get(m,k) );
		}
	}
}

void EXP_PHASE::CreateChildPopulation(MATRIX *targetData) {

	children = new TEST * [numTests];

	for (int i=0;i<numTests;i=i+2) {

		children[i] = new TEST(tests[i]);
		children[i]->fitness = 0.0;

		children[i+1] = new TEST(tests[i+1]);
		children[i+1]->fitness = 0.0;

		children[i]->Cross(children[i+1]);

		children[i]->Mutate(targetData);
		children[i+1]->Mutate(targetData);
	}
}

void EXP_PHASE::CreateHorizontalResults(RESULT **results, TEST **testsToEvaluate, int i) {

	for (int j=0;j<NUM_MODELS_FOR_TEST;j++) {
		
		results[j] = new RESULT;
		
		models[j]->Evaluate(testsToEvaluate[i],results[j],simParams->lengthOfExperiment);

		simParams->modelEvals++;
	}
}

void EXP_PHASE::CreateNewClimbers(MATRIX *targetData) {

	children = new TEST * [numTests];

	for (int i=0;i<numTests;i++) {

		children[i] = new TEST(tests[i]);
		children[i]->fitness = 0.0;

		children[i]->Mutate(targetData);
	}
}

void EXP_PHASE::CreateResultMatrix(int numPreviousTests, TEST **previousTests) {

	resultMatrix = new RESULT * [numModels * numPreviousTests];
	
	for (int i=0;i<numModels;i++)

		for (int j=0;j<numPreviousTests;j++) {

			resultMatrix[i*numPreviousTests + j] = new RESULT;
				
			models[i]->Evaluate(previousTests[j],resultMatrix[i*numPreviousTests + j],simParams->lengthOfExperiment);
		}
}

void EXP_PHASE::DestroyResultMatrix(int numPreviousTests) {

	for (int i=0;i<numModels;i++)

		for (int j=0;j<numPreviousTests;j++) {

			delete resultMatrix[i*numPreviousTests + j];
			resultMatrix[i*numPreviousTests + j] = NULL;
		}

	delete[] resultMatrix;
}

void EXP_PHASE::DestroyResults(RESULT **results) {

	for (int j=0;j<numModels;j++) {

		delete results[j];
		results[j] = NULL;
	}
}

int  EXP_PHASE::Done(void) {

	return( simParams->currGen == simParams->expPhaseMaxGens );
}

void EXP_PHASE::Evaluate(int useChildren, int numT1, TEST **t1, int numT2, TEST **t2) {

	RESULT **results       = new RESULT * [NUM_MODELS_FOR_TEST];
	MATRIX *finalValues    = new MATRIX(NUM_MODELS_FOR_TEST,simParams->numEqs,0.0);
	MATRIX *variances      = new MATRIX(1,simParams->numEqs,0.0);

	TEST **testsToEvaluate = WhichTestsToUse(useChildren);

	for (int i=0;i<numTests;i++) {

		if ( testsToEvaluate[i]->TestExists(numT1,t1) ||
			 testsToEvaluate[i]->TestExists(numT2,t2) )
			
			 testsToEvaluate[i]->fitness = pow(10.0,10.0);

		else {

			CreateHorizontalResults(results,testsToEvaluate,i);

			ComputeFinalValues(results,finalValues);

			SeekExplosions(i,testsToEvaluate,finalValues);

			if ( !simParams->resultDist ) {

				if ( !simParams->fullDis )
					SeekFinalDisagreement(i,testsToEvaluate,finalValues,variances,numT1,t1);
				else
					SeekTotalDisagreement(i,testsToEvaluate,results,numT1,t1,numT2,t2);
			}
			else
				SeekPastDisagreement(i,testsToEvaluate,results,numT1);

			DestroyResults(results);
		}
	}

	testsToEvaluate = NULL;

	delete variances;
	variances = NULL;

	delete finalValues;
	finalValues = NULL;

	delete[] results;
	results = NULL;
}

int  EXP_PHASE::ParentChildSimilarity(int index) {

	double p1c1 = tests[index]->Similarity(children[index]);
	double p2c2 = tests[index+1]->Similarity(children[index+1]);
	double p1c2 = tests[index]->Similarity(children[index+1]);
	double p2c1 = tests[index+1]->Similarity(children[index]);

	return( ( p1c1 + p2c2 ) < ( p1c2 + p2c1 ) );
}

void EXP_PHASE::Print(void) {

	printf("[G: %d]\t",simParams->currGen);

	for (int i=0;i<numTests;i++)

		printf("%3.3f ",tests[i]->fitness);

	printf("\n");
}

void EXP_PHASE::ReplaceClimbers(void) {

	for (int i=0;i<numTests;i++) {

		if ( children[i]->fitness < tests[i]->fitness ) {

			delete tests[i];
			tests[i] = children[i];
		}
		else
			delete children[i];
		children[i] = NULL;
	}

	delete[] children;
	children = NULL;
}

void EXP_PHASE::ReplaceParentsWithChildren(void) {

	int p1c1;
	TEST *temp;

	for (int i=0;i<numTests;i=i+2) {

		p1c1 = ParentChildSimilarity(i);

		if ( !p1c1 ) {
			temp = children[i];
			children[i] = children[i+1];
			children[i+1] = temp;
			temp = NULL;
		}

		if ( children[i]->fitness < tests[i]->fitness ) {

			delete tests[i];
			tests[i] = children[i];
		}
		else
			delete children[i];
		children[i] = NULL;

		if ( children[i+1]->fitness < tests[i+1]->fitness ) {

			delete tests[i+1];
			tests[i+1] = children[i+1];
		}
		else
			delete children[i+1];
		children[i+1] = NULL;
	}

	delete[] children;
	children = NULL;
}

void EXP_PHASE::SeekFinalDisagreement(int i, TEST **testsToEvaluate, MATRIX *finalValues, MATRIX *variances, 
								 int numPreviousTests, TEST **previousTests) {

	if ( testsToEvaluate[i]->fitness == simParams->numEqs * NUM_MODELS_FOR_TEST ) {

		finalValues->ComputeVariances(variances);

		for (int k=0;k<simParams->numEqs;k++)

			if ( !simParams->testDist )
				testsToEvaluate[i]->fitness = testsToEvaluate[i]->fitness + 
											  (1.0/variances->Get(0,k));
			else
				testsToEvaluate[i]->fitness = testsToEvaluate[i]->fitness + 
											  (1.0/
											     (variances->Get(0,k)*
											      testsToEvaluate[i]->ClosestTestDistance(numPreviousTests,previousTests))
											  );
	}
}

void EXP_PHASE::SeekPastDisagreement(int i, TEST **testsToEvaluate, RESULT **results, int numPreviousTests) {
	
	MATRIX *variances = new MATRIX(simParams->lengthOfExperiment,simParams->numEqs,0.0);

	double totalDisagreement = 0.0;

	double valMean;
	double valVariance;

	for (int j=0;j<numModels;j++) {

		for (int k=0;k<simParams->numEqs;k++) {

			for (int l=0;l<simParams->lengthOfExperiment;l++) {

				valMean = 0.0;
				for (int m=0;m<numPreviousTests;m++)

					valMean = valMean + resultMatrix[j*numPreviousTests+m]->outputData->Get(l,k);

				valMean = valMean + results[j]->outputData->Get(l,k);
				valMean = valMean / double(numPreviousTests+1);

				valVariance = 0.0;
				for (m=0;m<numPreviousTests;m++)

					valVariance = valVariance + pow(resultMatrix[j*numPreviousTests+m]->outputData->Get(l,k) - valMean,2.0);

				valVariance = valVariance + pow(results[j]->outputData->Get(l,k) - valMean,2.0);
				valVariance = valVariance / double(numPreviousTests+1);
				
				variances->Set(l,k,valVariance);
			}
		}

		totalDisagreement = totalDisagreement + variances->Mean();
	}

	testsToEvaluate[i]->fitness =	testsToEvaluate[i]->fitness + 
									(1.0/totalDisagreement);

	delete variances;
	variances = NULL;
}

void EXP_PHASE::SeekTotalDisagreement(int i, TEST **testsToEvaluate, RESULT **results, 
									  int numPreviousTests, TEST **previousTests,
									  int numBankedTests,   TEST **bankedTests) {
	
	MATRIX *variances = new MATRIX(simParams->lengthOfExperiment,simParams->numEqs);

	double valMean;
	double valVar;

	for (int j=0;j<simParams->lengthOfExperiment;j++) {

		for (int k=0;k<simParams->numEqs;k++) {

			valMean = 0.0;
			for (int l=0;l<numModels;l++) {

				valMean = valMean + results[l]->outputData->Get(j,k);
			}
			valMean = valMean / double(numModels);

			valVar = 0.0;
			for (l=0;l<numModels;l++) {

				valVar = valVar + pow(results[l]->outputData->Get(j,k) - valMean,2.0);
			}

			valVar = valVar / double(numModels);

			variances->Set(j,k,valVar);
		}
	}

	if ( (simParams->fullTestDist) || (simParams->testDist) ) {

		double distFromClosestPreviousTest = testsToEvaluate[i]->ClosestTestDistance(numPreviousTests,previousTests);
		double distFromClosestBankedTest;
		double closestDist;

		if ( numBankedTests > 0 ) {

			distFromClosestBankedTest = testsToEvaluate[i]->ClosestTestDistance(numBankedTests,bankedTests);
			if ( distFromClosestBankedTest < distFromClosestPreviousTest )
				closestDist = distFromClosestBankedTest;
			else
				closestDist = distFromClosestPreviousTest;
		}
		else
			closestDist = distFromClosestPreviousTest;

		if ( simParams->testDist )
			testsToEvaluate[i]->fitness = testsToEvaluate[i]->fitness + 
										  (1.0/variances->Mean()) * 
										  (1.0/closestDist);
		else
			testsToEvaluate[i]->fitness = testsToEvaluate[i]->fitness + 
										  (1.0/closestDist);
	}
	else

		testsToEvaluate[i]->fitness = testsToEvaluate[i]->fitness + 
									  (1.0/variances->Mean());

	delete variances;
	variances = NULL;
}

void EXP_PHASE::SeekExplosions(int i, TEST **testsToEvaluate, MATRIX *finalValues) {

	testsToEvaluate[i]->fitness = simParams->numEqs * NUM_MODELS_FOR_TEST;

	for (int k=0;k<simParams->numEqs;k++) {
		for (int j=0;j<NUM_MODELS_FOR_TEST;j++) {

			if ( (fabs(finalValues->Get(j,k)) >= pow(10.0,6.0)) )

				testsToEvaluate[i]->fitness--;
		}
	}
}

void EXP_PHASE::ShuffleTests(void) {

	MATRIX *perm = new MATRIX(1,numTests,0.0);
	perm->CreatePermutation(0,numTests-1);

	TEST **temp = new TEST * [numTests];

	for (int i=0;i<numTests;i++) {
		temp[i] = tests[ int(perm->Get(0,i)) ];
		tests[ int(perm->Get(0,i)) ] = NULL;
	}

	delete[] tests;
	tests = temp;
	temp = NULL;

	delete perm;
	perm = NULL;
}

void EXP_PHASE::SortTests(void) {

	int done = false;
	int i;
	TEST *temp;

	while (!done) {

		done = true;

		for (i=0;i<(numTests-1);i++) {

			if ( tests[i+1]->fitness < tests[i]->fitness ) {

				temp = tests[i];
				tests[i] = tests[i+1];
				tests[i+1] = temp;
				temp = NULL;
				done = false;
			}
		}
	}
}

TEST **EXP_PHASE::WhichTestsToUse(int useChildren) {

	if ( useChildren )
		return( children );
	else
		return( tests );
}

#endif