/* ---------------------------------------------------
   FILE:     bodyPlan.cpp
	AUTHOR:   Josh Bongard
	DATE:     October 2, 2000
	FUNCTION: This class contains all information for
				 a single physical segment of an organism
				 in the MathEngine environment.
 -------------------------------------------------- */

#include "stdafx.h"

#ifndef _EST_PHASE_CPP
#define _EST_PHASE_CPP

#include "stdlib.h"
#include "math.h"

#include "constants.h"
#include "estPhase.h"
#include "simParams.h"

extern SIM_PARAMS *simParams;

EST_PHASE::EST_PHASE(MODEL *target, int c, int testsSoFar, TEST **t, RESULT **r) {

	numOfTests = testsSoFar;
	tests = t;
	results = r;
	children = NULL;
	simParams->currGen = 0;
	numModels = NUM_MODELS;
	currCycle = c;

	InitModels(target);
}

EST_PHASE::EST_PHASE(MODEL **m, int c, int testsSoFar, TEST **t, RESULT **r) {

	numOfTests = testsSoFar;
	tests = t;
	results = r;
	children = NULL;
	simParams->currGen = 0;
	numModels = NUM_MODELS;
	currCycle = c;

	models = m;
}

EST_PHASE::~EST_PHASE(void) {

	tests = NULL;
	results = NULL;
}

MODEL **EST_PHASE::Evolve(void) {

	if ( simParams->partitioning )

		return( Evolve_Partitioning() );

	else

		return( Evolve_NoPartitioning() );

}

double EST_PHASE::GetBestFitness(void) {

	double bestFitness = 10^8;

	for (int i=0;i<numModels;i++)
		if ( models[i]->fitness < bestFitness )
			bestFitness = models[i]->fitness;

	return( bestFitness );
}

void EST_PHASE::Print(int toModel) {

	if ( simParams->writeToScreen ) {

		printf("[T: %d] [I: %d] [G: %d] ",numOfTests,currCycle,simParams->currGen);

		for (int i=0;i<toModel;i++)
			if ( i<(toModel-1) )
				printf("%3.3f ",models[i]->fitness);
			else
				printf("%3.3f\n",models[i]->fitness);
	}
}

// ----------------- Private Functions ------------------

void EST_PHASE::CloseDataFiles(void) {

	fitFile->close();
	delete fitFile;
	fitFile = NULL;

	depthFile->close();
	delete depthFile;
	depthFile = NULL;

	sizeFile->close();
	delete sizeFile;
	sizeFile = NULL;

	modelFile->close();
	delete modelFile;
	modelFile = NULL;
}

void EST_PHASE::CreateChildPopulation(void) {

	children = new MODEL * [numModels];

	for (int i=0;i<numModels;i=i+2) {

		children[i] = new MODEL(models[i],true,true);
		children[i+1] = new MODEL(models[i+1],true,true);

		children[i]->Cross(children[i+1]);

		children[i]->Mutate();
		children[i+1]->Mutate();

		while ( (!simParams->softGP) &&
				((children[i]->GetMaxDepth()>(simParams->targetDepth)) || 
				 (children[i+1]->GetMaxDepth()>(simParams->targetDepth))) ) {

			delete children[i];
			delete children[i+1];

			children[i] = new MODEL(models[i],true,true);
			children[i+1] = new MODEL(models[i+1],true,true);

			children[i]->Cross(children[i+1]);

			children[i]->Mutate();
			children[i+1]->Mutate();
		}
	}
}

void EST_PHASE::CreateNewClimbers(void) {

	children = new MODEL * [numModels];

	for (int i=0;i<numModels;i++) {

		children[i] = new MODEL(models[i],true,true);

		if ( simParams->worstTest )
			children[i]->parentFitness = models[i]->fitness;
		else
			children[i]->parentFitness = models[i]->fitness * double(numOfTests);

		children[i]->Mutate();
	
		while ( (!simParams->softGP) &&
				(children[i]->GetMaxDepth()>(simParams->targetDepth)) ) {

			delete children[i];

			children[i] = new MODEL(models[i],true,true);

			if ( simParams->worstTest )
				children[i]->parentFitness = models[i]->fitness;
			else
				children[i]->parentFitness = models[i]->fitness * double(numOfTests);

			children[i]->Mutate();
		}
	}
}

int  EST_PHASE::Done(void) {

	return( simParams->currGen == (simParams->estPhaseMaxGens*simParams->numEqs) );

	//while ( (models[0]->fitness > MODEL_FIT_THRESHOLD) &&
	//	    (stagnationCounter<STAGNATION_LIMIT) ) {

}

void EST_PHASE::Evaluate(int evaluateParents) {

	int i,j;
	double diff;

	MODEL **modelsToEval;

	if ( evaluateParents )
		modelsToEval = models;
	else
		modelsToEval = children;

	RESULT *currResult;

	for (i=0;i<numModels;i++) {

		modelsToEval[i]->fitness = 0.0;
		modelsToEval[i]->exploding = false;

		j = numOfTests-1;

		while ( (j>=0) &&
				(!modelsToEval[i]->WorseThanParent()) &&
				(!modelsToEval[i]->exploding) ) {

			currResult = new RESULT;

			modelsToEval[i]->Evaluate(tests[j],currResult,results[j]->GetLengthOfTest(),results[j]);
			simParams->modelEvals++;

			diff = currResult->outputData->MSQColumn(	simParams->DiffEqUnderScrutiny(),
														results[j]->outputData,
														simParams->DiffEqUnderScrutiny());

			if ( simParams->worstTest ) {

				if ( diff > modelsToEval[i]->fitness )
					modelsToEval[i]->fitness = diff;
			}
			else
				if ( simParams->prune )
					modelsToEval[i]->fitness = modelsToEval[i]->fitness + 0.5*diff + 
						0.5*modelsToEval[i]->diffEqs[simParams->DiffEqUnderScrutiny()]->root->Count(0);
				else
					modelsToEval[i]->fitness = modelsToEval[i]->fitness + diff;

			delete currResult;
			currResult = NULL;

			j--;
		}

		if ( !simParams->worstTest )
			modelsToEval[i]->fitness = modelsToEval[i]->fitness / double(numOfTests);

		if ( (modelsToEval[i]->fitness > pow(10.0,10.0)) ||
			 (modelsToEval[i]->fitness < -pow(10.0,10.0)) ||
			 (modelsToEval[i]->exploding) )
			modelsToEval[i]->fitness = pow(10.0,10.0);
	}

	modelsToEval = NULL;
}

void EST_PHASE::Evaluate_NoPartitioning(int evaluateParents) {

	int i,j;
	double diff;

	MODEL **modelsToEval;

	if ( evaluateParents )
		modelsToEval = models;
	else
		modelsToEval = children;

	RESULT *currResult;

	for (i=0;i<numModels;i++) {

		modelsToEval[i]->fitness = 0.0;
		modelsToEval[i]->exploding = false;

		j = numOfTests-1;

		while ( (j>=0) &&
				(!modelsToEval[i]->WorseThanParent()) &&
				(!modelsToEval[i]->exploding) ) {

			currResult = new RESULT;

			modelsToEval[i]->Evaluate(tests[j],currResult,results[j]->GetLengthOfTest());
			simParams->modelEvals++;

			diff = 0.0;
			for (int currEq=0;currEq<simParams->numEqs;currEq++)
				diff = diff + currResult->outputData->MSQColumn(currEq,
															results[j]->outputData,
															currEq);
			diff = diff / double(simParams->numEqs);

			if ( simParams->worstTest ) {

				if ( diff > modelsToEval[i]->fitness )
					modelsToEval[i]->fitness = diff;
			}
			else
				if ( simParams->prune )
					modelsToEval[i]->fitness = modelsToEval[i]->fitness + 0.5*diff + 
						0.5*modelsToEval[i]->diffEqs[simParams->DiffEqUnderScrutiny()]->root->Count(0);
				else
					modelsToEval[i]->fitness = modelsToEval[i]->fitness + diff;

			delete currResult;
			currResult = NULL;

			j--;
		}

		if ( !simParams->worstTest )
			modelsToEval[i]->fitness = modelsToEval[i]->fitness / double(numOfTests);

		if ( (modelsToEval[i]->fitness > pow(10.0,10.0)) ||
			 (modelsToEval[i]->fitness < -pow(10.0,10.0)) ||
			 (modelsToEval[i]->exploding) )
			modelsToEval[i]->fitness = pow(10.0,10.0);
	}

	modelsToEval = NULL;
}

MODEL **EST_PHASE::Evolve_NoPartitioning(void) {

	simParams->currGen = 0;

	StoreAllEqs();
	RecoverAllEqs();

	Evaluate_NoPartitioning(true);

	while ( !Done() ) {

		CreateNewClimbers();
		Evaluate_NoPartitioning(false);

		ReplaceClimbers();
		SortModels();

		simParams->currGen++;
	}

	SaveData();

	return(models);
}

MODEL **EST_PHASE::Evolve_Partitioning(void) {

	simParams->currGen = 0;

	StoreAllEqs();

	while ( !Done() ) {

		if ( simParams->GenOfChange() ) {
			RecoverEqs();
			Evaluate(true);
		}

		CreateNewClimbers();
		Evaluate(false);

		ReplaceClimbers();
		SortModels();

		if ( simParams->writeToScreen )
			Print(NUM_MODELS_FOR_TEST);

		if ( simParams->GenOfChangeNext() ) {
			StoreEqs();
		}

		simParams->currGen++;
	}

	RestoreEqs();

	SaveData();

	return( models );
}

int EST_PHASE::FitnessExists(double fit) {

	int found = false;
	int i=0;

	while ( (!found) && (i<numModels) ) {

		if ( models[i]->fitness == fit )
			found = true;
		else
			i++;
	}

	return( found );
}

void EST_PHASE::InitModels(MODEL *target) {

	models = new MODEL * [numModels];

	for (int i=0;i<numModels;i++)

		models[i] = new MODEL(target,!simParams->evolveParams,!simParams->evolveEqs);
}

void EST_PHASE::OpenDataFiles(int append) {

	OpenDepthDataFile(append);
	OpenModelDataFile(append);
	OpenFitnessDataFile(append);
	OpenSizeDataFile(append);
}

void EST_PHASE::OpenFitnessDataFile(int append) {

	char fileName[50];

	sprintf(fileName,"data/run%d_%d_fits.dat",simParams->randSeed,simParams->GetExperimentRegime());

	if ( append )
		fitFile = new ofstream(fileName,ios::app);
	else
		fitFile = new ofstream(fileName);
}

void EST_PHASE::OpenModelDataFile(int append) {

	char fileName[50];

	sprintf(fileName,"data/run%d_%d_models.dat",simParams->randSeed,simParams->GetExperimentRegime());

	if ( append )
		modelFile = new ofstream(fileName,ios::app);
	else
		modelFile = new ofstream(fileName);
}

void EST_PHASE::OpenDepthDataFile(int append) {

	char fileName[50];

	sprintf(fileName,"data/run%d_%d_depths.dat",simParams->randSeed,simParams->GetExperimentRegime());

	if ( append )
		depthFile = new ofstream(fileName,ios::app);
	else
		depthFile = new ofstream(fileName);
}

void EST_PHASE::OpenSizeDataFile(int append) {

	char fileName[50];

	sprintf(fileName,"data/run%d_%d_sizes.dat",simParams->randSeed,simParams->GetExperimentRegime());

	if ( append )
		sizeFile = new ofstream(fileName,ios::app);
	else
		sizeFile = new ofstream(fileName);
}

int  EST_PHASE::ParentChildSimilarity(int index) {

	double p1c1 = models[index]->Similarity(children[index]);
	double p2c2 = models[index+1]->Similarity(children[index+1]);
	double p1c2 = models[index]->Similarity(children[index+1]);
	double p2c1 = models[index+1]->Similarity(children[index]);

	return( simParams->Rand(0.0,1.0) < 0.5 );
	//return( ( p1c1 + p2c2 ) < ( p1c2 + p2c1 ) );
}

void EST_PHASE::PrintBestModel(void) {

}

void EST_PHASE::PrintProgress(void) {

}

void EST_PHASE::RecoverAllEqs(void) {

	int totalSlots = numModels*(simParams->numEqs);

	int i;

	int diffEqUnderScrutiny = simParams->DiffEqUnderScrutiny();
	int startSlot = diffEqUnderScrutiny;
	int k=0;

	for (i=startSlot;i<totalSlots;i=i+(simParams->numEqs)) {

		int l=0;

		for (int j=i;j<i+(simParams->numEqs);j++) {

			models[k]->diffEqs[l]->root = simParams->storedEqs[j];

			simParams->storedEqs[j] = NULL;

			l++;
		}

		models[k]->fitness = pow(10.0,10.0);
		models[k]->parentFitness = pow(10.0,10.0);
		k++;
	}
}

void EST_PHASE::RecoverEqs(void) {

	int totalSlots = numModels*(simParams->numEqs);

	int i;

	int diffEqUnderScrutiny = simParams->DiffEqUnderScrutiny();
	int startSlot = diffEqUnderScrutiny;
	int k=0;

	for (i=startSlot;i<totalSlots;i=i+(simParams->numEqs)) {

		models[k]->diffEqs[diffEqUnderScrutiny]->root = simParams->storedEqs[i];

		simParams->storedEqs[i] = NULL;

		models[k]->fitness = pow(10.0,10.0);
		models[k]->parentFitness = pow(10.0,10.0);
		k++;
	}
}

void EST_PHASE::ReplaceClimbers(void) {

	for (int i=0;i<numModels;i++) {

		if ( children[i]->fitness <= models[i]->fitness ) {

			delete models[i];
			models[i] = children[i];
		}
		else
			delete children[i];

		children[i] = NULL;
	}

	delete[] children;
	children = NULL;
}

void EST_PHASE::ReplaceParentsWithChildren(void) {

	int p1c1;
	MODEL *temp;

	for (int i=0;i<numModels;i=i+2) {

		p1c1 = ParentChildSimilarity(i);

		if ( !p1c1 ) {
			temp = children[i];
			children[i] = children[i+1];
			children[i+1] = temp;
			temp = NULL;
		}

		if ( (children[i]->fitness < models[i]->fitness) &&
			 (!FitnessExists(children[i]->fitness)) ) {
			delete models[i];
			models[i] = children[i];
		}
		else
			delete children[i];

		children[i] = NULL;

		if ( (children[i+1]->fitness < models[i+1]->fitness) &&
			 (!FitnessExists(children[i+1]->fitness)) ) {
			
			delete models[i+1];
			models[i+1] = children[i+1];
		}
		else
			delete children[i+1];

		children[i+1] = NULL;
	}

	delete[] children;
	children = NULL;
}

void EST_PHASE::RestoreEqs(void) {

	int totalSlots = numModels*(simParams->numEqs);

	int i, j;

	int k=0;
	int l;

	for (i=0;i<totalSlots;i=i+(simParams->numEqs)) {

		l = 0;

		for (j=i;j<(i+(simParams->numEqs));j++) {

			models[k]->diffEqs[l]->root = simParams->storedEqs[j];

			simParams->storedEqs[j] = NULL;
			l++;
		}
		k++;
	}
}

void EST_PHASE::ShuffleModels(void) {

	MATRIX *perm = new MATRIX(1,numModels,0.0);
	perm->CreatePermutation(0,numModels-1);

	MODEL **temp = new MODEL * [numModels];

	for (int i=0;i<numModels;i++) {
		temp[i] = models[ int(perm->Get(0,i)) ];
		models[ int(perm->Get(0,i)) ] = NULL;
	}

	delete[] models;
	models = temp;
	temp = NULL;

	delete perm;
	perm = NULL;
}

void EST_PHASE::SaveData(void) {

	if ( currCycle == 1 )
		OpenDataFiles(false);
	else
		OpenDataFiles(true);

	SaveFits();
	SaveDepths();
	SaveSizes();

	//SaveHits(); Under construction

	SaveModels();

	CloseDataFiles();
}

void EST_PHASE::SaveFits(void) {

	(*fitFile) << currCycle << " ";

	for (int j=0;j<numModels;j++)
		for (int i=0;i<simParams->numEqs;i++)
			(*fitFile) << simParams->storedFits->Get(i,j) << " ";

	(*fitFile) << "\n";
}

void EST_PHASE::SaveHits(void) {

	for (int j=0;j<numModels;j++)
		for (int i=0;i<simParams->numEqs;i++)
			
			printf("%d ",models[j]->diffEqs[i]->IsAHit(i));

	printf("\n");
}

void EST_PHASE::SaveDepths(void) {

	(*depthFile) << currCycle << " ";

	for (int j=0;j<numModels;j++)
		for (int i=0;i<simParams->numEqs;i++)
			(*depthFile) << simParams->storedDepths->Get(i,j) << " ";

	(*depthFile) << "\n";
}

void EST_PHASE::SaveModels(void) {

	(*modelFile) << currCycle << " ";

	models[0]->Save(modelFile);
}

void EST_PHASE::SaveSizes(void) {

	(*sizeFile) << currCycle << " ";

	for (int j=0;j<numModels;j++)
		for (int i=0;i<simParams->numEqs;i++)
			(*sizeFile) << simParams->storedSizes->Get(i,j) << " ";

	(*sizeFile) << "\n";
}

void EST_PHASE::SortModels(void) {

	int done = false;
	int i;
	MODEL *temp;

	while (!done) {

		done = true;

		for (i=0;i<(numModels-1);i++) {

			if ( models[i+1]->fitness < models[i]->fitness ) {

				temp = models[i];
				models[i] = models[i+1];
				models[i+1] = temp;
				temp = NULL;
				done = false;
			}
		}
	}
}

void EST_PHASE::StoreAllEqs(void) {

	int totalSlots = numModels*(simParams->numEqs);

	int i, j, l;

	if ( simParams->storedEqs == NULL ) {
		simParams->storedEqs = new EQ_NODE * [totalSlots];
		simParams->storedFits = new MATRIX(simParams->numEqs,numModels,1000.0);
		simParams->previousFits = new MATRIX(simParams->numEqs,numModels,1000.0);
		simParams->storedDepths = new MATRIX(simParams->numEqs,numModels);
		simParams->storedSizes  = new MATRIX(simParams->numEqs,numModels);
	}

	int diffEqUnderScrutiny = simParams->DiffEqUnderScrutiny();
	int startSlot = diffEqUnderScrutiny;
	int k=0;

	for (i=0;i<totalSlots;i=i+(simParams->numEqs)) {

		l = 0;

		for (j=i;j<(i+(simParams->numEqs));j++) {

			simParams->storedEqs[j] = models[k]->diffEqs[l]->root;

			models[k]->diffEqs[l]->root = NULL;
			l++;
		}
		k++;
	}
}

void EST_PHASE::StoreEqs(void) {

	int totalSlots = numModels*(simParams->numEqs);

	int i;

	int diffEqUnderScrutiny = simParams->DiffEqUnderScrutiny();
	int startSlot = diffEqUnderScrutiny;
	int k=0;

	for (i=startSlot;i<totalSlots;i=i+(simParams->numEqs)) {

		if ( simParams->snip )
			models[k]->diffEqs[diffEqUnderScrutiny]->Clean();

		simParams->storedEqs[i] = models[k]->diffEqs[diffEqUnderScrutiny]->root;
		models[k]->diffEqs[diffEqUnderScrutiny]->root = NULL;

		simParams->previousFits->Set(diffEqUnderScrutiny,k,simParams->storedFits->Get(diffEqUnderScrutiny,k));
		simParams->storedFits->Set(diffEqUnderScrutiny,k,models[k]->fitness);
		simParams->storedDepths->Set(diffEqUnderScrutiny,k,simParams->storedEqs[i]->GetDepth());
		simParams->storedSizes->Set(diffEqUnderScrutiny,k,simParams->storedEqs[i]->Count(0));
		k++;
	}
}

#endif