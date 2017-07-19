/* ---------------------------------------------------
   FILE:     bodyPlan.cpp
	AUTHOR:   Josh Bongard
	DATE:     October 2, 2000
	FUNCTION: This class contains all information for
				 a single physical segment of an organism
				 in the MathEngine environment.
 -------------------------------------------------- */

#include "stdafx.h"
#include "math.h"

#ifndef _MODEL_CPP
#define _MODEL_CPP

#include "constants.h"
#include "model.h"
#include "simParams.h"

extern SIM_PARAMS *simParams;

MODEL::MODEL(MODEL *parent, int copyParams, int copyEqs) {

	fitness = 0.0;

	targetData = NULL;

	if ( copyEqs ) {

		diffEqs = new DIFF_EQ *	[simParams->numEqs];

		for (int i=0;i<simParams->numEqs;i++)
			diffEqs[i] = new DIFF_EQ(parent->diffEqs[i]);
	}
	else
		CreateDiffEqs();
}

MODEL::MODEL(ifstream *inFile) {

	targetData = new MATRIX(inFile);

	//targetData->Mult(1.0/targetData->Min());
	//targetData->Add(-1.0);

	simParams->dataPoints = targetData->length;

	params = NULL;
	diffEqs = NULL;
}

MODEL::~MODEL(void) {

	if ( targetData ) {
		delete targetData;
		targetData = NULL;
	}
	else {

		for (int i=0;i<simParams->numEqs;i++) {

			if ( diffEqs[i] ) {
				delete diffEqs[i];
				diffEqs[i] = NULL;
			}
		}

		delete[] diffEqs;
		diffEqs = NULL;
	}
}

void MODEL::Cross(MODEL *sister) {

	if ( simParams->evolveParams )
		params->CrossRow(0,sister->params,0);

	if ( simParams->evolveEqs && (simParams->Rand(0.0,1.0)<0.5) )
		CrossEqs(sister);
}

void MODEL::Evaluate(TEST *test, RESULT *result, int lengthOfTest) {

	result->outputData = new MATRIX(lengthOfTest,simParams->numEqs);

	if ( simParams->evaluatingTarget )
	
		EvaluateTarget(test,result,lengthOfTest);

	else

		EvaluateModel(test,result,lengthOfTest);
}

void MODEL::Evaluate(TEST *test, RESULT *result, int lengthOfTest, RESULT *otherResult) {

	result->outputData = new MATRIX(lengthOfTest,simParams->numEqs,0.0);

	if ( simParams->currY == NULL )
		simParams->currY = new MATRIX(test->y0);
	else
		simParams->currY->Copy(test->y0);

	if ( simParams->nextY == NULL )
		simParams->nextY = new MATRIX(simParams->currY);

	int k = 0;
	double currT = 0.0;
	int valid = true;

	while (	(valid) && (k < (lengthOfTest-1)) && (!exploding) ) {

		simParams->nextY->CopyRow(0,otherResult->outputData,k+1);

		result->outputData->Set(k,
								simParams->DiffEqUnderScrutiny(),
								simParams->currY->Get(0,simParams->DiffEqUnderScrutiny()));

		EvaluateDiffEqs(simParams->DiffEqUnderScrutiny(),currT,simParams->currY,simParams->nextY,&valid);

		if ( ((fabs(simParams->currY->Get(0,simParams->DiffEqUnderScrutiny())) > pow(10.0,6.0)) &&
			 (k<(lengthOfTest-1))) ||
			 (!valid) ) {

			exploding = true;
		}

		currT = currT + h;
		k++;
	}
}

int MODEL::GetMaxDepth(void) {

	int maxDepth = -1;

	return( diffEqs[simParams->DiffEqUnderScrutiny()]->GetMaxDepth() );
}

void MODEL::Mutate(void) {

	int numMuts = simParams->RandInt(1,4);

	for (int i=0;i<numMuts;i++)
		MutateEqs();
}

void MODEL::Print(void) {

	printf("--------------------------\n");
	printf("Fit: %8.8f\n",fitness);

	diffEqs[simParams->DiffEqUnderScrutiny()]->Print();
	
	printf("\n");
	char ch = getchar();
}

void MODEL::PrintInfix(ofstream *outStream) {

	for (int i=0;i<simParams->numEqs;i++) {

		diffEqs[i]->PrintInfix();
	}
}

void MODEL::Save(ofstream *outFile) {

	if ( simParams->evolveEqs ) {

		for (int i=0;i<simParams->numEqs;i++) {

			(*outFile) << "eval('yprime(" << i+1 << ") = ";

			diffEqs[i]->Save(outFile);

			(*outFile) << ";');\n";
		}

		(*outFile) << "\n";
	}
}

double MODEL::Similarity(MODEL *other) {

	if ( !simParams->evolveEqs )
		return( params->MSQRow(0,other->params,0) );
	else
		return( SimilarityBetEqs(other) );
}

int MODEL::Size(void) {

	int size = 0;

	for (int i=0;i<simParams->numEqs;i++)
		size = size + diffEqs[i]->Count();

	return( size );
}

int  MODEL::WorseThanParent(void) {

	return( fitness > parentFitness );
}

// ----------------- Private Functions ------------------

void MODEL::CreateDiffEqs(void) {

	diffEqs = new DIFF_EQ *	[simParams->numEqs];

	for (int i=0;i<simParams->numEqs;i++)
		diffEqs[i] = new DIFF_EQ;
}

void MODEL::CreateRandomParams(void) {

	params = new MATRIX(1,simParams->numParams);

	for (int i=0;i<simParams->numParams;i++)
		params->Set(0,i,simParams->Rand(VAL_MIN,VAL_MAX));
}

void MODEL::CrossEqs(MODEL *sister) {

	int myEq = simParams->RandInt(0,simParams->numEqs-1);
	int herEq = simParams->RandInt(0,simParams->numEqs-1);

	diffEqs[myEq]->Cross( sister->diffEqs[herEq] );
}

void MODEL::EvaluateDiffEqs(double currT, MATRIX *currY, int *valid) {

	int i;

	// ---------------------------- kn1 -----------------------------

	for (i=0;i<simParams->numEqs;i++)

		simParams->kn1->vals[i] = diffEqs[i]->Evaluate(currT,currY,params,&*valid);

	for (i=0;i<simParams->numEqs;i++)

		simParams->tempY->vals[i] = simParams->kn1->vals[i]*(h/2.0) + currY->vals[i];

	// ---------------------------- kn2 -----------------------------

	for (i=0;i<simParams->numEqs;i++)

		simParams->kn2->vals[i] = diffEqs[i]->Evaluate(currT+(h/2.0),simParams->tempY,params,&*valid);

	for (i=0;i<simParams->numEqs;i++)

		simParams->tempY->vals[i] = simParams->kn2->vals[i]*(h/2.0) + currY->vals[i];

	// ---------------------------- kn3 -----------------------------

	for (i=0;i<simParams->numEqs;i++)

		simParams->kn3->vals[i] = diffEqs[i]->Evaluate(currT+(h/2.0),simParams->tempY,params,&*valid);

	for (i=0;i<simParams->numEqs;i++)

		simParams->tempY->vals[i] = simParams->kn3->vals[i]*h + currY->vals[i];

	// ---------------------------- kn4 -----------------------------

	for (i=0;i<simParams->numEqs;i++)

		simParams->kn4->vals[i] = diffEqs[i]->Evaluate(currT+h,simParams->tempY,params,&*valid);

	// ---------------------------- Combining -----------------------------

	for (i=0;i<simParams->numEqs;i++)

		currY->vals[i] = currY->vals[i] + (h/6.0) * ( (simParams->kn1->vals[i]) + 
													  (2.0*simParams->kn2->vals[i]) + 
													  (2.0*simParams->kn3->vals[i]) + 
													  (simParams->kn4->vals[i]) );
}

void MODEL::EvaluateDiffEqs(int i, double currT, MATRIX *currY, MATRIX *nextY, int *valid) {

	int j;

	// ---------------------------- kn1 -----------------------------
	for (j=0;j<=simParams->numEqs-1;j++)

		if ( j==i )
			simParams->kn1->vals[j] = diffEqs[j]->Evaluate(currT,currY,params,&*valid);
		else
			simParams->kn1->vals[j] = (nextY->vals[j] - currY->vals[j])/h;

	for (j=0;j<=simParams->numEqs-1;j++)

		simParams->tempY->vals[j] = simParams->kn1->vals[j]*(h/2.0) + currY->vals[j];

	// ---------------------------- kn2 -----------------------------
	for (j=0;j<=simParams->numEqs-1;j++)

		if ( j==i )
			simParams->kn2->vals[j] = diffEqs[j]->Evaluate(currT+(h/2.0),simParams->tempY,params,&*valid);
		else
			simParams->kn2->vals[j] = (nextY->vals[j] - currY->vals[j])/h;

	for (j=0;j<=simParams->numEqs-1;j++)

		simParams->tempY->vals[j] = simParams->kn2->vals[j]*(h/2.0) + currY->vals[j];

	// ---------------------------- kn3 -----------------------------
	for (j=0;j<=simParams->numEqs-1;j++)

		if ( j==i )
			simParams->kn3->vals[j] = diffEqs[j]->Evaluate(currT+(h/2.0),simParams->tempY,params,&*valid);
		else
			simParams->kn3->vals[j] = (nextY->vals[j] - currY->vals[j])/h;

	for (j=0;j<=simParams->numEqs-1;j++)

		simParams->tempY->vals[j] = simParams->kn3->vals[j]*h + currY->vals[j];

	// ---------------------------- kn4 -----------------------------
	for (j=0;j<=simParams->numEqs-1;j++)

		if ( j==i )
			simParams->kn4->vals[j] = diffEqs[j]->Evaluate(currT+h,simParams->tempY,params,&*valid);
		else
			simParams->kn4->vals[j] = (nextY->vals[j] - currY->vals[j])/h;
	
	// ---------------------------- Combining -----------------------------
	for (j=0;j<=simParams->numEqs-1;j++)

		if ( j==i )
			currY->vals[j] = currY->vals[j] + (h/6.0) * ( (simParams->kn1->vals[j]) + 
														  (2.0*simParams->kn2->vals[j]) + 
														  (2.0*simParams->kn3->vals[j]) + 
														  (simParams->kn4->vals[j]) );
		else
			currY->vals[j] = nextY->vals[j];
}

void MODEL::EvaluateModel(TEST *test, RESULT *result, int lengthOfTest) {

	if ( simParams->currY == NULL )
		simParams->currY = new MATRIX(test->y0);
	else
		simParams->currY->Copy(test->y0);

	int k = 0;
	double currT = 0.0;
	int valid = true;

	while ( (valid) && (k < lengthOfTest) && (!simParams->currY->Exploding()) ) {

		result->outputData->CopyRow(k,simParams->currY,0);

		EvaluateDiffEqs(currT,simParams->currY,&valid);

		currT = currT + h;
		k++;
	}
}

void MODEL::EvaluateTarget(TEST *test, RESULT *result, int lengthOfTest) {

	for (int i=test->chosenIndexIntoData;
		 i<(test->chosenIndexIntoData+lengthOfTest);
		 i++) {

		for (int j=0;j<simParams->numEqs;j++)

			result->outputData->Set(i-test->chosenIndexIntoData,j,targetData->Get(i,j));
	}
}

void MODEL::MutateEqs(void) {

	//int whichEq = simParams->RandInt(0,simParams->numEqs-1);
	int whichEq = simParams->DiffEqUnderScrutiny();

	diffEqs[whichEq]->Mutate();
	//diffEqs[whichEq]->Clean();
}

double MODEL::SimilarityBetEqs(MODEL *other) {

	double sim = 0.0;

	for (int i=0;i<simParams->numEqs;i++)
		sim = sim + diffEqs[i]->Similarity(other->diffEqs[i]);

	return( sim / double(simParams->numEqs) );
}

#endif