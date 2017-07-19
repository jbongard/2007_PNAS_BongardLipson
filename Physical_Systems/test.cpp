
#include "stdafx.h"

#ifndef _TEST_CPP
#define _TEST_CPP

#include "constants.h"
#include "simParams.h"
#include "test.h"

extern SIM_PARAMS *simParams;

TEST::TEST(MATRIX *targetData) {

	chosenIndexIntoData = simParams->RandInt(0,(simParams->dataPoints-1-simParams->lengthOfExperiment));

	y0 = new MATRIX(1,simParams->numEqs,0.0);

	for (int i=0;i<simParams->numEqs;i++)

		y0->Set(0,i,targetData->Get(chosenIndexIntoData,i));

	fitness = 0.0;
}

TEST::TEST(MATRIX *targetData, int startingIndex) {

	chosenIndexIntoData = startingIndex;

	y0 = new MATRIX(1,simParams->numEqs,0.0);

	for (int i=0;i<simParams->numEqs;i++)

		y0->Set(0,i,targetData->Get(chosenIndexIntoData,i));

	fitness = 0.0;
}

TEST::TEST(TEST *otherTest) {

	chosenIndexIntoData = otherTest->chosenIndexIntoData;

	y0 = new MATRIX(otherTest->y0);

	fitness = otherTest->fitness;
}

TEST::~TEST(void) {

	delete y0;
	y0 = NULL;
}

double TEST::ClosestTestDistance(int numOtherTests, TEST **otherTests) {

	double closestDistance = 1000.0;
	double currentDistance;

	for (int i=0;i<numOtherTests;i++) {

		currentDistance = Similarity(otherTests[i]);

		if ( currentDistance < closestDistance )
			closestDistance = currentDistance;
	}

	return( closestDistance );
}

void TEST::Cross(TEST *otherTest) {

	y0->CrossRow(0,otherTest->y0,0);

}

void TEST::Mutate(MATRIX *targetData) {

	//y0->Mutate(simParams->testMin,simParams->testMax);

	if ( simParams->FlipCoin() ) {

		if ( simParams->FlipCoin() )
			chosenIndexIntoData++;
		else
			chosenIndexIntoData--;

		if ( chosenIndexIntoData<0 )
			chosenIndexIntoData = simParams->dataPoints-1-simParams->lengthOfExperiment;
		
		else if ( chosenIndexIntoData > (simParams->dataPoints-1-simParams->lengthOfExperiment) )
			chosenIndexIntoData = 0;
	}
	else {
	
		chosenIndexIntoData = simParams->RandInt(0,simParams->dataPoints-1-simParams->lengthOfExperiment);
	}

	for (int i=0;i<simParams->numEqs;i++)

		y0->Set(0,i,targetData->Get(chosenIndexIntoData,i));
}

void TEST::Print(int wait) {

	printf("[%5.5f]:\t",fitness);
	y0->Print(wait);
}

void TEST::Save(ofstream *outFile) {

	//(*outFile) << 1845 + chosenIndexIntoData << " ";
	(*outFile) << 0 + chosenIndexIntoData << " ";

	y0->Save(outFile);
}

double TEST::Similarity(TEST *other) {

	return( y0->MSQRow(0,other->y0,0) );
}

int  TEST::TestExists(int numT, TEST **t) {

	int found = false;
	int currT = 0;

	while ( (currT<numT) && (!found) ) {

		if ( Equal(t[currT]) )
			found = true;
		else
			currT++;
	}

	return( found );
}

void TEST::ZeroIt(void) {

	y0->ZeroIt();
}

// ----------------- Private Functions ------------------

int  TEST::Equal(TEST *other) {
	
	int equal = y0->Equal(other->y0);

	return( equal );
}

#endif