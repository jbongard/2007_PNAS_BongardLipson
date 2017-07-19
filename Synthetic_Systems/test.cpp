
#include "stdafx.h"

#ifndef _TEST_CPP
#define _TEST_CPP

#include "constants.h"
#include "simParams.h"
#include "test.h"

extern SIM_PARAMS *simParams;

TEST::TEST(void) {

	y0 = new MATRIX(1,simParams->numEqs);

	for (int i=0;i<simParams->numEqs;i++)
		y0->Set(0,i,simParams->Rand(simParams->testMin,simParams->testMax));

	fitness = 0.0;
}

TEST::TEST(TEST *otherTest) {

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

void TEST::Mutate(void) {

	y0->Mutate(simParams->testMin,simParams->testMax);
}

void TEST::Print(int wait) {

	printf("[%5.5f]:\t",fitness);
	y0->Print(wait);
}

void TEST::Save(ofstream *outFile) {

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