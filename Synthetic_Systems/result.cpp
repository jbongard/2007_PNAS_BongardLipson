

#include "stdafx.h"

#ifndef _RESULT_CPP
#define _RESULT_CPP

#include "constants.h"
#include "result.h"
#include "simParams.h"

extern SIM_PARAMS *simParams;

RESULT::RESULT(void) {

	outputData = NULL;
}

RESULT::~RESULT(void) {

	if ( outputData ) {
		delete outputData;
		outputData = NULL;
	}
}

void RESULT::AddNoise(void) {

	outputData->AddNoise(simParams->noiseAmt);
}

double RESULT::ComputeDifference(RESULT *otherResult) {

	double maxDiff = -1;

	double currDiff = 0.0;

	for (int j=0;j<simParams->numEqs;j++) {
		currDiff = outputData->MaxDiffCol(j,otherResult->outputData,j);
		if ( currDiff > maxDiff )
			maxDiff = currDiff;
	}
	return( maxDiff );

	//for (int j=0;j<simParams->numEqs;j++)
	//	currDiff = currDiff + outputData->MaxDiffCol(j,otherResult->outputData,j);
	//return( currDiff / double(simParams->numEqs) );

}

int  RESULT::GetLengthOfTest(void) {

	return( outputData->length );
}

void RESULT::Print(void) {

	outputData->Print(true);
}

void RESULT::UpdateVarRanges(void) {

	if ( simParams->minVarValues == NULL ) {
		simParams->minVarValues = new MATRIX(1,simParams->numEqs,1000000.0);
		simParams->maxVarValues = new MATRIX(1,simParams->numEqs,-1000000.0);
	}

	for (int j=0;j<simParams->numEqs;j++) {

		double currMin = outputData->MinValInColumn(j);
		double currMax = outputData->MaxValInColumn(j);

		if ( currMin < simParams->minVarValues->Get(0,j) )
			simParams->minVarValues->Set(0,j,currMin);

		if ( currMax > simParams->maxVarValues->Get(0,j) )
			simParams->maxVarValues->Set(0,j,currMax);
	}
}

#endif