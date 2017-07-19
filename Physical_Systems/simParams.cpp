/* ---------------------------------------------------
   FILE:     simParams.cpp
	AUTHOR:   Josh Bongard
	DATE:     October 5, 2000
	FUNCTION: This class contains all miscellaneous
				 data and functions for this simulation.
 -------------------------------------------------- */

#include "stdafx.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#ifndef _SIM_PARAMS_CPP
#define _SIM_PARAMS_CPP

#include "constants.h"
#include "simParams.h"

SIM_PARAMS::SIM_PARAMS(int argc, char **argv) {

	randSeed = RANDOM_SEED;
	evolvingTests = false;
	evolveEqs = false;
	evolveParams = false;
	lengthOfExperiment = STARTING_EXPERIMENT_LENGTH;
	writeToScreen = true;
	useTheBank = false;
	softGP = false;
	maxCycles = MAX_CYCLES;
	modelEvals = 0;
	targetEvals = 0;
	targetDepth = MAX_MODEL_DEPTH;
	estPhaseMaxGens = MAX_EST_PHASE_GENS;
	minVarValues = NULL;
	maxVarValues = NULL;
	useTimer = false;
	maxMinutes = 0;

	storedEqs = NULL;
	previousEqs = NULL;
	storedFits = NULL;
	storedDepths = NULL;
	storedSizes = NULL;

	currY = NULL;
	nextY = NULL;
	testMin = TEST_MIN;
	testMax = TEST_MAX;
	prune = false;
	dumbSnip = false;
	snip = false;
	creep = false;
	worstTest = false;
	expPhaseMaxGens = GENS_PER_EPOCH;
	testDist = false;
	mutateCrunch = false;
	mutateSplit = false;
	mutateMerge = false;
	mutateNeutrally = false;
	fullDis = false;
	resultDist = false;
	addNoise = false;
	noiseAmt = 0.0;
	fullTestDist = false;

	startTime = time(NULL);

	sprintf(targetFileName,"%s%s",TARGET_FILENAME,"test.dat");

	ParseParameters(argc,argv);

	srand(randSeed);
}

SIM_PARAMS::~SIM_PARAMS(void) {

	WriteEndFile();

	delete tempY;
	tempY = NULL;

	delete kn1;
	kn1 = NULL;

	delete kn2;
	kn2 = NULL;

	delete kn3;
	kn3 = NULL;

	delete kn4;
	kn4 = NULL;

	delete currY;
	currY = NULL;

	delete nextY;
	nextY = NULL;

	delete[] storedEqs;
	storedEqs = NULL;

	if ( previousEqs ) {

		for (int i=0;i<(numEqs*NUM_MODELS);i++) {

			if ( previousEqs[i] ) {
				delete previousEqs[i];
				previousEqs[i] = NULL;
			}
		}

		delete[] previousEqs;
		previousEqs = NULL;
	}

	delete storedDepths;
	storedDepths = NULL;

	delete storedFits;
	storedFits = NULL;

	delete storedSizes;
	storedSizes = NULL;

	delete previousFits;
	previousFits = NULL;

	delete minVarValues;
	minVarValues = NULL;

	delete maxVarValues;
	maxVarValues = NULL;

	/*
	for (int i=0;i<numEqs;i++) {
		delete targetMatrices[i];
		targetMatrices[i] = NULL;
	}

	delete[] targetMatrices;
	targetMatrices = NULL;
	*/
}

void  SIM_PARAMS::CreateTargetMatrices(MODEL *target) {

	targetMatrices = new MATRIX * [numEqs];

	for (int i=0;i<numEqs;i++)
		targetMatrices[i] = target->diffEqs[i]->root->AllExtractTermsAndDegrees();
}

int   SIM_PARAMS::DiffEqUnderScrutiny(void) {

	int genToChangeFocus = estPhaseMaxGens;

	if ( currGen < genToChangeFocus )
		return( 0 );
	else if ( currGen < genToChangeFocus*2 )
		return( 1 );
	else if ( currGen < genToChangeFocus*3 )
		return( 2 );
	else if ( currGen < genToChangeFocus*4 )
		return( 3 );
	else if ( currGen < genToChangeFocus*5 )
		return( 4 );
	else if ( currGen < genToChangeFocus*6 )
		return( 5 );
	else if ( currGen < genToChangeFocus*7 )
		return( 6 );
	else if ( currGen < genToChangeFocus*8 )
		return( 7 );
	else if ( currGen < genToChangeFocus*9 )
		return( 8 );
	else
		return( 9 );
}

void  SIM_PARAMS::FileDelete(char *fileName) {

	char command[100];

	sprintf(command,"del %s",fileName);

	system(command);
}

int   SIM_PARAMS::FileExists(char *fileName) {

	int exists;

	ifstream *inFile = new ifstream(fileName,ios::nocreate);

	if ( inFile->good() ) {
		exists = true;
		inFile->close();
		delete inFile;
		inFile = NULL;
	}
	else {
		exists = false;
	}

	return( exists );
}

void SIM_PARAMS::FileRename(char *src, char *dest) {

	char command[200];

	sprintf(command,"rename %s %s",src,dest);

	system(command);
}

int SIM_PARAMS::FindStrEnd(char *s, int pos, int count) {

	if ( count == 0 )
		return( pos );
	else
		if ( s[pos] == '(' )
			return( FindStrEnd(s,pos+1,count+1) );
		
		else if ( s[pos] == ')' )
			return( FindStrEnd(s,pos+1,count-1) );
		else
			return( FindStrEnd(s,pos+1,count) );
}

int SIM_PARAMS::FlipCoin(void) {

	return( Rand(0.0,1.0) < 0.5 );
}

int SIM_PARAMS::GenOfChange(void) {

	int gensPerEpoch = estPhaseMaxGens;

	int genOfChange = (currGen % gensPerEpoch) == 0;

	return( genOfChange );
}

int SIM_PARAMS::GenOfChangeNext(void) {

	int gensPerEpoch = estPhaseMaxGens;

	int genOfChange = (currGen>0) && (((currGen+1) % gensPerEpoch) == 0);

	return( genOfChange );
}

int SIM_PARAMS::GetExperimentRegime(void) {

	if ( snip )
		return( 2 );
	else if ( evolvingTests )
		return( 1 );
	else
		return( 0 );
	/*
	if ( fullTestDist )
		return( 12 );

	else if ( addNoise )
		return( 11 );

	else if ( mutateNeutrally )
		return( 10 );

	else if ( testDist )
		return( 9 );

	else if ( worstTest )
		return( 8 );

	else if ( evolvingTests )
		return( 7 );
	
	else if ( mutateCrunch )
		return( 6 );

	else if ( creep )
		return( 5 );

	else if ( snip )
		return( 4 );

	else if ( dumbSnip )
		return( 3 );
	
	else if ( prune )
		return( 2 );

	else
		return( 1 );
	*/
}

ofstream *SIM_PARAMS::GetOutFile(char *fileName) {

	ofstream *outFile = new ofstream(fileName);

	return( outFile );
}

void  SIM_PARAMS::ParseParameters(int argc, char **argv) {

	int currParam;

	for(currParam=0;currParam<argc;currParam++) {

		if ( strcmp(argv[currParam],"-r") == 0 )
			randSeed = atoi(argv[currParam+1]);

		if ( strcmp(argv[currParam],"-g") == 0 )
			estPhaseMaxGens = atoi(argv[currParam+1]);

		if ( strcmp(argv[currParam],"-c") == 0 )
			maxCycles = atoi(argv[currParam+1]) + 1;

		if ( strcmp(argv[currParam],"-timer") == 0 ) {
			useTimer = true;
			maxMinutes = atoi(argv[currParam+1]);
		}

		if ( strcmp(argv[currParam],"-ei") == 0 )
			evolvingTests = true;

		if ( strcmp(argv[currParam],"-ep") == 0 )
			evolveParams = true;

		if ( strcmp(argv[currParam],"-null") == 0 )
			writeToScreen = false;

		if ( strcmp(argv[currParam],"-b") == 0 )
			useTheBank = true;

		if ( strcmp(argv[currParam],"-s") == 0 )
			softGP = true;

		if ( strcmp(argv[currParam],"-dumbSnip") == 0 )
			dumbSnip = true;

		if ( strcmp(argv[currParam],"-prune") == 0 )
			prune = true;

		if ( strcmp(argv[currParam],"-snip") == 0 )
			snip = true;

		if ( strcmp(argv[currParam],"-split") == 0 )
			mutateSplit = true;

		if ( strcmp(argv[currParam],"-merge") == 0 )
			mutateMerge = true;

		if ( strcmp(argv[currParam],"-creep") == 0 )
			creep = true;

		if ( strcmp(argv[currParam],"-worstTest") == 0 )
			worstTest = true;

		if ( strcmp(argv[currParam],"-testDist") == 0 )
			testDist = true;

		if ( strcmp(argv[currParam],"-fullTestDist") == 0 )
			fullTestDist = true;

		if ( strcmp(argv[currParam],"-crunch") == 0 )
			mutateCrunch = true;

		if ( strcmp(argv[currParam],"-fullDis") == 0 )
			fullDis = true;

		if ( strcmp(argv[currParam],"-neut") == 0 )
			mutateNeutrally = true;

		if ( strcmp(argv[currParam],"-resultDist") == 0 )
			resultDist = true;

		if ( strcmp(argv[currParam],"-n") == 0 ) {
			addNoise = true;
			noiseAmt = atof(argv[currParam+1]);
		}

		if ( strcmp(argv[currParam],"-tmin") == 0 )
			testMin = atof(argv[currParam+1]);

		if ( strcmp(argv[currParam],"-tmax") == 0 )
			testMax = atof(argv[currParam+1]);

		if ( strcmp(argv[currParam],"-t") == 0 )
			sprintf(targetFileName,"%s%s",TARGET_FILENAME,argv[currParam+1]);

		if ( strcmp(argv[currParam],"-ee") == 0 )
			evolveEqs = true;

	}
}

double SIM_PARAMS::Rand(double min, double max) {

	double zeroToOne = ((double)rand()) / RAND_MAX;
	double returnVal;

	returnVal = (zeroToOne * (max-min)) + min;
	return returnVal;
}

int SIM_PARAMS::RandInt(int min, int max) {

	if ( min == max )
		return( min );
	else {

		int val = (rand() % (max-min+1)) + min;

		if ( val > max )
			val = max;

		return( val );
	}
}

char *SIM_PARAMS::StrCpy(char *s, int start, int end) {

	char *cpy = new char[end-start+2];

	int k=0;
	for (int i=start;i<=end;i++)
		cpy[k++] = s[i];

	cpy[k] = 0;

	return( cpy );
}

void SIM_PARAMS::TuneNoiseByData(void) {

	double maxRange = -1.0;

	for (int j=0;j<numEqs;j++)

		if ( (maxVarValues->Get(0,j) - minVarValues->Get(0,j)) > maxRange )
			maxRange = maxVarValues->Get(0,j) - minVarValues->Get(0,j);


	noiseAmt = (noiseAmt * maxRange)/2.0;
}

// ----------------- Private Functions ------------------

void SIM_PARAMS::WriteEndFile(void) {

	ofstream *outFile;
	char fileName[50];

	sprintf(fileName,"./run%d_%d_end.dat",randSeed,GetExperimentRegime());

	outFile = new ofstream(fileName);

	outFile->close();
	delete outFile;
	outFile = NULL;
}

#endif