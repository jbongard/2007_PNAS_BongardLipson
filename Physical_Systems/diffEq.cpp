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

#ifndef _DIFF_EQ_CPP
#define _DIFF_EQ_CPP

#include "constants.h"
#include "diffEq.h"
#include "simParams.h"

extern SIM_PARAMS *simParams;

DIFF_EQ::DIFF_EQ(void) {

	root = NULL;
	softTree = NULL;

	CreateRandom();
}

DIFF_EQ::DIFF_EQ(DIFF_EQ *parent) {

	root = NULL;
	softTree = NULL;

	if ( simParams->softGP )
		softTree = new MATRIX(parent->softTree);
	else
		if ( parent->root )
			root = new EQ_NODE(parent->root);
}

DIFF_EQ::DIFF_EQ(ifstream *inFile) {

	char stringEq[500];

	(*inFile) >> stringEq;

	int i=0;
	int found = false;

	while (!found)
		if ( stringEq[i]==0 )
			found = true;
		else
			i++;

	root = new EQ_NODE(stringEq,0,i-1);
	softTree = NULL;
}

DIFF_EQ::~DIFF_EQ(void) {

	if ( root ) {
		delete root;
		root = NULL;
	}

	if ( softTree ) {
		delete softTree;
		softTree = NULL;
	}
}

void DIFF_EQ::Clean(void) {

	root->Clean();
}

int  DIFF_EQ::Count(void) {

	return( root->Count(0) );
}

void DIFF_EQ::Cross(DIFF_EQ *sister) {

	if ( root )
		CrossRoot(sister);
	else
		CrossSoftTree(sister);
}

double DIFF_EQ::Evaluate(double currT, MATRIX *y, MATRIX *params, int *valid) {

	if ( root )
		return( root->Evaluate(currT,y,params,&*valid) );
	else
		return( softTree->Evaluate(0,currT,y,params) );
}

int  DIFF_EQ::GetMaxDepth(void) {

	return( root->GetDepth() );
}

int  DIFF_EQ::IsAHit(int currEq) {

	EQ_NODE *temp = new EQ_NODE(root);
	
	temp->AllSplit();

	MATRIX *termDegrees = temp->AllExtractTermsAndDegrees();

	delete temp;
	temp = NULL;

	int allModelTermsInTarget = termDegrees->AllRowsAlsoIn( simParams->targetMatrices[currEq] );

	int allTargetTermsInModel = simParams->targetMatrices[currEq]->AllRowsAlsoIn( termDegrees );

	delete termDegrees;
	termDegrees = NULL;

	return( allModelTermsInTarget && allTargetTermsInModel );
}

void DIFF_EQ::Mutate(void) {

	if ( root )
		MutateRoot();
	else
		MutateSoftTree();
}

void DIFF_EQ::Print(void) {

	if ( simParams->softGP )
		softTree->Print(false);
	else
		root->Print();
	
	printf("\n");
}

void DIFF_EQ::PrintInfix(void) {

	root->PrintInfix();
	printf("\n");
}

void DIFF_EQ::Save(ofstream *outFile) {

	if ( root )
		root->Save(outFile);
	
	else
		softTree->Save(outFile);
}

double DIFF_EQ::Similarity(DIFF_EQ *other) {

	if ( root ) {
		int mySize = root->Count(0);
		int hisSize = other->root->Count(0);

		if ( mySize > hisSize )
			return( fabs(mySize-hisSize) / double(mySize) );
		else
			return( fabs(mySize-hisSize) / double(hisSize) );
	}

	return( softTree->AbsoluteDifference(other->softTree) );
}

// ----------------- Private Functions ------------------

void DIFF_EQ::CreateRandom(void) {

	if ( simParams->softGP ) {
		softTree = new MATRIX(int(pow(simParams->targetDepth,2))-1,SOFT_NODE_WIDTH,0.0);
		softTree->InsertValSomewhereInEachRow(1.0);
	}
	else
		root = new EQ_NODE(0);
}

void DIFF_EQ::CrossRoot(DIFF_EQ *sister) {

	int myTotalNodes = root->Count(0);
	int herTotalNodes = sister->root->Count(0);
		
	int myNodeToCross = simParams->RandInt(0,myTotalNodes-1);
	int herNodeToCross = simParams->RandInt(0,herTotalNodes-1);

	EQ_NODE *myNode = root->Find(&myNodeToCross);
	EQ_NODE *herNode = sister->root->Find(&herNodeToCross);

	EQ_NODE *temp;

	if ( (myNode!=root) && (herNode!=sister->root) ) {

		if ( myNode == myNode->parent->left )
			myNode->parent->left = herNode;
		else
			myNode->parent->right = herNode;

		if ( herNode == herNode->parent->left )
			herNode->parent->left = myNode;
		else
			herNode->parent->right = myNode;

		temp = myNode->parent;

		myNode->parent = herNode->parent;
		herNode->parent = temp;
	}
	else if ( (myNode==root) && (herNode!=sister->root) ) {

		if ( herNode == herNode->parent->left )
			herNode->parent->left = myNode;
		else
			herNode->parent->right = myNode;

		myNode->parent = herNode->parent;

		herNode->parent = NULL;
		root = herNode;
	}
	else if ( (myNode!=root) && (herNode==sister->root) ) {

		if ( myNode == myNode->parent->left )
			myNode->parent->left = herNode;
		else
			myNode->parent->right = herNode;

		herNode->parent = myNode->parent;

		myNode->parent = NULL;
		sister->root = myNode;
	}

	myNode = NULL;
	herNode = NULL;
	temp = NULL;
}

void DIFF_EQ::MutateNeutrally(void) {

	int numNodes = root->Count(0);

	int nodeToMutate = simParams->RandInt(0,numNodes-1);

	EQ_NODE *targettedNode = root->Find(&nodeToMutate);

	targettedNode->MutateNeutrally();
}

void DIFF_EQ::CrossSoftTree(DIFF_EQ *sister) {

	softTree->CrossSubTrees(simParams->RandInt(0,softTree->length-1),sister->softTree);
}

void DIFF_EQ::MutateRoot(void) {

	if ( simParams->mutateNeutrally ) {

		root->Clean();

		if ( simParams->FlipCoin() )
			MutateNeutrally();
	}

	int numNodes = root->Count(0);

	int nodeToMutate = simParams->RandInt(0,numNodes-1);

	EQ_NODE *targettedNode = root->Find(&nodeToMutate);

	targettedNode->Mutate();

	root->ResetNodeRanges();

	targettedNode = NULL;
}

void DIFF_EQ::MutateSoftTree(void) {

	int numMuts = simParams->RandInt(1,10);

	for (int i=0;i<numMuts;i++)
		softTree->ChangePairInRow(simParams->RandInt(0,softTree->length-1));
}

#endif