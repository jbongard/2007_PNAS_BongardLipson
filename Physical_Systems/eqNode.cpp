
#include "stdafx.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"

#ifndef _EQ_NODE_CPP
#define _EQ_NODE_CPP

#include "constants.h"
#include "eqNode.h"
#include "simParams.h"

extern SIM_PARAMS *simParams;

EQ_NODE::EQ_NODE(int depth) {

	parent = NULL;

	minVal = 1000000.0;
	maxVal = -1000000.0;

	if ( depth < (simParams->targetDepth-1) )
		nodeType = simParams->RandInt(0,NUM_OPS-1);
	else
		nodeType = simParams->RandInt(0,NUM_ZERO_PARITY-1);

	if ( nodeType == VAR )
		index = simParams->RandInt(0,simParams->numEqs-1);

	else if ( nodeType == PARAM )
		index = simParams->Rand(VAL_MIN,VAL_MAX);
		//index = double( simParams->RandInt(int(VAL_MIN),(VAL_MAX)) );

	left = NULL;
	right = NULL;

	if ( nodeType >= PLUS ) {

		left = new EQ_NODE(depth+1);
		left->parent = this;

		right = new EQ_NODE(depth+1);
		right->parent = this;
	}
}

EQ_NODE::EQ_NODE(int nt, double i) {

	nodeType = nt;
	index = i;

	minVal = 1000000.0;
	maxVal = -1000000.0;

	left = NULL;
	right = NULL;
}

EQ_NODE::EQ_NODE(EQ_NODE *parentEq) {

	nodeType = parentEq->nodeType;
	index = parentEq->index;
	
	minVal = parentEq->minVal;
	maxVal = parentEq->maxVal;

	parent = NULL;

	if ( parentEq->left ) {
		left = new EQ_NODE(parentEq->left);
		left->parent = this;
	}
	else
		left = NULL;

	if ( parentEq->right ) {
		right = new EQ_NODE(parentEq->right);
		right->parent = this;
	}
	else
		right = NULL;
}

EQ_NODE::EQ_NODE(char *s, int start, int end) {

	parent = NULL;
	left = NULL;
	right = NULL;

	char *temp;
	char *leftTemp;
	char *rightTemp;

	int leftStart;
	int leftEnd;
	int rightStart;
	int rightEnd;

	ExtractNodeOp(s,start+1);

	if ( (nodeType==VAR) || (nodeType==PARAM) ) {

		temp = simParams->StrCpy(s,start+3,end-2);
		index = atoi(temp)-1;

		delete temp;
		temp = NULL;
	}
	
	else if ( (nodeType>=PLUS) && (nodeType<=MULT) ) {

		leftStart = start+2;
		leftEnd = simParams->FindStrEnd(s,leftStart+1,1)-1;

		leftTemp = simParams->StrCpy(s,leftStart,leftEnd);

		left = new EQ_NODE(leftTemp,0,leftEnd-leftStart+1);
		left->parent = this;

		rightStart = leftEnd+1;
		rightEnd = simParams->FindStrEnd(s,rightStart+1,1)-1;

		rightTemp = simParams->StrCpy(s,rightStart,rightEnd);

		right = new EQ_NODE(rightTemp,0,rightEnd-rightStart+1);
		right->parent = this;

		delete leftTemp;
		leftTemp = NULL;

		delete rightTemp;
		rightTemp = NULL;
	}
}

EQ_NODE::~EQ_NODE(void) {

	parent = NULL;

	if ( left ) {
		delete left;
		left = NULL;
	}

	if ( right ) {
		delete right;
		right = NULL;
	}
}

void EQ_NODE::AllCalculateTermDegrees(int currTerm, MATRIX *termDegrees) {

	if ( nodeType==VAR )
		termDegrees->Add(currTerm,int(index),1.0);

	if ( left )
		left->AllCalculateTermDegrees(currTerm,termDegrees);

	if ( right )
		right->AllCalculateTermDegrees(currTerm,termDegrees);

}

MATRIX *EQ_NODE::AllExtractTermsAndDegrees(void) {

	EQ_NODE **terms = new EQ_NODE *[100];
	int numTerms = 0;

	AllFindTerms(&numTerms,terms);

	/*
	for (int i=0;i<numTerms;i++) {
		terms[i]->Print();
		printf("\n");
	}
	*/

	MATRIX *termDegrees = new MATRIX(numTerms,simParams->numEqs,0.0);

	for (int i=0;i<numTerms;i++) {
		terms[i]->AllCalculateTermDegrees(i,termDegrees);
		terms[i] = NULL;
	}

	delete[] terms;
	terms = NULL;

	//termDegrees->Print(true);

	return( termDegrees );
}

void EQ_NODE::AllFindTerms(int *numTerms, EQ_NODE **terms) {

	if ( AllMultsBelow() ) {
		terms[(*numTerms)] = this;
		(*numTerms) = (*numTerms) +1;
	}
	else {
		left->AllFindTerms(&*numTerms,terms);
		right->AllFindTerms(&*numTerms,terms);
	}
}

int  EQ_NODE::AllMultsBelow(void) {

	if ( IsZeroParity() )
		return( true );
	else
		if ( nodeType==MULT )
			return( left->AllMultsBelow() && right->AllMultsBelow() );
		else
			return( false );
}

void EQ_NODE::AllSplit(void) {

	MutateSplit();

	if ( nodeType==MULT ) {

		MutateSwapSubtrees();

		MutateSplit();
	}

	if ( left )
		left->AllSplit();

	if ( right )
		right->AllSplit();
}

void EQ_NODE::Clean(void) {

	if ( simParams->FlipCoin() )
		MutateSwapSubtrees();

	if ( simParams->FlipCoin() )
		MutateSwapOperators();

	if ( left )
		left->Clean();

	if ( right )
		right->Clean();

	if ( nodeType==PLUS )
		CollapsePlus();

	if ( nodeType==SUB )
		CollapseSubtraction();

	if ( nodeType==MULT )
		CollapseMult();

	if ( IsTwoParity() && left->IsParam() && right->IsParam() )
		CollapseParameters();
}

int  EQ_NODE::Count(int currCount) {

	int newCount = currCount + 1;

	if ( left )
		newCount = newCount + left->Count(0);

	if ( right )
		newCount = newCount + right->Count(0);

	return( newCount );
}

double EQ_NODE::Evaluate(double currT, MATRIX *y, MATRIX *params, int *valid) {

	double val;

	switch ( nodeType ) {

		case VAR:
			val = y->vals[int(index)];
			break;

		case PARAM:
			if ( simParams->evaluatingTarget )
				val = params->vals[int(index)];
			else
				val = index;
			break;

		case PLUS:
			val = left->Evaluate(currT,y,params,&*valid) + right->Evaluate(currT,y,params,&*valid);
			break;
		case SUB:
			val = left->Evaluate(currT,y,params,&*valid) - right->Evaluate(currT,y,params,&*valid);
			break;
		case MULT:
			val = left->Evaluate(currT,y,params,&*valid) * right->Evaluate(currT,y,params,&*valid);
			break;
		case SIN:
			val = sin( left->Evaluate(currT,y,params,&*valid) );
			break;
		case COS:
			val = cos( left->Evaluate(currT,y,params,&*valid) );
			break;
	}

	if ( val < minVal )
		minVal = val;

	if ( val > maxVal )
		maxVal = val;

	return( val );
}

EQ_NODE *EQ_NODE::Find(int *counter) {

	if ( (*counter) == 0 ) {
	
		return( this );
	}

	else {

		(*counter) = (*counter) - 1;

		if ( left ) {

			EQ_NODE *tmp = left->Find( &*counter );

			if ( tmp ) 
				return( tmp );

			else
				if ( right ) {

					EQ_NODE *tmp2 = right->Find( &*counter );

					if ( tmp2 )
						return( tmp2 );
					else
						return( NULL );
				}

				else
					return( NULL );
		}
		else
			return( NULL );
	}

	return( NULL );
}

int  EQ_NODE::GetDepth(void) {

	if ( (!left) && (!right) )
		return( 1 );
	else
		if ( !right )
			return( 1 + left->GetDepth() );
		else {
			int leftDepth = left->GetDepth();
			int rightDepth = right->GetDepth();

			if ( leftDepth > rightDepth )
				return( 1 + leftDepth );
			else
				return( 1 + rightDepth );
		}

	return( 0 );
}

int  EQ_NODE::IsEqualTo(EQ_NODE *otherTree) {

	if ( nodeType != otherTree->nodeType )
		return( false );

	else {

		if ( (nodeType==VAR) || (nodeType==PARAM) )
			return( index == otherTree->index );

		else
			return( (left->IsEqualTo(otherTree->left)) && 
				    (right->IsEqualTo(otherTree->right)) );
	}
}

void EQ_NODE::Mutate(void) {

//	if ( simParams->FlipCoin() )
//		MutateSwapSubtrees();

	double mutAction = simParams->Rand(0.0,1.0);

	if ( simParams->dumbSnip ) {

		if ( mutAction < 0.5 )
			MutateDumbSnip();
		else
			MutateNode();
	}
	else {

		if ( !simParams->snip )
			MutateNode();

		else if ( !simParams->creep) {
			if ( mutAction < 0.5 )
				MutateSnip();
			else
				MutateNode();
		}
		else if ( !simParams->mutateCrunch ) {
		
			if ( mutAction < 0.33 )
				MutateCreep();

			else if ( mutAction < 0.66 )
				MutateSnip();
			else
				MutateNode();
		}

		else if ( !simParams->mutateSplit ) {

			if ( mutAction < 0.25 )
				MutateCrunch();

			else if ( mutAction < 0.5 )
				MutateCreep();
			
			else if ( mutAction < 0.75 )
				MutateSnip();

			else
				MutateNode();
		}
		else if ( !simParams->mutateMerge ) {

			if ( mutAction < 0.2 )
				MutateSplit();

			else if ( mutAction < 0.4 )
				MutateCrunch();

			else if ( mutAction < 0.6 )
				MutateCreep();
			
			else if ( mutAction < 0.8 )
				MutateSnip();

			else
				MutateNode();
		}
		else  {

			if ( mutAction < 0.1667 )
				MutateMerge();

			else if ( mutAction < (0.1667*2.0) )
				MutateSplit();

			else if ( mutAction < (0.1667*3.0) )
				MutateCrunch();

			else if ( mutAction < (0.1667*4.0) )
				MutateCreep();
			
			else if ( mutAction < (0.1667*5.0) )
				MutateSnip();

			else
				MutateNode();
		}

	}
}

void EQ_NODE::MutateNeutrally(void) {

	/*
	int mutAction = simParams->RandInt(0,3);

	switch ( mutAction ) {
		case 0:
			CollapseDivision();
			break;

		case 1:
			CollapseDoubleDiv();
			break;

		case 2:
			CollapseBigDenominator();
			break;

		case 3:
			MutateSquish();
			break;
	}
	*/

	MutateSquish();

	//AddANewTerm();
	/*
	if ( nodeType==PARAM )
		ExpandParam();
	else if ( (nodeType==POW) || (nodeType==PLUS) || (nodeType==SUB) ) {

		if ( nodeType==PLUS )
			CollapsePlus();

		if ( (nodeType==SUB) && (left->IsEqualTo(right)) )
			CollapseSubtraction();

		else if ( (left->nodeType==PARAM) && (right->nodeType==PARAM) )
			CollapseParameters();
	}

	else if ( (nodeType==MULT) || (nodeType==DIV) ) {

		if ( (left->nodeType==PARAM) && (right->nodeType==PARAM) )
			CollapseParameters();

		else {

			if ( left && (left->IsPlusMinusMultDiv()) )

				MutateNeutrallyLeft();

			else if ( right && (right->IsPlusMinusMultDiv()) )

				MutateNeutrallyRight();
		}
	}
	else if ( (nodeType==SIN) || (nodeType==COS) ) {

		if ( left->nodeType==PARAM )
			CollapseParameters();
		else
			MutateNeutrallyTrig();
	}
	*/
}

void EQ_NODE::Print(void) {

	printf("(");

	PrintNode();

	if ( left )
		left->Print();

	if ( right && (nodeType<=MULT) )
		right->Print();

	printf(")");
}

void EQ_NODE::PrintInfix(void) {

	printf("(");

	if ( (nodeType==SIN) || (nodeType==COS) ) {

		PrintNode();
		left->PrintInfix();
	}

	else {

		if ( left )
			left->PrintInfix();

		PrintNode();

		if ( right )
			right->PrintInfix();
	}

	printf(")");
}

void EQ_NODE::ResetNodeRanges(void) {

	minVal = 1000000.0;
	maxVal = -1000000.0;
}

void EQ_NODE::Save(ofstream *outFile) {

	(*outFile) << "(";

	if ( (nodeType==SIN) || (nodeType==COS) ) {

		SaveNode(outFile);

		left->Save(outFile);
	}

	else {

		if ( left )
			left->Save(outFile);

		SaveNode(outFile);

		if ( right )
			right->Save(outFile);
	}

	(*outFile) << ")";
}

// ----------------- Private Functions ------------------

void EQ_NODE::AddANewTerm(void) {

	EQ_NODE *temp = new EQ_NODE(PARAM,0.0);
	temp->nodeType = nodeType;
	temp->index = index;

	temp->left = left;
	if ( left ) {
		left = NULL;
		temp->left->parent = temp;
	}

	temp->right = right;
	if ( right ) {
		right = NULL;
		temp->right->parent = temp;
	}

	right = temp;
	temp = NULL;
	right->parent = this;

	nodeType = PLUS;

	left = new EQ_NODE(simParams->targetDepth-1);
	left->parent = this;
}

void EQ_NODE::CollapseMult(void) {

	if ( nodeType == MULT ) {

		if ( (left->nodeType  == PARAM) &&
			 (right->nodeType == MULT) &&
			 (right->left->nodeType == PARAM) ) {

			left->index = left->index * right->left->index;
			delete right->left;
			right->left = NULL;

			EQ_NODE *A = right->right;
			right->right = NULL;
			
			delete right;

			right = A;
			right->parent = this;
			A = NULL;
		}
	}
}

void EQ_NODE::CollapseParameters(void) {

	if ( nodeType==PLUS )
		index = left->index + right->index;
	
	else if ( nodeType==SUB )
		index = left->index - right->index;

	else if ( nodeType==MULT )
		index = left->index * right->index;

	nodeType = PARAM;

	if ( left ) {

		delete left;
		left = NULL;
	}

	if ( right ) {
		delete right;
		right = NULL;
	}
}

void EQ_NODE::CollapsePlus(void) {

	if ( left->IsEqualTo(right) ) {

		delete left;
		left = new EQ_NODE(PARAM,2.0);
		left->parent = this;
		nodeType = MULT;
	}
	/*
	if ( (left->nodeType==PLUS) && (right->nodeType==PARAM) ) {

		if ( left->left->nodeType==PARAM ) {

			right->index = right->index + left->left->index;
			EQ_NODE *temp = new EQ_NODE(left->right);

			delete left;
			left = temp;
			temp = NULL;
			left->parent = this;
		}
		else if ( left->right->nodeType==PARAM ) {

			right->index = right->index + left->right->index;
			EQ_NODE *temp = new EQ_NODE(left->left);

			delete left;
			left = temp;
			temp = NULL;
			left->parent = this;
		}
	}
	else if ( (right->nodeType==PLUS) && (left->nodeType==PARAM) ) {

		if ( right->left->nodeType == PARAM ) {

			left->index = left->index + right->left->index;
			EQ_NODE *temp = new EQ_NODE(right->right);

			delete right;
			right = temp;
			temp = NULL;
			right->parent = this;
		}
		else if ( right->right->nodeType == PARAM ) {

			left->index = left->index + right->right->index;
			EQ_NODE *temp = new EQ_NODE(right->left);

			delete right;
			right = temp;
			temp = NULL;
			right->parent = this;
		}
	}
	*/
}

void EQ_NODE::CollapseSubtraction(void) {

	if ( left->IsEqualTo(right) ) {

		delete left;
		left = NULL;

		delete right;
		right = NULL;

		nodeType = PARAM;
		index = 0.0;
	}
}

void EQ_NODE::ExpandParam(void) {

	double originalVal = index;
	double split = simParams->Rand(0.0,1.0);

	left = new EQ_NODE(PARAM,split * originalVal);
	left->parent = this;

	right = new EQ_NODE(PARAM,(1.0-split) * originalVal);
	right->parent = this;

	nodeType = PLUS;
	index = 0.0;
}

void EQ_NODE::ExtractNodeOp(char *s, int start) {

	if ( (s[start]=='y') )
		nodeType = VAR;
	
	else if ( (s[start]=='p') )
		nodeType = PARAM;

	else if ( (s[start]=='+') )
		nodeType = PLUS;

	else if ( (s[start]=='-') && (s[start+1]=='(') )
		nodeType = SUB;

	else if ( (s[start]=='*') )
		nodeType = MULT;
}

void EQ_NODE::InheritParams(EQ_NODE *donatingNode) {

	nodeType	= donatingNode->nodeType;
	index		= donatingNode->index;
	minVal		= donatingNode->minVal;
	maxVal		= donatingNode->maxVal;
}

int  EQ_NODE::IsParam(void) {

	return( nodeType == PARAM );
}

int  EQ_NODE::IsPlusMinusMult(void) {

	return( (nodeType == PLUS) || 
		    (nodeType == SUB)  || 
			(nodeType == MULT) );
}

int  EQ_NODE::IsTerminal(void) {

	return( (!left) && (!right) );
}

int  EQ_NODE::IsTwoParity(void) {

	return( (nodeType == PLUS) || 
		    (nodeType == SUB)  || 
			(nodeType == MULT) ||
			(nodeType == SIN) ||
			(nodeType == COS) );
}

int  EQ_NODE::IsVar(void) {

	return( nodeType == VAR );
}

int  EQ_NODE::IsZeroParity(void) {

	return( (nodeType==PARAM) || (nodeType==VAR) );

}

void EQ_NODE::MutateCreep(void) {

	// A

	EQ_NODE *newLeft  = new EQ_NODE(PARAM,1.0);
	EQ_NODE *newRight = new EQ_NODE(PARAM,1.0);

	newRight->nodeType = nodeType;
	newRight->index = index;

	newRight->left = left;
	if ( left ) {
		left = NULL;
		newRight->left->parent = newRight;
	}

	newRight->right = right;
	if ( right ) {
		right = NULL;
		newRight->right->parent = newRight;
	}

	left = newLeft;
	newLeft = NULL;
	left->parent = this;

	right = newRight;
	newRight = NULL;
	right->parent = this;

	/*
	if ( simParams->FlipCoin() ) {

		nodeType = MULT;
		left->nodeType = PARAM;
		left->index = 1.0;
	}
	else {
		nodeType = PLUS;
		left->nodeType = PARAM;
		left->index = 0.0;
	}
	*/

	nodeType = simParams->RandInt(PLUS,MULT);
	left->nodeType = simParams->RandInt(VAR,PARAM);
	
	if ( left->nodeType == PARAM )
		left->index = simParams->Rand(VAL_MIN,VAL_MAX);
		//index = double( simParams->RandInt(int(VAL_MIN),(VAL_MAX)) );
	else
		left->index = double(simParams->RandInt(0,simParams->numEqs-1));

	// TERM OP A
}

void EQ_NODE::MutateCrunch(void) {

	// A(B)(C)

	if ( IsTwoParity() ) {

		EQ_NODE *temp;

		if ( simParams->FlipCoin() ) {
			delete left;
			left = NULL;

			InheritParams(right);

			left = right->left;
			right->left = NULL;
			if ( left )
				left->parent = this;

			temp = right->right;
			right->right = NULL;

			right->parent = NULL;
			delete right;
			
			right = temp;
			temp = NULL;
			
			if ( right )
				right->parent = this;
		}
		else {
			delete right;
			right = NULL;

			InheritParams(left);

			right = left->right;
			left->right = NULL;
			if ( right )
				right->parent = this;

			temp = left->left;
			left->left = NULL;

			left->parent = NULL;
			delete left;
			
			left = temp;
			temp = NULL;
			
			if ( left )
				left->parent = this;
		}
	}

	//(B) or (C)
}

void EQ_NODE::MutateDumbSnip(void) {

	if ( left ) {
		delete left;
		left = NULL;
	}

	if ( right ) {
		delete right;
		right = NULL;
	}

	double mutAction = simParams->Rand(0.0,1.0);

	if ( mutAction < 0.5 ) {

		nodeType = PARAM;
		index = simParams->Rand(VAL_MIN,VAL_MAX);
	}
	else {

		nodeType = VAR;
		index = double( simParams->RandInt(0,simParams->numEqs-1) );
	}
}

void EQ_NODE::MutateGrow(void) {

	EQ_NODE *newBranch = new EQ_NODE(nodeType,index);

	if ( left ) {
		newBranch->left = left;
		newBranch->left->parent = newBranch;
		left = NULL;
	}
	
	if ( right ) {
		newBranch->right = right;
		newBranch->right->parent = newBranch;
		right = NULL;
	}

	nodeType = simParams->RandInt(NUM_ZERO_PARITY,NUM_OPS-1);

	if ( simParams->FlipCoin() ) {

		left = newBranch;
		newBranch = NULL;
		left->parent = this;

		right = new EQ_NODE(simParams->targetDepth-1);
		right->parent = this;
	}
	else {
		right = newBranch;
		newBranch = NULL;
		right->parent = this;

		left = new EQ_NODE(simParams->targetDepth-1);
		left->parent = this;
	}
}

void EQ_NODE::MutateMerge(void) {

	if ( ((nodeType==PLUS) || (nodeType==SUB)) &&
		 ((left->nodeType==MULT) && (right->nodeType==MULT)) ) {

		EQ_NODE *A;
		EQ_NODE *A2;
		EQ_NODE *B;
		EQ_NODE *C;

		int valid = false;

		if ( left->left->IsEqualTo(right->left) ) {

			A  = left->left;
			B  = left->right;

			A2 = right->left;
			C  = right->right;

			valid = true;
		}
		else if ( left->left->IsEqualTo(right->right) ) {

			A  = left->left;
			B  = left->right;

			C  = right->left;
			A2 = right->right;

			valid = true;
		}
		else if ( left->right->IsEqualTo(right->left) ) {

			B  = left->left;
			A  = left->right;

			A2 = right->left;
			C  = right->right;

			valid = true;
		}
		else if ( left->right->IsEqualTo(right->right) ) {

			B  = left->left;
			A  = left->right;

			C  = right->left;
			A2 = right->right;

			valid = true;
		}

		if ( valid ) {

			// AB op AC
		
			left->left = NULL;
			left->right = NULL;
			right->left = NULL;
			right->right = NULL;

			right->nodeType = nodeType;
			nodeType = MULT;
			
			delete left;
			delete A2;

			right->right = C;
			C = NULL;
			right->right->parent = right;

			right->left = B;
			B = NULL;
			right->left->parent = right;

			left = A;
			A = NULL;
			left->parent = this;

			// A(B op C)
		}
	}
}

void EQ_NODE::MutateNeutrallyLeft(void) {

	EQ_NODE *newRight  = new EQ_NODE(nodeType,0.0);
	
	newRight->right = new EQ_NODE(right);
	newRight->right->parent = newRight;

	newRight->left = left->right;
	newRight->left->parent = newRight;
	
	left->right = right;
	left->right->parent = left;

	right = newRight;
	right->parent = this;
	newRight = NULL;

	nodeType = left->nodeType;
	left->nodeType = right->nodeType;
}

void EQ_NODE::MutateNeutrallyRight(void) {

	EQ_NODE *newLeft  = new EQ_NODE(nodeType,0.0);
	
	newLeft->left = new EQ_NODE(left);
	newLeft->left->parent = newLeft;

	newLeft->right = right->left;
	newLeft->right->parent = newLeft;
	
	right->left = left;
	right->left->parent = right;

	left = newLeft;
	left->parent = this;
	newLeft = NULL;

	nodeType = right->nodeType;
	right->nodeType = left->nodeType;
}

void EQ_NODE::MutateNode(void) {

	nodeType = simParams->RandInt(0,NUM_OPS-1);

	if ( nodeType <= PARAM )
		MutateZeroParity();

	else
		MutateTwoParity();
}

void EQ_NODE::MutateOneParity(void) {

	if ( !left )
		left = new EQ_NODE(simParams->targetDepth-1);

	if ( right ) {

		delete right;
		right = NULL;
	}
}

void EQ_NODE::MutateSnip(void) {

	if ( (minVal < maxVal) && 
		 (fabs(minVal) < pow(10.0,10.0)) && 
		 (fabs(maxVal) < pow(10.0,10.0)) &&
		 (fabs(minVal) > 0) &&
		 (fabs(maxVal) > 0) ) {

		if ( left ) {
			delete left;
			left = NULL;
		}

		if ( right ) {
			delete right;
			right = NULL;
		}

		nodeType = PARAM;
		index = simParams->Rand(minVal,maxVal);
	}
	else
		MutateNode();
}

void EQ_NODE::MutateSplit(void) {

	if (  (nodeType == MULT) &&
		 ((right->nodeType==PLUS) || (right->nodeType==SUB)) ) {

		// A(B op C)

		EQ_NODE *A = left;
		EQ_NODE *B = right->left;

		left = NULL;
		right->left = NULL;

		left = new EQ_NODE(MULT,0.0);
		left->parent = this;

		nodeType = right->nodeType;
		right->nodeType = MULT;

		left->left = A;
		left->left->parent = left;
		A = NULL;
		
		left->right = B;
		left->right->parent = left;
		B = NULL;

		right->left = new EQ_NODE(left->left);
		right->left->parent = right;

		// AB op AC
	}
}

void EQ_NODE::MutateSquish(void) {

	if ( IsTwoParity() ) {

		int currSize = Count(0);

		EQ_NODE *newBranch = new EQ_NODE(0);
		int attemptNum = 0;
		int withinBounds = false;

		MATRIX *y = new MATRIX(1,simParams->numEqs,0.0);
		int valid;
		double minResult, maxResult;

		for (int j=0;j<simParams->numEqs;j++)
			y->Set(0,j,simParams->minVarValues->Get(0,j));
		minResult = newBranch->Evaluate(0,y,NULL,&valid);
		for (j=0;j<simParams->numEqs;j++)
			y->Set(0,j,simParams->maxVarValues->Get(0,j));
		maxResult = newBranch->Evaluate(0,y,NULL,&valid);

//		if ( (newBranch->Count(0)<currSize) && (minResult>=minVal) && (maxResult<=maxVal) )
		if ( newBranch->Count(0)<currSize )
			withinBounds = true;

		while ( (attemptNum<100) && (!withinBounds) ) {

			delete newBranch;
			
			newBranch = new EQ_NODE(0);

//			for (int j=0;j<simParams->numEqs;j++)
//				y->Set(0,j,simParams->minVarValues->Get(0,j));
//			minResult = newBranch->Evaluate(0,y,NULL,&valid);
//			for (j=0;j<simParams->numEqs;j++)
//				y->Set(0,j,simParams->maxVarValues->Get(0,j));
//			maxResult = newBranch->Evaluate(0,y,NULL,&valid);

//			if ( (newBranch->Count(0)<currSize) && (minResult>=minVal) && (maxResult<=maxVal) )
			if ( newBranch->Count(0)<currSize )
				withinBounds = true;

			attemptNum++;
		}

		delete y;

		if ( withinBounds ) {

			nodeType = newBranch->nodeType;
			index = newBranch->index;

			if ( left )
				delete left;

			if ( right )
				delete right;

			left = newBranch->left;
			if ( left )
				left->parent = this;
			newBranch->left = NULL;

			right = newBranch->right;
			if ( right )
				right->parent = this;
			newBranch->right = NULL;
		}

		delete newBranch;
		newBranch = NULL;
	}
}

void EQ_NODE::MutateSwapOperators(void) {

	if ( ((nodeType==MULT) && (left->nodeType==MULT)) || 
		 ((nodeType==PLUS) && (left->nodeType==PLUS)) ) {

		// ABC or A+B+C

		EQ_NODE *B = left->right;
		EQ_NODE *C = right;

		if ( C->GetDepth() <= B->GetDepth() ) {

			left->right = C;
			right = B;

			left->right->parent = left;
			right->parent = this;
		}

		B = NULL;
		C = NULL;

		// ACB or A+C+B
	}
}

void EQ_NODE::MutateSwapSubtrees(void) {

	if ( (nodeType==MULT) || (nodeType==PLUS) ) {

		// A op B

		EQ_NODE *temp = right;
		right = left;
		left = temp;
		temp = NULL;

		// B op A
	}
}

void EQ_NODE::MutateTerminal(int allowNudge) {

	if ( nodeType == VAR )
		index = simParams->RandInt(0,simParams->numEqs-1);

	else if ( nodeType == PARAM )
		MutateVal();
}

void EQ_NODE::MutateTwoParity(void) {

	if ( !left )
		left = new EQ_NODE(simParams->targetDepth-1);

	if ( !right )
		right = new EQ_NODE(simParams->targetDepth-1);
}

void EQ_NODE::MutateVal(void) {

	if ( simParams->FlipCoin() ) {

		double bias = simParams->Rand(0.0,5.0) - 5.0;

		if ( simParams->Rand(0.0,1.0) < 0.5 )
			index = index * (1.0 - exp(bias) );
		else
			index = index * (1.0 + exp(bias) );
		/*
		if ( simParams->FlipCoin() )
			index = index + 1;
		else
			index = index - 1;
		*/
	}
	else
		index = simParams->Rand(VAL_MIN,VAL_MAX);
		//index = double( simParams->RandInt(int(VAL_MIN),(VAL_MAX)) );

}

void EQ_NODE::MutateZeroParity(void) {

	if ( left ) {
		delete left;
		left = NULL;
	}

	if ( right ) {
		delete right;
		right = NULL;
	}

	MutateTerminal(false);
}

void EQ_NODE::PrintNode(void) {

	switch (nodeType) {
	case VAR:
		printf("y(%d)",int(index+1));
		break;
	case PARAM:
		printf("%3.3f",index);
		break;
	case PLUS:
		printf("+");
		break;
	case SUB:
		printf("-");
		break;
	case MULT:
		printf("*");
		break;
	case SIN:
		printf("sin");
		break;
	case COS:
		printf("cos");
		break;
	}
}

void EQ_NODE::SaveNode(ofstream *outFile) {

	switch (nodeType) {
	case VAR:
		(*outFile) << "y(" << int(index+1) << ")";
		break;
	case PARAM:
		(*outFile) << index;
		break;
	case PLUS:
		(*outFile) << "+";
		break;
	case SUB:
		(*outFile) << "-";
		break;
	case MULT:
		(*outFile) << "*";
		break;
	case SIN:
		(*outFile) << "sin";
		break;
	case COS:
		(*outFile) << "cos";
		break;
	}
}

#endif