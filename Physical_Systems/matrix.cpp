/* ---------------------------------------------------
   FILE:     bodyPlan.cpp
	AUTHOR:   Josh Bongard
	DATE:     October 2, 2000
	FUNCTION: This class contains all information for
				 a single physical segment of an organism
				 in the MathEngine environment.
 -------------------------------------------------- */

#include "stdafx.h"
#include "stdlib.h"
#include "math.h"

#ifndef _MATRIX_CPP
#define _MATRIX_CPP

#include "constants.h"
#include "matrix.h"
#include "simParams.h"

extern SIM_PARAMS   *simParams;

MATRIX::MATRIX(int ln, int wd) {

	length = ln;
	width = wd;

	vals = new double[length*width];

	InitRandomly();
}

MATRIX::MATRIX(int ln, int wd, double val) {

	length = ln;
	width = wd;

	vals = new double[length*width];

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++)

			Set(i,j,val);
}

MATRIX::MATRIX(MATRIX *m) {

	length = m->length;
	width  = m->width;

	vals = new double[length*width];

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++)

			Set(i,j,m->Get(i,j));
}

MATRIX::MATRIX(ifstream *inFile) {

	ReadFromFile(inFile);
}

MATRIX::MATRIX(char *fileName) {

	ifstream *inFile = new ifstream(fileName);

	ReadFromFile(inFile);

	inFile->close();
	delete inFile;
	inFile = NULL;
}

MATRIX::~MATRIX(void) {

	delete[] vals;
	vals = NULL;
}

double MATRIX::AbsoluteDifference(MATRIX *m) {

	return( AbsoluteDifference(m,width) );
}

double MATRIX::AbsoluteDifference(MATRIX *m, int startWd, int endWd) {

	if ( (length==m->length) && (width==m->width) ) {

		double sum = 0.0;
		
		for (int i=0;i<length;i++)

			for (int j=startWd;j<=endWd;j++)
				
				sum = sum + fabs( double(Get(i,j)) - double(m->Get(i,j)) );

		return( sum / double(length*(endWd-startWd+1)) );
	}
	else
		return( 0.0 );
}

double MATRIX::AbsoluteDifference(MATRIX *m, int wd) {

	if ( (length==m->length) && (width==m->width) ) {

		double sum = 0.0;
		
		for (int i=0;i<length;i++)

			for (int j=0;j<wd;j++)
				
				sum = sum + fabs( double(Get(i,j)) - double(m->Get(i,j)) );

		return( sum / double(length*wd) );
	}
	else
		return( 0.0 );
}

double MATRIX::AbsoluteDifference(MATRIX *m, int wd, double normFactor) {

	return( AbsoluteDifference(m,wd) * double(length*wd) / normFactor );
}

double MATRIX::AbsoluteDifferenceOfTails(MATRIX *m) {

	if ( length==m->length ) {

		double sum = 0.0;
		
		for (int i=0;i<length;i++)
				
			sum = sum + fabs( double(Get(i,width-1)) - double(m->Get(i,m->width-1)) );

		return( sum / double(length) );
	}
	else
		return( 0.0 );
}

void MATRIX::Add(int i, int j, double val) {

	Set(i,j,Get(i,j)+val);
}

void MATRIX::Add(double val) {

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++)

			Set(i,j,Get(i,j)+val);
}

void MATRIX::Add(MATRIX *him) {

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++)

			Set(i,j,Get(i,j)+him->Get(i,j));
}

void MATRIX::AddNoise(double amt) {

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++)

			Add(i,j,simParams->Rand(-amt,amt));
}

int  MATRIX::AllRowsAlsoIn(MATRIX *other) {

	int uniqueRowFound = false;

	int i=0;

	MATRIX *myRow;
	MATRIX *hisRow;

	while ( (!uniqueRowFound) && (i<length) ) {

		myRow = GetRow(i);

		int foundAMatch = false;

		int i2 = 0;

		while ( (!foundAMatch) && (i2<other->length) ) {

			hisRow = other->GetRow(i2);

			if ( myRow->Equal(hisRow) )
				foundAMatch = true;

			delete hisRow;

			i2++;
		}

		if ( !foundAMatch )
			uniqueRowFound = true;

		delete myRow;

		i++;
	}

	return( !uniqueRowFound );
}

void MATRIX::ChangePairInRow(int i) {

	if ( simParams->Rand(0.0,1.0) < 0.5 )

		NudgePairInRow(i);

	else
		ReplacePairInRow(i);
}

int  MATRIX::ComputeFSM(MATRIX *str, MATRIX *states, int i, int strIndex, int strLength, int q0, MATRIX *F) {

	return( 0 );
}

void MATRIX::ComputeVariances(MATRIX *vars) {

	double var;

	for (int j=0;j<width;j++) {

		var = VarOfColumn(j);

		vars->Set(0,j,var);
	}
}

void MATRIX::Copy(MATRIX *him) {

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++)

			Set(i,j,him->Get(i,j));
}

void MATRIX::CopyRow(int myI, MATRIX *him, int hisI) {

	for (int j=0;j<width;j++)

		Set(myI,j, him->Get(hisI,j) );
}

MATRIX *MATRIX::CountValues(int i1, int i2, int j1, int j2, int maxVal) {

	return( 0 );
}

void MATRIX::CreateChain(int depth) {

	int i;
	int depthCounter;

	for (int j=0;j<width;j++) {

		depthCounter = depth;
		i = j;

		while ( depthCounter > 0 ) {

			Set(i,j,1);
			i = i + 1;
			if ( i==length )
				i = 0;
			depthCounter--;
		}
	}
}

void MATRIX::CreateIdentity(void) {

	for (int i=0;i<length;i++)
		for (int j=0;j<width;j++)
			if ( i == j )
				Set(i,j,1);
			else
				Set(i,j,0);
}

void MATRIX::CreateParity(void) {

}

void MATRIX::CreatePermutation(int min, int max) {

	int j;
	int newVal;

	for (j=0;j<width;j++)
		Set(0,j,min-1);

	for (j=0;j<width;j++) {

		newVal = simParams->RandInt(min,max);

		while ( In(newVal) )
			newVal = simParams->RandInt(min,max);

		Set(0,j,newVal);
	}

}

void MATRIX::CrossRow(int myRow, MATRIX *him, int hisRow) {

	int crossPoint = simParams->RandInt(1,width-1);
	double temp;

	for (int j=crossPoint;j<width;j++) {

		temp = Get(myRow,j);
		Set(myRow,j, him->Get(hisRow,j) );
		him->Set(hisRow,j,temp);
	}
}

void MATRIX::CrossSubTrees(int i, MATRIX *him) {

	SwapRow(i,him);

	if ( (2*i+2) < length ) {
		CrossSubTrees(2*i+1,him);
		CrossSubTrees(2*i+2,him);
	}
}

void MATRIX::DecreaseFullColumns(void) {

	for (int j=0;j<width;j++) {

		if ( SumOfColumn(j) == length )
			Set(simParams->RandInt(0,length-1),j,0);
	}
}

int MATRIX::Equal(MATRIX *other) {

	int equal = true;
	int i=0;
	int j=0;

	while ( (i<length) && (equal) ) {

		if ( Get(i,j) != other->Get(i,j) )
			equal = false;
		else {
			j++;
			if ( j==width ) {
				j = 0;
				i++;
			}
		}
	}

	return( equal );
}

double MATRIX::EqualColumnVals(int col1, int col2) {

	double sumOfEquals = 0.0;

	for (int i=0;i<length;i++)

		if ( Get(i,col1) == Get(i,col2) )
			sumOfEquals++;

	return( sumOfEquals / double(length) );
}

int    MATRIX::Exploding(void) {

	int exploding = false;

	int i=0;
	int j;

	while ( (!exploding) && (i<length) ) {

		j = 0;

		while ( (!exploding) && (j<width) ) {

			if ( fabs(Get(i,j)) > pow(10.0,6.0) )
				exploding = true;

			j++;
		}
		i++;
	}

	return( exploding );
}

double MATRIX::Evaluate(int i, double currT, MATRIX *y, MATRIX *params) {

	if ( (2*i+2) < length )
		return( EvaluateNonTerminal(i,currT,y,params) );
	
	else
		return( EvaluateTerminal(i,currT,y,params) );
}

void MATRIX::FillColumn(int myCol, int hisCol, MATRIX *m) {

	for (int i=0;i<length;i++)
		Set(i,myCol,m->Get(i,hisCol));
}

void MATRIX::FillRow(int myRow, int hisRow, MATRIX *m) {

	for (int j=0;j<width;j++)
		Set(myRow,j,m->Get(hisRow,j));
}

void MATRIX::Flip(int i, int j) {

	if ( Get(i,j) )
		Set(i,j,0);
	else
		Set(i,j,1);
}

void MATRIX::FlipRandomBit(void) {

	Flip(simParams->RandInt(0,length-1),simParams->RandInt(0,width-1));
}

void MATRIX::FlipRandomBitInRow(int i) {

	Flip(i,simParams->RandInt(0,width-1));
}


double MATRIX::Get(int i, int j) {

	return( vals[i*width + j] );
}

MATRIX *MATRIX::GetColumn(int j) {

	MATRIX *column = new MATRIX(length,1);

	for (int currRow=0;currRow<length;currRow++)
		column->Set(currRow,0,Get(currRow,j));

	return( column );
}

void MATRIX::GetMaxValsInColumn(int j, MATRIX *maxVals) {

	double maxVal = -1000.0;

	int i;
	int k;

	for (i=0;i<length;i++)

		if ( Get(i,j) > maxVal )
			maxVal = Get(i,j);

	maxVals->Set(0,0,maxVal);

	for (k=1;k<maxVals->width;k++) {

		maxVal = -1000.0;

		for (i=0;i<length;i++) {

			if ( (Get(i,j) > maxVal) && !(maxVals->In(Get(i,j))) )
				maxVal = Get(i,j);
		}

		maxVals->Set(0,k,maxVal);
	}
}

void MATRIX::GetMaxValsInColumn(int j, int iFirst, int iLast, MATRIX *maxVals) {

	double maxVal = -1000.0;

	int i;
	int k;

	for (i=iFirst;i<=iLast;i++)

		if ( Get(i,j) > maxVal )
			maxVal = Get(i,j);

	maxVals->Set(0,0,maxVal);

	for (k=1;k<maxVals->width;k++) {

		maxVal = -1000.0;

		for (i=iFirst;i<=iLast;i++) {

			if ( (Get(i,j) > maxVal) && !(maxVals->In(Get(i,j))) )
				maxVal = Get(i,j);
		}

		maxVals->Set(0,k,maxVal);
	}
}

MATRIX *MATRIX::GetRow(int i) {

	MATRIX *row = new MATRIX(1,width);

	for (int currColumn=0;currColumn<width;currColumn++)
		row->Set(0,currColumn,Get(i,currColumn));

	return( row );
}

int  MATRIX::In(double val) {

	int found = false;

	int i=0;
	int j;

	while ( (i<length) && (!found) ) {

		j = 0;

		while ( (j<width) && (!found) ) {

			if ( Get(i,j) == val )
				found = true;

			j++;
		}

		i++;
	}

	return( found );
}

void MATRIX::IncreaseEmptyColumns(void) {

	for (int j=0;j<width;j++) {

		if ( SumOfColumn(j) == 0 )
			Set(simParams->RandInt(0,length-1),j,1);
	}
}

void MATRIX::InitColumns(int colSum) {

	SetAllTo(0);

	for (int j=0;j<width;j++)
		InitColumn(j,colSum);
}

void MATRIX::InsertValSomewhereInEachRow(double val) {

	for (int i=0;i<length;i++)
		Set(i,simParams->RandInt(0,width-1),val);
}

double MATRIX::MaxDiffCol(int myCol, MATRIX *him, int hisCol) {

	double max = -1.0;
	double diff;

	for (int i=0;i<length;i++) {
	
		diff = fabs( Get(i,myCol) - him->Get(i,hisCol) );

		if ( diff > max )
			max = diff;
		else if ( diff < 0 )
			max = pow(10.0,10.0);
		}

	return( max );
}

double MATRIX::MaxDiffColChange(int myCol, MATRIX *him, int hisCol) {

	double max = -1.0;
	double diff;

	for (int i=0;i<length-1;i++) {
	
		diff = fabs( (Get(i,myCol)-Get(i+1,myCol)) - (him->Get(i,hisCol)-him->Get(i+1,hisCol)) );

		if ( diff > max )
			max = diff;
		else if ( diff < 0 )
			max = pow(10.0,10.0);
		}

	return( max );
}

double MATRIX::MaxDiff(MATRIX *other) {

	double max = -1.0;
	double diff;

	int i,j;

	for (i=0;i<length;i++)
		for (j=0;j<width;j++) {
	
			diff = fabs( Get(i,j) - other->Get(i,j) );

			if ( diff > max )
				max = diff;
			else if ( diff < 0 )
				max = pow(10.0,10.0);
		}

	return( max );
}

double MATRIX::MaxValInColumn(int j) {

	double maxVal = -1000.0;

	for (int i=0;i<length;i++)

		if ( Get(i,j) > maxVal )

			maxVal = Get(i,j);

	return( maxVal );
}

double MATRIX::Mean(void) {

	double mean = 0.0;

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++)

			mean = mean + Get(i,j);

	return( mean / double(length*width) );
}

double MATRIX::MeanOfColumn(int j) {

	double mean = 0.0;

	for (int i=0;i<length;i++)
		mean = mean + Get(i,j);

	return( mean / double(length) );
}

double MATRIX::Min(void) {

	double minVal = 1000000.0;

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++)

			if ( Get(i,j) < minVal )

				minVal = Get(i,j);

	return( minVal );
}

double MATRIX::MinValInColumn(int j) {

	double minVal = 1000.0;

	for (int i=0;i<length;i++)

		if ( Get(i,j) < minVal )

			minVal = Get(i,j);

	return( minVal );
}

int  MATRIX::MostSimilarRow(MATRIX *r, int i1, int i2) {

	int maxSimilarity = 10000;
	int rowDiff;

	for (int i=i1;i<i2;i++) {

		rowDiff = RowDifference(i,0,r);

		if ( rowDiff < maxSimilarity )
			maxSimilarity = rowDiff;
	}

	return( maxSimilarity );
}

double MATRIX::MSQColumn(int myCol, MATRIX *him, int hisCol) {

	double msq = 0.0;
	double val1, val2;

	for (int i=0;i<length-1;i++) {

		val1 = Get(i,myCol);
		val2 = him->Get(i,hisCol);

		msq = msq + pow( val1 - val2 ,2.0);
	}

	return( msq / double(length-1) );
}

double MATRIX::MSQRow(int myRow, MATRIX *him, int hisRow) {

	double msq = 0.0;
	double val1, val2;

	for (int j=0;j<width;j++) {

		val1 = Get(myRow,j);
		val2 = him->Get(hisRow,j);

		msq = msq + pow( val1 - val2 ,2.0);
	}

	return( msq );
}

void MATRIX::Mult(double val) {

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++)

			Set(i,j,Get(i,j) * val);
}

void MATRIX::Mutate(double min, double max) {

	if ( simParams->Rand(0.0,1.0) < 0.5 )
		Perturb(min,max);
	else
		Nudge(min,max);
}

void MATRIX::Nudge(double min, double max) {

	int i, j;
	double change;
	
	i = simParams->RandInt(0,length-1);
	j = simParams->RandInt(0,width-1);

	change = pow(10.0, simParams->RandInt(-7,2) );

	if ( simParams->Rand(0.0,1.0) < 0.5 )

		Set(i,j, Get(i,j) + change );

	else

		Set(i,j, Get(i,j) - change );

	if ( Get(i,j) < min )
		Set(i,j,min);

	if ( Get(i,j) > max )
		Set(i,j,max);	
}

void MATRIX::Perturb(double min, double max) {

	int i, j;
	double newVal;
	
	i = simParams->RandInt(0,length-1);
	j = simParams->RandInt(0,width-1);

	newVal = simParams->Rand(min,max);

	Set(i,j,newVal);
}

void MATRIX::Print(int wait) {

	for (int i=0;i<length;i++) {

		for (int j=0;j<width;j++) {
			printf("%3.3f ",vals[i*width+j]);
		}
		printf("\n");
	}

	if ( wait )
		char ch = getchar();
}

void MATRIX::Randomize(int maxVal) {

	for (int j=0;j<width;j++)
		RandomizeColumn(j,maxVal);
}

void MATRIX::RandomizeColumn(int j, int maxVal) {

	for (int i=0;i<length;i++)
		Set(i,j,simParams->RandInt(0,maxVal));
}

void MATRIX::RandomizeRow(int i, int maxVal) {

	for (int j=0;j<width;j++)
		Set(i,j,simParams->RandInt(0,maxVal));
}

void MATRIX::Replace(MATRIX *m) {

	if ( (length==m->length) && (width==m->width) ) {

		for (int i=0;i<length;i++)

			for (int j=0;j<width;j++)

				m->Set(i,j,Get(i,j));
	}
}

double MATRIX::RollingMean(MATRIX *other, int h, int w) {

	double rm = 0.0;

	for (int j=0;j<width;j++)
		rm = rm + RollingMean(other,j,h,w);

	return( rm / double(width) );
}

int  MATRIX::RowDifference(int myRow, int hisRow, MATRIX *m) {

	return( 0 );
}

void MATRIX::Save(ofstream *outFile) {

	for (int i=0;i<length;i++) {

		for (int j=0;j<width;j++) {

			(*outFile) << vals[i*width+j] << " ";
		}

		(*outFile) << "\n";
	}
}

void MATRIX::SelectUniquelyFrom(int maxVal) {

	MATRIX *chosen = new MATRIX(1,maxVal,0);
	int j;
	int chosenVal;

	for (j=0;j<width;j++) {

		chosenVal = simParams->RandInt(0,maxVal-1);

		while ( chosen->Get(0,chosenVal) )
			chosenVal = simParams->RandInt(0,maxVal-1);

		Set(0,j,chosenVal);
		chosen->Set(0,chosenVal,1);
	}

	delete chosen;
	chosen = NULL;
}

void MATRIX::Set(int i, int j, double val) {

	vals[i*width + j] = val;
}

void MATRIX::SetAllTo(int val) {

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++)

			Set(i,j,val);
}

void MATRIX::SetRow(int i, int val) {

	for (int j=0;j<width;j++)
		Set(i,j,val);
}


double MATRIX::SumOfColumn(int j) {

	double sum = 0.0;

	for (int i=0;i<length;i++)
		sum = sum + Get(i,j);

	return( sum );
}

double MATRIX::SumOfColumn(int i1, int i2, int j) {

	double sum = 0;

	for (int i=i1;i<i2;i++)
		sum = sum + Get(i,j);

	return( sum );
}

double MATRIX::SumOfIndices(MATRIX *indices, int i, int j1, int j2) {

	double sum = 0;

	for (int j=j1;j<j2;j++)
		sum = sum + Get(0,int(indices->Get(i,j)));

	return( sum );
}

double MATRIX::SumOfRow(int i) {

	double sum = 0;

	for (int j=0;j<width;j++)
		sum = sum + Get(i,j);

	return( sum );
}

void MATRIX::SwapRow(int i, MATRIX *him) {

	double temp;

	for (int j=0;j<width;j++) {
		
		temp = Get(i,j);
		Set(i,j,him->Get(i,j));
		him->Set(i,j,temp);
	}
}

double MATRIX::VarOfColumn(int j) {

	double mean = MeanOfColumn(j);
	double var = 0.0;

	for (int i=0;i<length;i++)
		var = var + pow(Get(i,j) - mean,2.0);

	return( var / double(length) );
}

void MATRIX::Write(ofstream *outFile) {

	(*outFile) << length << " " << width << "\n";

	for (int i=0;i<length;i++) {

		for (int j=0;j<width;j++)

			(*outFile) << Get(i,j) << " ";

		(*outFile) << "\n";
	}
}

void MATRIX::WriteAndRename(char *fileName) {

}

void MATRIX::ZeroColumn(int j) {

	for (int i=0;i<length;i++)
		Set(i,j,0);
}

void MATRIX::ZeroIt(void) {

	for (int j=0;j<width;j++)
		ZeroColumn(j);
}

// ----------------------------------------------------------------
//                           Private methods
// ----------------------------------------------------------------

int  MATRIX::BinaryToDecimal(int j) {

	return( 0 );
}

int MATRIX::Contains(int val) {

	int found = false;
	int i = 0;
	int j;

	while ( (i<length) && (!found) ) {

		j = 0;

		while ( (j<width) && (!found) ) {

			if ( Get(i,j) == val )
				found = true;

			j++;
		}

		i++;
	}

	return( found );
}

MATRIX *MATRIX::DecimalToBinary(int val, int maxValue) {

	int power = 0;

	while ( (maxValue - pow(2,power)) >= 0 ) {

		power = power + 1;
	}

	MATRIX *binaryVal = new MATRIX(1,power,0);

	while ( val > 0 ) {

		power = 0;

		while ( (val - pow(2,power)) >= 0 ) {

			power = power + 1;
		}

		binaryVal->Set(0,binaryVal->width-power,1);
		val = val - int(pow(2,power-1));
	}

	return( binaryVal );
}

double MATRIX::EvaluateNonTerminal(int i, double currT, MATRIX *y, MATRIX *params) {

	double left = Evaluate(2*i+1,currT,y,params);
	double right;

	if ( (vals[i*width+2] > 0.0) ||
		 (vals[i*width+3] > 0.0) ||
		 (vals[i*width+4] > 0.0) ||
		 (vals[i*width+5] > 0.0) ||
		 (vals[i*width+6] > 0.0) )
		right = Evaluate(2*i+2,currT,y,params);

	double val = 0.0;

	if ( vals[i*width] > 0.0 )
		val = val + vals[i*width]*sin(left);

	if ( vals[i*width+1] > 0.0 )
		val = val + vals[i*width+1]*cos(left);

	if ( vals[i*width+2] > 0.0 )
		val = val + vals[i*width+2]*(left+right);

	if ( vals[i*width+3] > 0.0 )
		val = val + vals[i*width+3]*(left-right);

	if ( vals[i*width+4] > 0.0 )
		val = val + vals[i*width+4]*(left*right);

	if ( vals[i*width+5] > 0.0 )
		val = val + vals[i*width+5]*(left/right);

	if ( vals[i*width+6] > 0.0 )
		val = val + vals[i*width+6]*pow(left,right);

	if ( vals[i*width+7] > 0.0 )
		val = val + vals[i*width+7]*left;

	return( val );
}

double MATRIX::EvaluateTerminal(int i, double currT, MATRIX *y, MATRIX *params) {

	double val = 0.0;
	double index;

	if ( vals[i*width] > 0.0 ) { // Var
		
		index = floor(vals[i*width+1]*simParams->numEqs);

		if ( index == simParams->numEqs )
			index = simParams->numEqs-1;

		val = val + vals[i*width]*y->vals[int(index)];
	}

	if ( vals[i*width+2] > 0.0 ) { // Param
		
		index = floor(vals[i*width+3]*simParams->numParams);

		if ( index == simParams->numParams )
			index = simParams->numParams-1;

		val = val + vals[i*width+2]*params->vals[int(index)];
	}

	return( val );
}

int  MATRIX::FindFirstValue(int j) {

	int found = false;
	int index = 0;

	while ( (index<length) && (!found) ) {

		if ( Get(index,j) )
			found = true;
		else
			index++;
	}

	return( index );
}

void MATRIX::InitColumn(int j, int colSum) {

	int i;

	for (int c=0;c<colSum;c++) {

		i = simParams->RandInt(0,length-1);

		while ( Get(i,j) > 0 )
			i = simParams->RandInt(0,length-1);

		Add(i,j,1);
	}
}

void MATRIX::InitRandomly(void) {

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++)

			Set(i,j,simParams->Rand(0,1));
}

double MATRIX::Mean(int i, int startJ, int endJ) {

	double mean = 0.0;

	for (int j=startJ;j<=endJ;j++)
		mean = mean + Get(i,j);

	return( mean/double(endJ-startJ+1) );
}

void MATRIX::NudgePairInRow(int i) {

	int j1 = simParams->RandInt(0,width-1);
	int j2 = simParams->RandInt(0,width-1);
	double overshoot;

	while ( (Get(i,j1)==0.0) && (Get(i,j2)==0.0) ) {
		
		j1 = simParams->RandInt(0,width-1);
		j2 = simParams->RandInt(0,width-1);
	}

	Add(i,j1,0.01);
	Add(i,j2,-0.01);

	if ( Get(i,j1) > 1.0 ) {
		
		overshoot = Get(i,j1) - 1.0;

		Add(i,j1,-overshoot);
		Add(i,j2,overshoot);
	}

	if ( Get(i,j2) > 1.0 ) {
		
		overshoot = Get(i,j2) - 1.0;

		Add(i,j2,-overshoot);
		Add(i,j1,overshoot);
	}

	if ( Get(i,j1) < 0.0 ) {

		overshoot = Get(i,j1);

		Add(i,j1,-overshoot);
		Add(i,j2,overshoot);
	}

	if ( Get(i,j2) < 0.0 ) {

		overshoot = Get(i,j2);

		Add(i,j2,-overshoot);
		Add(i,j1,overshoot);
	}
}

void MATRIX::ReadFromFile(ifstream *inFile) {

	double temp;

	(*inFile) >> length;
	(*inFile) >> width;

	vals = new double[length*width];

	for (int i=0;i<length;i++)

		for (int j=0;j<width;j++) {

			(*inFile) >> temp;
			Set(i,j,temp);
		}
}

void MATRIX::ReplacePairInRow(int i) {

	int j1 = simParams->RandInt(0,width-1);
	int j2 = simParams->RandInt(0,width-1);
	double total;
	double newVal1, newVal2;

	while ( (Get(i,j1)==0.0) && (Get(i,j2)==0.0) ) {
		
		j1 = simParams->RandInt(0,width-1);
		j2 = simParams->RandInt(0,width-1);
	}

	total = Get(i,j1) + Get(i,j2);
	
	newVal1 = simParams->Rand(0.0,1.0)*total;
	newVal2 = total - newVal1;

	Set(i,j1,newVal1);
	Set(i,j2,newVal2);
}

double MATRIX::RollingMean(MATRIX *other, int j, int h, int w) {

	double rm = 0.0;
	
	/*
  	double hisSum;
	double mySum;

	int start = int( double(w-1) / 2.0 );
	int end   = start + h;

	int leftOffset  = -start;
	int rightOffset = start;

	for ( int i=start ; i<end ; i++ ) {

		hisSum = 0.0;
		mySum  = 0.0;

		for ( int k=(i+leftOffset) ; k<=(i+rightOffset) ; k++ ) {

			hisSum = hisSum + other->Get(k,j);
			mySum  = mySum  + Get(k,j);
		}

		hisSum = hisSum / double(w);
		mySum  = mySum / double(w);

		rm = rm + fabs( hisSum - mySum );

	}

  	return( rm / double(h) );

	*/

	for (int i=0;i<length;i++) {

		rm = rm + fabs( Get(i,j) - other->Get(i,j) );
	}

	return( rm / double(length) );
}

#endif