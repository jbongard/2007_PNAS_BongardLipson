// DiffEq_Inference.cpp : Infer a model of a system using differential equations.
//

#include "stdafx.h"
#include "stdlib.h"

#include "EE.h"
#include "simParams.h"

SIM_PARAMS *simParams;

void main(int argc, char* argv[])
{

	simParams = new SIM_PARAMS(argc,argv);

	system("process -p DiffEq_Inference.exe low");
	//system("process -p Prune_DumbSnip.exe low");

	EE *ee = new EE;

	ee->PerformInference();
	//ee->PerformBatch();

	delete ee;
	ee = NULL;

	delete simParams;
	simParams = NULL;
}

