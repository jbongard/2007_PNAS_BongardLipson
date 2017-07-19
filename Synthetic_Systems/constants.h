
#include "stdafx.h"

#ifndef _CONSTANTS_H
#define _CONSTANTS_H

const int		RANDOM_SEED			= 0;

// ----------------------------------------------------------------
//                           EE constants
// ----------------------------------------------------------------

const int		MAX_CYCLES			= 40 + 1;

const int		NUM_MODELS_FOR_TEST = 5;

const int		NUM_MODELS			= NUM_MODELS_FOR_TEST;
const int		NUM_TESTS			= NUM_MODELS;

const int		STAGNATION_LIMIT    = 20;

const int		EVALS_FOR_OBJ_FITNESS = 1000;

const double	BANKED_TEST_CUTOFF = 4.0;

const double	TERMINATION_THRESHOLD = 0.0001;

// ----------------------------------------------------------------
//                           Model constants
// ----------------------------------------------------------------

const char		TARGET_FILENAME[100] = "Targets/";

const double	TEST_MIN	= 0.0;
const double	TEST_MAX	= 1.0;

const double	VAL_MIN		= -10.0;
const double	VAL_MAX		= 10.0;

//const double	h			= 0.01;
const double	h			= 0.05;

//const int	MAX_MODEL_DEPTH = 9;
const int	MAX_MODEL_DEPTH = 6;

// ----------------------------------------------------------------
//                  Differential equation constants
// ----------------------------------------------------------------

const int		VAR		= 0;
const int		PARAM	= 1;
const int		PLUS	= 2;
const int		SUB		= 3;
const int		MULT	= 4;

const int		NUM_ZERO_PARITY = 2;
const int       NUM_OPS = 5;

const int		SOFT_NODE_WIDTH = NUM_OPS - NUM_ZERO_PARITY + 1;

// ----------------------------------------------------------------
//						Genetic algorithm constants
// ----------------------------------------------------------------

const int		GENS_PER_EPOCH = 100;

const int		STARTING_EXPERIMENT_LENGTH = 20;
//const int		STARTING_EXPERIMENT_LENGTH = 5;

const int		MAX_EST_PHASE_GENS = 1000;

#endif