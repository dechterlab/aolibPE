// main.cpp -- The main file of the library.

/*
 * Copyright (C) 2004 Robert Mateescu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

using namespace std;


#include "header/problem.h"
#include "header/bayes.h"
#include "header/ProblemMixed.h"
#include "header/defs.h"
#include "header/timers.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

long deadEnds;   // to count the deadends due to the constraints
long zeroDeadEnds; // to count the deadends due to 0's in CPT's
double resultsAOsearch[100][20][20];   //ibound x algorithm x measure
double resultsSecondTree[100][20][20];   //ibound x algorithm x measure
FILE* fp, fpsum;

bool CACHE_AT_AND_NODES		= false;
bool CACHE_AT_OR_NODES		= false;
bool CUTSET_CACHE			= false;
bool BRUTE_FORCE_W_CUTSET	= false; // to explore the cutset with linear space
bool SWITCH_TO_BE			= false; // this is to switch to bucket elimination when space permits

bool COUNTING				= false;
bool g_useSuperlinkOrder	= false;

vector<CProblemMixed*>	independentComponentProblems; // to keep the connected components of a problem
extern double globalCostFromEvidence;
variable_v topChainVariables;
int percentageDoneVariable;
double scalingFactor = 1.0;





// reads the command line arguments and sets up the problem parameters
//void initExperiment(int& N, int& C, int& R, int& P_bayes, int& K, int& N_Const, int& P_Const, double& tightness,
//					int& Algorithms, int& iStart, int& iEnd, int& iStep, int& COUNT, int argc, char* argv[],
//					char inputFile[])
//{
//	int i ;
//	for (i = 1 ; i < argc ; i++)
//	{
//		char *arg = argv[i] ;
//		//printf("\n[%s]", arg) ;
//
//		// read in switches
//		if (arg[0] == '-')
//		{
//			switch (arg[1])
//			{
//				case 'n':
//				case 'N':
//					// format is -N=20
//					// format is -Nc=10
//					{
//						char *iPar ;
//						iPar = strstr(arg, "Nc") ;
//						if (!iPar) iPar = strstr(arg, "NC") ;
//						if (!iPar) iPar = strstr(arg, "nC") ;
//						if (!iPar) iPar = strstr(arg, "nc") ;
//						if (iPar)
//						{
//							// format -Nc=10
//							iPar += 3 ;
//							N_Const = atoi(iPar) ;
//							break ;
//						}
//
//						// format is -N=10
//						N = atoi(arg+3) ;
//						break;
//					}
//
//				case 'p':
//				case 'P':
//					// format is -P=20
//					// format is -Pc=10
//					{
//						char *iPar ;
//						iPar = strstr(arg, "Pc") ;
//						if (!iPar) iPar = strstr(arg, "PC") ;
//						if (!iPar) iPar = strstr(arg, "pC") ;
//						if (!iPar) iPar = strstr(arg, "pc") ;
//						if (iPar)
//						{
//							// format -P_bayes=20
//							iPar += 3 ;
//							P_Const = atoi(iPar) ;
//							break;
//						}
//
//						// format is -P=20
//						P_bayes = atoi(arg+3) ;
//						break;
//					}
//
//				case 'k':
//				case 'K':
//					// format is -K=10
//					{
//						K = atoi(arg+3) ;
//						break;
//					}
//
//				case 't':
//				case 'T':
//					{
//						tightness = ((double) atol(arg+3)) / 100.0;
//						break;
//					}
//
//				case 'r':
//				case 'R':
//					{
//						R = atoi(arg+3);
//						break;
//					}
//
//				case 'a':
//				case 'A':
//					// define the algorithms
//					// AO_NO_PRUNING		1
//					// AO_CONSTRAINTS		2
//					// AO_FC				4
//					// AO_FC_PROJ			8
//					// AO_AUX				16
//					{
//					// format is -a=10
//					Algorithms = atoi(arg+3) ;
//					break;
//					}
//
//				case 'i':
//				case 'I':
//					// format is -iStart=2
//					// format is -iEnd=10
//					// format is -iStep=2
//					{
//						char *iPar ;
//						iPar = strstr(arg, "Start") ;
//						if (iPar)
//						{
//							iPar += 6 ;
//							iStart = atoi(iPar) ;
//							break ;
//						}
//						iPar = strstr(arg, "End") ;
//						if (iPar)
//						{
//							iPar += 4 ;
//							iEnd = atoi(iPar) ;
//							break ;
//						}
//						iPar = strstr(arg, "Step") ;
//						if (iPar)
//						{
//							iPar += 5 ;
//							iStep = atoi(iPar) ;
//							break ;
//						}
//						// format is -i=10
//						COUNT = atoi(arg+3) ;
//						break ;
//					}
//				case 'f':
//				case 'F':
//					// format is -file=file_name
//					{
//						char *iPar;
//						iPar = strstr(arg, "file") ;
//						if (iPar)
//						{
//							iPar += 5 ;
//							strcpy(inputFile, iPar);
//							break ;
//						}
//					}
//			}
//		}
//	}
//
//	if(R == 0)
//		C = 0; //this means belief part is empty (problem is pure constraint)
//	else
//		C = N - R;
//}
//
// outputs some temporary results
//void outputTempResults(FILE* fp, double tm, long exps, long andN, long orN, double prob, long deadEnd, long zeroDeadEnds)
//{
//	fprintf(fp, "\nCPUtime: %f", tm);
//	fprintf(fp, "\nNnodes: %.0f", (double) exps);
//	fprintf(fp, "\nANDnodes: %.0f", (double) andN);
//	fprintf(fp, "\nORnodes: %.0f", (double) orN);
//	fprintf(fp, "\nDead-ends: %.0f", (double) deadEnds);
//	fprintf(fp, "\nZero dead-ends: %.0f", (double) zeroDeadEnds);
//	fprintf(fp, "\nProbability of Query: %e", prob);
//	fflush(fp);
//}
//

// run instances Robert
//int runExperiment(int argc, char* argv[])
//{
//	int		N = 10;					//number of variables
//	int		C = 8;					//number of CPTs
//	int		R = 2;					//number of root nodes
//	int		P_bayes = 2;			//number of parents for bayes
//	int		K = 2;					//number of values per variable
//
//	int		N_Const = 0;			//number of constraints
//	int		P_Const = 3;			//number of variables per constraint
//	double	tightness = 0.10;		//tightness of constraints
//	int		Algorithms = 4;
//
//	int iStart =0;
//	int iEnd = 30;
//	int iStep = 1;
//
//	int COUNT = 10;
//	//=======================================================================
//	double avgWidth = 0;
//	double avgHeight = 0;
//	double totalAvgWidth = 0;
//	double totalAvgHeight = 0;
//
//
//	int BE_Count = 0;
//	int i;
//	int scopeSize;
//
//	char inputFile[256];
//
//	if (argc <= 1)
//		return 0;
//
//	// begin initialization
//	if (argc >= 1)
//	{
//		initExperiment(	N, C, R,  P_bayes,  K,  N_Const,  P_Const,  tightness,
//						Algorithms,  iStart,  iEnd,  iStep,  COUNT, argc, argv, inputFile);
//	}
//
//
//	int result;
//	result = strncmp(inputFile, "", 1);
//	if(!result) cout << "no file";
//
//
//	CProblemMixed* prob;
//	// load file if necessary
//	if (result)
//	{
//		prob = new CProblemMixed();
//
//		// CPCS
//		char *file_name;
//		file_name = strstr(inputFile, "fileE") ;
//		if (file_name)
//			prob->load2(INP_SIMPLE, inputFile);
//		else
//		{
//			file_name = strstr(inputFile, "simple");
//			if (file_name)
//				prob->load(INP_SIMPLE, inputFile);
//		}
//		N = prob->getN();
//		K = prob->getK();
//		C = prob->getC();
//	}
//
//
//	// this will keep the total, then average by pcnt
//	// first 12 are the algorithm types (see defs.h)
//	// the next 10 are in order: 0 time; 1 Nnodes; 2 NANDnodes; 3 NORnodes; 4 ProbOfQuery;
//	//					5 deadEnds 6zeroDeadEnds 7 w-cutset depth
//
//	for (int k = 0; k < COUNT; k++)
//		for (i = 0; i < NUMBER_OF_ALGORITHMS; i++)
//			for (int j = 0; j < NUMBER_OF_MEASURES; j++)
//			{
//				resultsAOsearch[k][i][j] = 0.0;
//				resultsSecondTree[k][i][j] = 0.0;
//			}
//
//	// this keeps the problem specific measures; total, then average by pcnt
//	//  0 induced-width 1 tree-height
//	int measures = 10;
//	int problems = 3;
//
//	int mixed = 0;
//	int auxiliary = 1;
//	int aux_BE = 2;
//	int problemMeasure[3][15];
//	for (int j = 0; j < problems; j++)
//		for (i = 0; i < measures; i++)
//			problemMeasure[j][i] = 0;
//
//	// Write statistics to a file; create file name
//	char OutputFile[256] = {0} ;
//	strcpy(OutputFile, "ANDORcache_") ;
//	char buffer[256] = {0};
//	strcat(OutputFile, "N") ; itoa(N, buffer, 10); strcat(OutputFile, buffer); strcat(OutputFile, "_") ;
//	strcat(OutputFile, "K") ; itoa(K, buffer, 10); strcat(OutputFile, buffer); strcat(OutputFile, "_") ;
//	strcat(OutputFile, "C") ; itoa(C, buffer, 10); strcat(OutputFile, buffer); strcat(OutputFile, "_") ;
//	if (!result)
//	{
//		strcat(OutputFile, "P") ; itoa(P_bayes, buffer, 10); strcat(OutputFile, buffer); strcat(OutputFile, "_") ;
//		strcat(OutputFile, "Nc") ; itoa(N_Const, buffer, 10); strcat(OutputFile, buffer); strcat(OutputFile, "_") ;
//		strcat(OutputFile, "Pc") ; itoa(P_Const, buffer, 10); strcat(OutputFile, buffer); strcat(OutputFile, "_") ;
//		strcat(OutputFile, "T") ; ltoa((long) floor(tightness*100), buffer, 10); strcat(OutputFile, buffer); strcat(OutputFile, "_") ;
//	}
//	strcat(OutputFile, "alg") ; ltoa(Algorithms, buffer, 10); strcat(OutputFile, buffer); strcat(OutputFile, "_") ;
//	strcat(OutputFile, "inst") ; itoa(COUNT, buffer, 10); strcat(OutputFile, buffer); strcat(OutputFile, ".txt") ;
//
////	FILE* fp = fopen(OutputFile, "w");
//	fp = fopen(OutputFile, "w");
//	if (NULL == fp) return 0;
//
//	fprintf(fp, "\n N = %d", N);
//	fprintf(fp, "\n K = %d", K);
//	fprintf(fp, "\n C = %d", C);
//	fprintf(fp, "\n P_bayes = %d", P_bayes);
//	fprintf(fp, "\n N_Const = %d", N_Const);
//	fprintf(fp, "\n P_Const = %d", P_Const);
//	fprintf(fp, "\n tightness = %f", tightness);
//	fprintf(fp, "\n Number of Instances = %d", COUNT);
//	fprintf(fp, "\n\n Algorithms run:");
//	if(Algorithms & AO_NO_PRUNING)		fprintf(fp, "\n    AND/OR NO PRUNING");
//	if(Algorithms & AO_CONSTRAINTS)		fprintf(fp, "\n    AND/OR CONSTRAINTS ONLY");
//	if(Algorithms & AO_FC)				fprintf(fp, "\n    AND/OR FC");
//	if(Algorithms & AO_FC_PROJ)			fprintf(fp, "\n    AND/OR RFC");
//	if(Algorithms & AO_AUX)				fprintf(fp, "\n    AUXILIARY AND/OR");
//	if(Algorithms & BE_AUX)				fprintf(fp, "\n    AUXILIARY BUCKET ELIMINATION");
//	if(Algorithms & OR_CONSTRAINTS)		fprintf(fp, "\n    OR CONSTRAINTS ONLY");
//	if(Algorithms & OR_FC)				fprintf(fp, "\n    OR FC");
//	if(Algorithms & OR_FC_PROJ)			fprintf(fp, "\n    OR RFC");
//	if(Algorithms & AO_ADVANCED_ADVANCED)				fprintf(fp, "\n    AO CUTSET CACHE");
//	if(Algorithms & AO_LINEAR_ADVANCED)			fprintf(fp, "\n    AO BRUTE FORCE CUTSET");
//
//
//
//	int Width = 0;
//	if (!result)
//	{
//		// create problems of the same induced width:
//		for (int pcnt = 0; pcnt<10; )
//		{
//			CProblemMixed* prob = new CProblemMixed(GRAPH_RANDOM_FIXED, N, K, C, N_Const, P_Const, tightness);
//			prob->setParamParents(P_bayes);
//			prob->setParamConnectivity(UNKNOWN);
//			prob->setParamWidth(UNKNOWN);
//			prob->setParamHeight(UNKNOWN);
//			prob->setParamSigma(UNKNOWN);
//			prob->create();
//
//			// create problem
//			if (prob->preprocess(VO_MINFILL, true))
//			{
//				++pcnt;
//				Width += prob->getInducedWidth();
//				delete prob;
//			}
//			else
//			{
//				if (prob) delete prob;
//				continue;
//			}
//		}
//
//		Width = Width / 10;
//
//		if (iEnd > Width) iEnd = Width;
//		if (iEnd * (log((double) K)/log(2.0)) > 20) iEnd = floor (20 * (log(2.0)/log((double) K)) );
//
//		//iStart = iEnd;
//
//
//		//if (Width * (log(K)/log(2)) > 25 ) return 0;
//	}
//
//
//
//
//	//start running the algorithms
//	BE_Count = 0;
//
//	for (int pcnt = 0; pcnt < COUNT; )
//	{
////		CProblemMixed* prob;
//		// load file
//		if (result)
//		{
//			pcnt = COUNT;
////			prob = new CProblemMixed();
//
//			// CPCS
////			char *file_name;
////			file_name = strstr(inputFile, "fileE") ;
////			if (file_name)
////				prob->load2(INP_SIMPLE, inputFile);
////			else
////			{
////				file_name = strstr(inputFile, "cpcs");
////				if (file_name)
////					prob->load(INP_SIMPLE, inputFile);
////			}
//
//			// CPCS
//			// 	prob->load(INP_SIMPLE, "C:\\Robert\\Research\\Benchmarks\\CPCS\\cpcs54.simple");
//
//			// GENETIC LINKAGE
//		//	prob->load2(INP_SIMPLE, "C:\\Robert\\Research\\Benchmarks\\GENETIC_LINKAGE\\BayesXml\\EA\\simple\\comp01_fileEA7.simple");
//
//
//		//	prob->extractDeterminism();
//
//			if (!prob->preprocess(VO_MINFILL, true))
//				exit(1);
//			if ( iEnd > prob->getInducedWidth() )
//				iEnd = prob->getInducedWidth();
//
//			problemMeasure[mixed][0] += prob->getInducedWidth();
//			problemMeasure[mixed][1] += prob->getTreeHeight();
//		}
//
//		else
//		{
//			prob = new CProblemMixed(GRAPH_RANDOM_FIXED, N, K, C, N_Const, P_Const, tightness);
//			prob->setParamParents(P_bayes);
//			prob->setParamConnectivity(UNKNOWN);
//			prob->setParamWidth(UNKNOWN);
//			prob->setParamHeight(UNKNOWN);
//			prob->setParamSigma(UNKNOWN);
//			prob->create();
//
//			int width = 0;
//			int height = 0;
//
//			// create problem
//			if (prob->preprocess(VO_MINFILL, true))
//			{
//				//++pcnt;
//
//				width = prob->getInducedWidth();
//				height = prob->getTreeHeight();
//			}
//			else
//			{
//				if (prob) delete prob;
//				continue;
//			}
//
//			// if((width <= 0) || (width * (log(K)/log(2)) > 25 ))
//			if (width != Width)
//			{
//				if(prob) delete prob;
//				continue;
//			}
//
//			++pcnt;
//			cout << "\nInstance " << pcnt
//				 << " induced width: " << width
//				 << " rooted tree height: " << height
//			;
//			//prob->print();
//
//			problemMeasure[mixed][0] += width;
//			problemMeasure[mixed][1] += height;
//			fprintf(fp, "\n\n== Instance %d ===================================================================", pcnt);
//			fprintf(fp, "\n Induced width: %d\n Tree height: %d\n", width, height);
//			fflush(fp);
//		}
//
//		double tm;
//		long exps = 0;
//		long andN = 0, orN = 0;
//
//		//for (scopeSize = iStart; scopeSize <= iEnd; scopeSize += ((iStep < (iEnd-scopeSize)) || (scopeSize == iEnd)) ? iStep : (iEnd-scopeSize) )
//		for (scopeSize = iEnd; scopeSize >= iStart; scopeSize -= ((iStep < (scopeSize-iStart)) || (scopeSize == iStart)) ? iStep : (scopeSize-iStart) )
//		//scopeSize = Width/2;
//		//scopeSize = -1;
//		//while(scopeSize <= Width+1)
//		{
//		//	scopeSize += iStep ;
//
//			fprintf(fp, "\n\n Maximum cache size = %d ---------------- \n", scopeSize);
//			fflush(fp);
//
//			// alg1: AO(i)
//			if(Algorithms & AO_NO_PRUNING)
//			{
//				CACHE_AT_AND_NODES		= true;
//				CACHE_AT_OR_NODES		= true;
//				CUTSET_CACHE			= false;
//				BRUTE_FORCE_W_CUTSET	= false;
//				SWITCH_TO_BE			= false;
//
//				fprintf(fp, "\n AO(i) - OLD");
//				if(scopeSize >=0) prob->cacheInit(scopeSize);
//				prob->execAO(-1, tm, exps, andN, orN, false, PRUNING_NO, scopeSize);
//				if(scopeSize >=0) prob->cacheDestroy();
//				resultsAOsearch[scopeSize][PRUNING_NO][0] += tm;
//				resultsAOsearch[scopeSize][PRUNING_NO][1] += exps;
//				resultsAOsearch[scopeSize][PRUNING_NO][2] += andN;
//				resultsAOsearch[scopeSize][PRUNING_NO][3] += orN;
//				resultsAOsearch[scopeSize][PRUNING_NO][4] += prob->getProbabilityOfQuery();
//				resultsAOsearch[scopeSize][PRUNING_NO][5] += deadEnds;
//				resultsAOsearch[scopeSize][PRUNING_NO][6] += zeroDeadEnds;
//				outputTempResults(fp, tm, exps, andN, orN, prob->getProbabilityOfQuery(), deadEnds, zeroDeadEnds);
//			}
//
//			// alg2: linear cutset (brute force) + advanced caching for the rest (full caching)
//			if(Algorithms & AO_LINEAR_ADVANCED)
//			{
//				CACHE_AT_AND_NODES		= true;
//				CACHE_AT_OR_NODES		= true;
//				CUTSET_CACHE			= true;
//				BRUTE_FORCE_W_CUTSET	= true;
//				SWITCH_TO_BE			= false;
//
//				prob->setAdvancedCacheFlag(scopeSize);
//				prob->wCutsetInit();
//
//				fprintf(fp, "\n AO LINEAR - ADVANCED");
//				if(scopeSize >=0) prob->cacheInit(scopeSize);
//				prob->execAO(-1, tm, exps, andN, orN, false, PRUNING_NO, scopeSize);
//				if(scopeSize >=0) prob->cacheDestroy();
//				resultsAOsearch[scopeSize][AO_LINEAR_ADVANCED_ID][0] += tm;
//				resultsAOsearch[scopeSize][AO_LINEAR_ADVANCED_ID][1] += exps;
//				resultsAOsearch[scopeSize][AO_LINEAR_ADVANCED_ID][2] += andN;
//				resultsAOsearch[scopeSize][AO_LINEAR_ADVANCED_ID][3] += orN;
//				resultsAOsearch[scopeSize][AO_LINEAR_ADVANCED_ID][4] += prob->getProbabilityOfQuery();
//				resultsAOsearch[scopeSize][AO_LINEAR_ADVANCED_ID][5] += deadEnds;
//				resultsAOsearch[scopeSize][AO_LINEAR_ADVANCED_ID][6] += zeroDeadEnds;
//				outputTempResults(fp, tm, exps, andN, orN, prob->getProbabilityOfQuery(), deadEnds, zeroDeadEnds);
//			}
//
//			// alg2-bis: linear cutset (brute force) + advanced caching for the rest (full caching)
//			if(Algorithms & AO_LINEAR_BE)
//			{
//				CACHE_AT_AND_NODES		= true;
//				CACHE_AT_OR_NODES		= true;
//				CUTSET_CACHE			= true;
//				BRUTE_FORCE_W_CUTSET	= true;
//				SWITCH_TO_BE			= true;
//
//				prob->setAdvancedCacheFlag(scopeSize);
//				prob->wCutsetInit();
//
//				fprintf(fp, "\n AO LINEAR - BE");
//				if(scopeSize >=0) prob->cacheInit(scopeSize);
//				prob->execAO(-1, tm, exps, andN, orN, false, PRUNING_NO, scopeSize);
//				if(scopeSize >=0) prob->cacheDestroy();
//				resultsAOsearch[scopeSize][AO_LINEAR_BE_ID][0] += tm;
//				resultsAOsearch[scopeSize][AO_LINEAR_BE_ID][1] += exps;
//				resultsAOsearch[scopeSize][AO_LINEAR_BE_ID][2] += andN;
//				resultsAOsearch[scopeSize][AO_LINEAR_BE_ID][3] += orN;
//				resultsAOsearch[scopeSize][AO_LINEAR_BE_ID][4] += prob->getProbabilityOfQuery();
//				resultsAOsearch[scopeSize][AO_LINEAR_BE_ID][5] += deadEnds;
//				resultsAOsearch[scopeSize][AO_LINEAR_BE_ID][6] += zeroDeadEnds;
//				outputTempResults(fp, tm, exps, andN, orN, prob->getProbabilityOfQuery(), deadEnds, zeroDeadEnds);
//			}
//
//			// alg3: advanced cutset cache everywhere
//			if(Algorithms & AO_ADVANCED_ADVANCED)
//			{
//				CACHE_AT_AND_NODES		= true;
//				CACHE_AT_OR_NODES		= true;
//				CUTSET_CACHE			= true;
//				BRUTE_FORCE_W_CUTSET	= false;
//				SWITCH_TO_BE			= false;
//
//				prob->setAdvancedCacheFlag(scopeSize);
//				prob->wCutsetInit();
//
//				fprintf(fp, "\n AO ADVANCED - ADVANCED");
//				if(scopeSize >=0) prob->cacheInit(scopeSize);
//				prob->execAO(-1, tm, exps, andN, orN, false, PRUNING_NO, scopeSize);
//				if(scopeSize >=0) prob->cacheDestroy();
//				resultsAOsearch[scopeSize][AO_ADVANCED_ADVANCED_ID][0] += tm;
//				resultsAOsearch[scopeSize][AO_ADVANCED_ADVANCED_ID][1] += exps;
//				resultsAOsearch[scopeSize][AO_ADVANCED_ADVANCED_ID][2] += andN;
//				resultsAOsearch[scopeSize][AO_ADVANCED_ADVANCED_ID][3] += orN;
//				resultsAOsearch[scopeSize][AO_ADVANCED_ADVANCED_ID][4] += prob->getProbabilityOfQuery();
//				resultsAOsearch[scopeSize][AO_ADVANCED_ADVANCED_ID][5] += deadEnds;
//				resultsAOsearch[scopeSize][AO_ADVANCED_ADVANCED_ID][6] += zeroDeadEnds;
//				outputTempResults(fp, tm, exps, andN, orN, prob->getProbabilityOfQuery(), deadEnds, zeroDeadEnds);
//			}
//
//			// alg4: advanced cache on cutset + switch to BE
//			// need to count the number of tuples in the separators of BE
//			if(Algorithms & AO_ADVANCED_BE)
//			{
//				CACHE_AT_AND_NODES		= true;
//				CACHE_AT_OR_NODES		= true;
//				CUTSET_CACHE			= true;
//				BRUTE_FORCE_W_CUTSET	= false;
//				SWITCH_TO_BE			= true;
//
//				prob->setAdvancedCacheFlag(scopeSize);
//				prob->wCutsetInit();
////				prob->print();
//
////				cout << endl << scopeSize << " - cutset depth: ";
////				cout << prob->getWCutsetHeight()[scopeSize];
//
//				fprintf(fp, "\n %d - cutset depth: %d", scopeSize, prob->getWCutsetHeight()[scopeSize]);
//				resultsAOsearch[scopeSize][AO_LINEAR_ADVANCED_ID][7] += prob->getWCutsetHeight()[scopeSize];
//
//				fprintf(fp, "\n AO ADVANCED - BE");
//				if(scopeSize >=0) prob->cacheInit(scopeSize);
//				prob->execAO(-1, tm, exps, andN, orN, false, PRUNING_NO, scopeSize);
//				if(scopeSize >=0) prob->cacheDestroy();
//				resultsAOsearch[scopeSize][AO_ADVANCED_BE_ID][0] += tm;
//				resultsAOsearch[scopeSize][AO_ADVANCED_BE_ID][1] += exps;
//				resultsAOsearch[scopeSize][AO_ADVANCED_BE_ID][2] += andN;
//				resultsAOsearch[scopeSize][AO_ADVANCED_BE_ID][3] += orN;
//				resultsAOsearch[scopeSize][AO_ADVANCED_BE_ID][4] += prob->getProbabilityOfQuery();
//				resultsAOsearch[scopeSize][AO_ADVANCED_BE_ID][5] += deadEnds;
//				resultsAOsearch[scopeSize][AO_ADVANCED_BE_ID][6] += zeroDeadEnds;
//				outputTempResults(fp, tm, exps, andN, orN, prob->getProbabilityOfQuery(), deadEnds, zeroDeadEnds);
//			}
//
//			if(Algorithms & AO_CONSTRAINTS)
//			{
//				fprintf(fp, "\n\n AND/OR CONSTRAINTS ONLY");
//				if(scopeSize >=0) prob->cacheInit(scopeSize);
//				prob->execAO(-1, tm, exps, andN, orN, false, PRUNING_CONSTRAINTS_ONLY, scopeSize);
//				if(scopeSize >=0) prob->cacheDestroy();
//				resultsAOsearch[scopeSize][PRUNING_CONSTRAINTS_ONLY][0] += tm;
//				resultsAOsearch[scopeSize][PRUNING_CONSTRAINTS_ONLY][1] += exps;
//				resultsAOsearch[scopeSize][PRUNING_CONSTRAINTS_ONLY][2] += andN;
//				resultsAOsearch[scopeSize][PRUNING_CONSTRAINTS_ONLY][3] += orN;
//				resultsAOsearch[scopeSize][PRUNING_CONSTRAINTS_ONLY][4] += prob->getProbabilityOfQuery();
//				resultsAOsearch[scopeSize][PRUNING_CONSTRAINTS_ONLY][5] += deadEnds;
//				resultsAOsearch[scopeSize][PRUNING_CONSTRAINTS_ONLY][6] += zeroDeadEnds;
//				outputTempResults(fp, tm, exps, andN, orN, prob->getProbabilityOfQuery(), deadEnds, zeroDeadEnds);
//			}
//
//			if(Algorithms & AO_FC)
//			{
//				fprintf(fp, "\n\n AND/OR FC");
//				if(scopeSize >=0) prob->cacheInit(scopeSize);
//				prob->execAO(-1, tm, exps, andN, orN, false, PRUNING_FORWARD_CHECKING, scopeSize);
//				if(scopeSize >=0) prob->cacheDestroy();
//				resultsAOsearch[scopeSize][PRUNING_FORWARD_CHECKING][0] += tm;
//				resultsAOsearch[scopeSize][PRUNING_FORWARD_CHECKING][1] += exps;
//				resultsAOsearch[scopeSize][PRUNING_FORWARD_CHECKING][2] += andN;
//				resultsAOsearch[scopeSize][PRUNING_FORWARD_CHECKING][3] += orN;
//				resultsAOsearch[scopeSize][PRUNING_FORWARD_CHECKING][4] += prob->getProbabilityOfQuery();
//				resultsAOsearch[scopeSize][PRUNING_FORWARD_CHECKING][5] += deadEnds;
//				resultsAOsearch[scopeSize][PRUNING_FORWARD_CHECKING][6] += zeroDeadEnds;
//				outputTempResults(fp, tm, exps, andN, orN, prob->getProbabilityOfQuery(), deadEnds, zeroDeadEnds);
//			}
//
//			if(Algorithms & AO_FC_PROJ)
//				{
//				fprintf(fp, "\n\n AND/OR RFC");
//				if(scopeSize >=0) prob->cacheInit(scopeSize);
//				prob->execAO(-1, tm, exps, andN, orN, false, PRUNING_FC_PROJECTION, scopeSize);
//				if(scopeSize >=0) prob->cacheDestroy();
//				resultsAOsearch[scopeSize][PRUNING_FC_PROJECTION][0] += tm;
//				resultsAOsearch[scopeSize][PRUNING_FC_PROJECTION][1] += exps;
//				resultsAOsearch[scopeSize][PRUNING_FC_PROJECTION][2] += andN;
//				resultsAOsearch[scopeSize][PRUNING_FC_PROJECTION][3] += orN;
//				resultsAOsearch[scopeSize][PRUNING_FC_PROJECTION][4] += prob->getProbabilityOfQuery();
//				resultsAOsearch[scopeSize][PRUNING_FC_PROJECTION][5] += deadEnds;
//				resultsAOsearch[scopeSize][PRUNING_FC_PROJECTION][6] += zeroDeadEnds;
//				outputTempResults(fp, tm, exps, andN, orN, prob->getProbabilityOfQuery(), deadEnds, zeroDeadEnds);
//			}
//		}
//
//
//
//		// now run the same problem and algorithms but with a different pseudo tree
////=================================================================================================//
////*
//
//		int scope = 1;
//		cout << "====== min depth pseudo tree ================================";
//		// create a tree decomposition
//		prob->makeTreeDecomposition(VO_MINFILL);
////			CTreeDecomposition* treeDecomp = prob->getTreeDecomposition();
////			treeDecomp->computeSeparatorWeights();
//
////		prob->getTreeDecomposition()->print();
//
//		for (int j = 0; j < 25; j++)
//		{
//			int temp = 0;
//			for (int i = 0 ; i < prob->getN(); i++)
//			{
//				if (prob->getWCutset()[i] >= j) temp++;
//			}
//			cout << "\n nr de noduri in " << j << "-cutset:" << temp;
//		}
//
//		prob->makeAndOrCutsetTree(scope);
//
////		prob->makeGWCutsetTree(scope);
//
////		prob->makeCutsetOr(scope);
//
//
////		CLegalTree* w_cutsetTree = prob->getTreeCutset();
////		int* w_values = w_cutsetTree->getWCutsetHeight();
//
////		for (int i = 30; i>0 ; i--)
////			cout << "\n w=" << i << " - cutset has depth: " << w_values[i-1] + 1 << " ---- f(" << i << ") = " << w_values[i-1] + i + 1 ;
//
////		prob->getTreeCutset()->print();
////		prob->getTreeDecomposition()->print();
//
//		CACHE_AT_AND_NODES = true;
//		CACHE_AT_OR_NODES = true;
//		prob->preprocessCutsetTree();
//
//		prob->wCutsetInit();
//
//		int temp = 0;
//		for ( i = 0 ; i < prob->getN(); i++)
//		{
//			if (prob->getWCutset()[i] >= scope) temp++;
//		}
//		cout << "nr de noduri in " << scope << "-cutset:" << temp << "\n";
//
//
////*/
////=================================================================================================//
//
//		//for (scopeSize = iStart; scopeSize <= iEnd; scopeSize += ((iStep < (iEnd-scopeSize)) || (scopeSize == iEnd)) ? iStep : (iEnd-scopeSize) )
//		for (scopeSize = iEnd; scopeSize >= iStart; scopeSize -= ((iStep < (scopeSize-iStart)) || (scopeSize == iStart)) ? iStep : (scopeSize-iStart) )
//		{
//			// alg1: AO(i)
//			if(Algorithms & AO_NO_PRUNING)
//			{
//				CACHE_AT_AND_NODES		= true;
//				CACHE_AT_OR_NODES		= true;
//				CUTSET_CACHE			= false;
//				BRUTE_FORCE_W_CUTSET	= false;
//				SWITCH_TO_BE			= false;
//
//				fprintf(fp, "\n AO(i) - OLD");
//				if(scopeSize >=0) prob->cacheInit(scopeSize);
//				prob->execAO(-1, tm, exps, andN, orN, false, PRUNING_NO, scopeSize);
//				if(scopeSize >=0) prob->cacheDestroy();
//				resultsSecondTree[scopeSize][PRUNING_NO][0] += tm;
//				resultsSecondTree[scopeSize][PRUNING_NO][1] += exps;
//				resultsSecondTree[scopeSize][PRUNING_NO][2] += andN;
//				resultsSecondTree[scopeSize][PRUNING_NO][3] += orN;
//				resultsSecondTree[scopeSize][PRUNING_NO][4] += prob->getProbabilityOfQuery();
//				resultsSecondTree[scopeSize][PRUNING_NO][5] += deadEnds;
//				resultsSecondTree[scopeSize][PRUNING_NO][6] += zeroDeadEnds;
//				outputTempResults(fp, tm, exps, andN, orN, prob->getProbabilityOfQuery(), deadEnds, zeroDeadEnds);
//			}
//
//			// alg2: linear cutset (brute force) + advanced caching for the rest (full caching)
//			if(Algorithms & AO_LINEAR_ADVANCED)
//			{
//				CACHE_AT_AND_NODES		= true;
//				CACHE_AT_OR_NODES		= true;
//				CUTSET_CACHE			= true;
//				BRUTE_FORCE_W_CUTSET	= true;
//				SWITCH_TO_BE			= false;
//
//				prob->setAdvancedCacheFlag(scopeSize);
//				prob->wCutsetInit();
////				prob->extractDeterminism();
//
////				prob->print();
//				fprintf(fp, "\n AO LINEAR - ADVANCED");
//				if(scopeSize >=0) prob->cacheInit(scopeSize);
//				prob->execAO(-1, tm, exps, andN, orN, false, PRUNING_NO, scopeSize);
//				if(scopeSize >=0) prob->cacheDestroy();
//				resultsSecondTree[scopeSize][AO_LINEAR_ADVANCED_ID][0] += tm;
//				resultsSecondTree[scopeSize][AO_LINEAR_ADVANCED_ID][1] += exps;
//				resultsSecondTree[scopeSize][AO_LINEAR_ADVANCED_ID][2] += andN;
//				resultsSecondTree[scopeSize][AO_LINEAR_ADVANCED_ID][3] += orN;
//				resultsSecondTree[scopeSize][AO_LINEAR_ADVANCED_ID][4] += prob->getProbabilityOfQuery();
//				resultsSecondTree[scopeSize][AO_LINEAR_ADVANCED_ID][5] += deadEnds;
//				resultsSecondTree[scopeSize][AO_LINEAR_ADVANCED_ID][6] += zeroDeadEnds;
//				outputTempResults(fp, tm, exps, andN, orN, prob->getProbabilityOfQuery(), deadEnds, zeroDeadEnds);
//			}
//
//			// alg2 - bis: linear cutset (brute force) + advanced caching for the rest (full caching)
//			if(Algorithms & AO_LINEAR_BE)
//			{
//				CACHE_AT_AND_NODES		= true;
//				CACHE_AT_OR_NODES		= true;
//				CUTSET_CACHE			= true;
//				BRUTE_FORCE_W_CUTSET	= true;
//				SWITCH_TO_BE			= true;
//
//				prob->setAdvancedCacheFlag(scopeSize);
//				prob->wCutsetInit();
//
//				fprintf(fp, "\n AO LINEAR - BE");
//				if(scopeSize >=0) prob->cacheInit(scopeSize);
//				prob->execAO(-1, tm, exps, andN, orN, false, PRUNING_NO, scopeSize);
//				if(scopeSize >=0) prob->cacheDestroy();
//				resultsSecondTree[scopeSize][AO_LINEAR_BE_ID][0] += tm;
//				resultsSecondTree[scopeSize][AO_LINEAR_BE_ID][1] += exps;
//				resultsSecondTree[scopeSize][AO_LINEAR_BE_ID][2] += andN;
//				resultsSecondTree[scopeSize][AO_LINEAR_BE_ID][3] += orN;
//				resultsSecondTree[scopeSize][AO_LINEAR_BE_ID][4] += prob->getProbabilityOfQuery();
//				resultsSecondTree[scopeSize][AO_LINEAR_BE_ID][5] += deadEnds;
//				resultsSecondTree[scopeSize][AO_LINEAR_BE_ID][6] += zeroDeadEnds;
//				outputTempResults(fp, tm, exps, andN, orN, prob->getProbabilityOfQuery(), deadEnds, zeroDeadEnds);
//			}
//
//			// alg3: advanced cutset cache everywhere
//			if(Algorithms & AO_ADVANCED_ADVANCED)
//			{
//				CACHE_AT_AND_NODES		= true;
//				CACHE_AT_OR_NODES		= true;
//				CUTSET_CACHE			= true;
//				BRUTE_FORCE_W_CUTSET	= false;
//				SWITCH_TO_BE			= false;
//
//				prob->setAdvancedCacheFlag(scopeSize);
//				prob->wCutsetInit();
//
//				fprintf(fp, "\n AO ADVANCED - ADVANCED");
//				if(scopeSize >=0) prob->cacheInit(scopeSize);
//				prob->execAO(-1, tm, exps, andN, orN, false, PRUNING_NO, scopeSize);
//				if(scopeSize >=0) prob->cacheDestroy();
//				resultsSecondTree[scopeSize][AO_ADVANCED_ADVANCED_ID][0] += tm;
//				resultsSecondTree[scopeSize][AO_ADVANCED_ADVANCED_ID][1] += exps;
//				resultsSecondTree[scopeSize][AO_ADVANCED_ADVANCED_ID][2] += andN;
//				resultsSecondTree[scopeSize][AO_ADVANCED_ADVANCED_ID][3] += orN;
//				resultsSecondTree[scopeSize][AO_ADVANCED_ADVANCED_ID][4] += prob->getProbabilityOfQuery();
//				resultsSecondTree[scopeSize][AO_ADVANCED_ADVANCED_ID][5] += deadEnds;
//				resultsSecondTree[scopeSize][AO_ADVANCED_ADVANCED_ID][6] += zeroDeadEnds;
//				outputTempResults(fp, tm, exps, andN, orN, prob->getProbabilityOfQuery(), deadEnds, zeroDeadEnds);
//			}
//
//
//			// alg4: advanced cache on cutset + switch to BE
//			// need to count the number of tuples in the separators of BE
//			if(Algorithms & AO_ADVANCED_BE)
//			{
//				CACHE_AT_AND_NODES		= true;
//				CACHE_AT_OR_NODES		= true;
//				CUTSET_CACHE			= true;
//				BRUTE_FORCE_W_CUTSET	= false;
//				SWITCH_TO_BE			= true;
//
//				prob->setAdvancedCacheFlag(scopeSize);
//				prob->wCutsetInit();
////				prob->print();
//
//				cout << endl << scopeSize << " - cutset depth: ";
//				cout << prob->getWCutsetHeight()[scopeSize];
//
//				fprintf(fp, "\n %d - cutset depth: %d", scopeSize, prob->getWCutsetHeight()[scopeSize]);
//				resultsSecondTree[scopeSize][AO_LINEAR_ADVANCED_ID][7] += prob->getWCutsetHeight()[scopeSize];
//
//				fprintf(fp, "\n AO ADVANCED - BE");
//				if(scopeSize >=0) prob->cacheInit(scopeSize);
//				prob->execAO(-1, tm, exps, andN, orN, false, PRUNING_NO, scopeSize);
//				if(scopeSize >=0) prob->cacheDestroy();
//				resultsSecondTree[scopeSize][AO_ADVANCED_BE_ID][0] += tm;
//				resultsSecondTree[scopeSize][AO_ADVANCED_BE_ID][1] += exps;
//				resultsSecondTree[scopeSize][AO_ADVANCED_BE_ID][2] += andN;
//				resultsSecondTree[scopeSize][AO_ADVANCED_BE_ID][3] += orN;
//				resultsSecondTree[scopeSize][AO_ADVANCED_BE_ID][4] += prob->getProbabilityOfQuery();
//				resultsSecondTree[scopeSize][AO_ADVANCED_BE_ID][5] += deadEnds;
//				resultsSecondTree[scopeSize][AO_ADVANCED_BE_ID][6] += zeroDeadEnds;
//				outputTempResults(fp, tm, exps, andN, orN, prob->getProbabilityOfQuery(), deadEnds, zeroDeadEnds);
//			}
//		}
//
//
//
//		// now create the OR space (same ordering as before, only legal tree is a chain
//		if((Algorithms & OR_CONSTRAINTS) || (Algorithms & OR_FC) || (Algorithms & OR_FC_PROJ) )
//		{
//			// erase the previous links from m_tree
//			CLegalTree* m_tree = prob->getTree();
//			m_tree->eraseLinks();
//			m_tree->createChain(prob->getOrdering());
//
//			for (int var = 0; var < prob->getN(); ++var)
//				m_tree->descendants(var, ((CProblemMixed*) prob)->getTreeDescendants(var));
//
//			prob->cacheDestroy();
//			prob->parentsetDestroy();
//			prob->parentsetCreate();
//			prob->parentsetInit();
//			prob->setHighestVarInScope();
//
//			//prob->print();
//
//			int width = prob->getInducedWidth();
//			int height = prob->getTreeHeight();
//
//			cout << "\nInstance " << pcnt
//				 << " induced width: " << width
//				 << " rooted tree height: " << height
//			;
//
//			// problemMeasure[mixed][0] += width;
//			// problemMeasure[mixed][1] += height;
//
//			fprintf(fp, "\n\n -- OR SEARCH-------------------------------------");
//			fprintf(fp, "\n Induced width: %d\n Tree height: %d\n", width, height);
//			fflush(fp);
//
//			double tm;
//			long exps = 0;
//			long andN = 0, orN = 0;
//
//			// prob->print();
//
//			for (scopeSize = iStart; scopeSize <= iEnd; scopeSize += ((iStep < (iEnd-scopeSize)) || (scopeSize == iEnd)) ? iStep : (iEnd-scopeSize))
//			//scopeSize = (Width/2);
//			//scopeSize = -1;
//			//while(scopeSize <= Width+1)
//			{
//			//	scopeSize += iStep ;
//
//				fprintf(fp, "\n\n Maximum cache size = %d ----------\n", scopeSize);
//				fflush(fp);
//
//				if(Algorithms & OR_NO_PRUNING)
//				{
//					fprintf(fp, "\n\n OR NO PRUNING");
//					if(scopeSize >=0) prob->cacheInit(scopeSize);
//					prob->execAO(-1, tm, exps, andN, orN, false, PRUNING_NO, scopeSize);
//					if(scopeSize >=0) prob->cacheDestroy();
//					resultsAOsearch[scopeSize][OR_PRUNING_NO][0] += tm;
//					resultsAOsearch[scopeSize][OR_PRUNING_NO][1] += exps;
//					resultsAOsearch[scopeSize][OR_PRUNING_NO][2] += andN;
//					resultsAOsearch[scopeSize][OR_PRUNING_NO][3] += orN;
//					resultsAOsearch[scopeSize][OR_PRUNING_NO][4] += prob->getProbabilityOfQuery();
//					resultsAOsearch[scopeSize][OR_PRUNING_NO][5] += deadEnds;
//					resultsAOsearch[scopeSize][OR_PRUNING_NO][6] += zeroDeadEnds;
//					outputTempResults(fp, tm, exps, andN, orN, prob->getProbabilityOfQuery(), deadEnds, zeroDeadEnds);
//				}
//
//
//				if(Algorithms & OR_CONSTRAINTS)
//				{
//					fprintf(fp, "\n\n OR CONSTRAINTS ONLY");
//					if(scopeSize >=0) prob->cacheInit(scopeSize);
//					prob->execAO(-1, tm, exps, andN, orN, false, PRUNING_CONSTRAINTS_ONLY, scopeSize);
//					if(scopeSize >=0) prob->cacheDestroy();
//					resultsAOsearch[scopeSize][OR_PRUNING_CONSTRAINTS_ONLY][0] += tm;
//					resultsAOsearch[scopeSize][OR_PRUNING_CONSTRAINTS_ONLY][1] += exps;
//					resultsAOsearch[scopeSize][OR_PRUNING_CONSTRAINTS_ONLY][2] += andN;
//					resultsAOsearch[scopeSize][OR_PRUNING_CONSTRAINTS_ONLY][3] += orN;
//					resultsAOsearch[scopeSize][OR_PRUNING_CONSTRAINTS_ONLY][4] += prob->getProbabilityOfQuery();
//					resultsAOsearch[scopeSize][OR_PRUNING_CONSTRAINTS_ONLY][5] += deadEnds;
//					resultsAOsearch[scopeSize][OR_PRUNING_CONSTRAINTS_ONLY][6] += zeroDeadEnds;
//					outputTempResults(fp, tm, exps, andN, orN, prob->getProbabilityOfQuery(), deadEnds, zeroDeadEnds);
//				}
//
//				if(Algorithms & OR_FC)
//				{
//					fprintf(fp, "\n\n OR FC");
//					if(scopeSize >=0) prob->cacheInit(scopeSize);
//					prob->execAO(-1, tm, exps, andN, orN, false, PRUNING_FORWARD_CHECKING, scopeSize);
//					if(scopeSize >=0) prob->cacheDestroy();
//					resultsAOsearch[scopeSize][OR_PRUNING_FORWARD_CHECKING][0] += tm;
//					resultsAOsearch[scopeSize][OR_PRUNING_FORWARD_CHECKING][1] += exps;
//					resultsAOsearch[scopeSize][OR_PRUNING_FORWARD_CHECKING][2] += andN;
//					resultsAOsearch[scopeSize][OR_PRUNING_FORWARD_CHECKING][3] += orN;
//					resultsAOsearch[scopeSize][OR_PRUNING_FORWARD_CHECKING][4] += prob->getProbabilityOfQuery();
//					resultsAOsearch[scopeSize][OR_PRUNING_FORWARD_CHECKING][5] += deadEnds;
//					resultsAOsearch[scopeSize][OR_PRUNING_FORWARD_CHECKING][6] += zeroDeadEnds;
//					outputTempResults(fp, tm, exps, andN, orN, prob->getProbabilityOfQuery(), deadEnds, zeroDeadEnds);
//				}
//
//				if(Algorithms & OR_FC_PROJ)
//				{
//					fprintf(fp, "\n\n OR RFC");
//					if(scopeSize >=0) prob->cacheInit(scopeSize);
//					prob->execAO(-1, tm, exps, andN, orN, false, PRUNING_FC_PROJECTION, scopeSize);
//					if(scopeSize >=0) prob->cacheDestroy();
//					resultsAOsearch[scopeSize][OR_PRUNING_FC_PROJECTION][0] += tm;
//					resultsAOsearch[scopeSize][OR_PRUNING_FC_PROJECTION][1] += exps;
//					resultsAOsearch[scopeSize][OR_PRUNING_FC_PROJECTION][2] += andN;
//					resultsAOsearch[scopeSize][OR_PRUNING_FC_PROJECTION][3] += orN;
//					resultsAOsearch[scopeSize][OR_PRUNING_FC_PROJECTION][4] += prob->getProbabilityOfQuery();
//					resultsAOsearch[scopeSize][OR_PRUNING_FC_PROJECTION][5] += deadEnds;
//					resultsAOsearch[scopeSize][OR_PRUNING_FC_PROJECTION][6] += zeroDeadEnds;
//					outputTempResults(fp, tm, exps, andN, orN, prob->getProbabilityOfQuery(), deadEnds, zeroDeadEnds);
//				}
//			}
//		}
//
//
//
//		// now create auxiliary problem and run on it
//
//		CProblemBayes* prob_auxiliary;
//		if((Algorithms & AO_AUX) || (Algorithms & BE_AUX))
//		{
//
//			// run bucket elimination on the auxiliary problem
//			if(Algorithms & BE_AUX)
//			{
//				prob_auxiliary = new CProblemBayes(GRAPH_RANDOM_FIXED, N + N_Const, K, C + N_Const);
//				prob_auxiliary->setParamParents(UNKNOWN);
//				prob_auxiliary->setParamConnectivity(UNKNOWN);
//				prob_auxiliary->setParamWidth(UNKNOWN);
//				prob_auxiliary->setParamHeight(UNKNOWN);
//				prob_auxiliary->setParamSigma(UNKNOWN);
//
//				int width = 0;
//				int height = 0;
//				if(prob_auxiliary->initialize(prob->functions(), prob->functionsConstraint(), prob, BE_AUX) )
//				{
//					width = prob_auxiliary->getInducedWidth();
//					height = prob_auxiliary->getTreeHeight();
//
//					problemMeasure[aux_BE][0] += width;
//					problemMeasure[aux_BE][1] += height;
//
//					fprintf(fp, "\n\n -- AUXILIARY BE ---------------------------------");
//					fprintf(fp, "\n Induced width auxiliary: %d\n Tree height auxiliary: %d\n", width, height);
//
////					prob_auxiliary->print();
//				}
//
//				if((width !=0) && (width * (log((double) K)/log(2.0)) < 27 ))
//				{
//					BE_Count += 1;
//
//					double tm;
//					double result;
//					printf("\n BUCKET ELIMINATION AUXILIARY NETWORK\n");
//					prob_auxiliary->execBE_AUX(tm, result);
//
//					fprintf(fp, "\nCPUtime: %f", tm);
//					fprintf(fp, "\nProbability of Query: %e", result);
//					printf("\nCPUtime: %f", tm);
//					printf("\nProbability of Query: %7g", result);
//					fflush(fp);
//					resultsAOsearch[0][BE_AUXILIARY][0] += tm;
//					resultsAOsearch[0][BE_AUXILIARY][4] += result;
//				}
//
//				if (prob_auxiliary)	delete prob_auxiliary;
//			}
//
///*
//			// run search on the auxiliary problem
//			if(Algorithms & AO_AUX)
//			{
//				prob_auxiliary = new CProblemBayes(GRAPH_RANDOM_FIXED, N + N_Const, K, C + N_Const);
//				prob_auxiliary->setParamParents(UNKNOWN);
//				prob_auxiliary->setParamConnectivity(UNKNOWN);
//				prob_auxiliary->setParamWidth(UNKNOWN);
//				prob_auxiliary->setParamHeight(UNKNOWN);
//				prob_auxiliary->setParamSigma(UNKNOWN);
//
//				printf("\n\n AND/OR SEARCH ON AUXILIARY NETWORK\n");
//				if(prob_auxiliary->initialize(prob->functions(), prob->functionsConstraint(), prob, AO_AUX) )
//				{
//					int width = prob_auxiliary->getInducedWidth();
//					int height = prob_auxiliary->getTreeHeight();
//
//					problemMeasure[auxiliary][0] += width;
//					problemMeasure[auxiliary][1] += height;
//
//					fprintf(fp, "\n\n -- AUXILIARY BE ---------------------------------");
//					fprintf(fp, "\n Induced width auxiliary: %d\n Tree height auxiliary: %d\n", width, height);
//
////					prob_auxiliary->print();
//
//					double tm;
//					long exps = 0;
//					long andN = 0, orN = 0;
//					printf("\n AUXILIARY NETWORK\n");
//					prob_auxiliary->execAO(-1, tm, exps, andN, orN);
//					fprintf(fp, "\nCPUtime: %f", tm);
//					fprintf(fp, "\nNnodes: %.0f", (double) exps);
//					fprintf(fp, "\nANDnodes: %.0f", (double) andN);
//					fprintf(fp, "\nORnodes: %.0f", (double) orN);
//					fprintf(fp, "\nZero dead-ends: %.0f", (double) zeroDeadEnds);
//					fprintf(fp, "\nProbability of Query: %e", prob->getProbabilityOfQuery());
//					fflush(fp);
//					resultsAOsearch[scopeSize][AUXILIARY][0] += tm;
//					resultsAOsearch[scopeSize][AUXILIARY][1] += exps;
//					resultsAOsearch[scopeSize][AUXILIARY][2] += andN;
//					resultsAOsearch[scopeSize][AUXILIARY][3] += orN;
//					resultsAOsearch[scopeSize][AUXILIARY][4] += prob->getProbabilityOfQuery();
//
//					resultsAOsearch[scopeSize][AUXILIARY][6] += zeroDeadEnds;
//				}
//
//				if(prob_auxiliary) delete prob_auxiliary;
//			}
//*/
//
//		}
//
//		if(prob) delete prob;
//	}
//
//////================================================================================================================================================
//////=========== print at the end of all results file ===============================================================================================
//////================================================================================================================================================
//
//	fprintf(fp, "\n\n================================================================================");
//	fprintf(fp, "\n Total results");
//	fprintf(fp, "\n N = %d", N);
//	fprintf(fp, "\n K = %d", K);
//	fprintf(fp, "\n C = %d", C);
//	fprintf(fp, "\n P_bayes = %d", P_bayes);
//	fprintf(fp, "\n N_Const = %d", N_Const);
//	fprintf(fp, "\n P_Const = %d", P_Const);
//	fprintf(fp, "\n tightness = %f", tightness);
//	fprintf(fp, "\n Number of Instances = %d", COUNT);
//	fprintf(fp, "\n\n Algorithms run:");
//	if(Algorithms & AO_NO_PRUNING)					fprintf(fp, "\n    AND/OR NO PRUNING");
//	if(Algorithms & AO_CONSTRAINTS)					fprintf(fp, "\n    AND/OR CONSTRAINTS ONLY");
//	if(Algorithms & AO_FC)							fprintf(fp, "\n    AND/OR FC");
//	if(Algorithms & AO_FC_PROJ)						fprintf(fp, "\n    AND/OR RFC");
//	if(Algorithms & AO_AUX)							fprintf(fp, "\n    AUXILIARY AND/OR");
//	if(Algorithms & BE_AUX)							fprintf(fp, "\n    AUXILIARY BUCKET ELIMINATION");
//	if(Algorithms & OR_CONSTRAINTS)					fprintf(fp, "\n    OR CONSTRAINTS ONLY");
//	if(Algorithms & OR_FC)							fprintf(fp, "\n    OR FC");
//	if(Algorithms & OR_FC_PROJ)						fprintf(fp, "\n    OR RFC");
//	if(Algorithms & AO_ADVANCED_ADVANCED)				fprintf(fp, "\n    AO CUTSET CACHE");
//	if(Algorithms & AO_LINEAR_ADVANCED)			fprintf(fp, "\n    AO BRUTE FORCE CUTSET");
//
//	fprintf(fp, "\n -----------------------------------");
//
//	fprintf(fp, "\n\n MIXED NETWORK RESULTS:");
//	fprintf(fp, "\n Average induced width: %.0f", (double) problemMeasure[0][0] / COUNT);
//	fprintf(fp, "\n Average tree height: %.0f", (double) problemMeasure[0][1] / COUNT);
//
//
////	for (scopeSize = iStart; scopeSize <= iEnd; scopeSize += ((iStep < (iEnd-scopeSize)) || (scopeSize == iEnd)) ? iStep : (iEnd-scopeSize))
//	for (scopeSize = iEnd; scopeSize >= iStart; scopeSize -= ((iStep < (scopeSize-iStart)) || (scopeSize == iStart)) ? iStep : (scopeSize-iStart) )
//
//	//scopeSize = (Width/2);
//	//scopeSize = -1;
//	//while(scopeSize <= Width+1)
//	{
//	//	scopeSize ++ ;
//
//		fprintf(fp, "\n\n Maximum cache size = %d-----------", scopeSize);
//
//		// ==========================================ALGORITHMS AND/OR =====================================
//
//
//		// AND/OR PRUNING_NO
//		if(Algorithms & AO_NO_PRUNING)
//		{
//			fprintf(fp, "\n\n AO(i) - OLD");
//			fprintf(fp, "\nAverage CPUtime: %f", resultsAOsearch[scopeSize][PRUNING_NO][0] / COUNT);
//			fprintf(fp, "\nAverage Nnodes: %.0f", resultsAOsearch[scopeSize][PRUNING_NO][1] / COUNT);
//			fprintf(fp, "\nAverage ANDnodes: %.0f", resultsAOsearch[scopeSize][PRUNING_NO][2] / COUNT);
//			fprintf(fp, "\nAverage ORnodes: %.0f", resultsAOsearch[scopeSize][PRUNING_NO][3] / COUNT);
//			fprintf(fp, "\nNumber of deadends: %.0f", resultsAOsearch[scopeSize][PRUNING_NO][5] / COUNT);
//			fprintf(fp, "\nNumber of zero deadends: %.0f", resultsAOsearch[scopeSize][PRUNING_NO][6] / COUNT);
//			fprintf(fp, "\nAverage Probability of Query: %e", resultsAOsearch[scopeSize][PRUNING_NO][4] / COUNT);
//		}
//
//		if(Algorithms & AO_LINEAR_ADVANCED)
//		{
//			fprintf(fp, "\n\n AO LINEAR - ADVANCED");
//			fprintf(fp, "\nAverage CPUtime: %f", resultsAOsearch[scopeSize][AO_LINEAR_ADVANCED_ID][0] / COUNT);
//			fprintf(fp, "\nAverage Nnodes: %.0f", resultsAOsearch[scopeSize][AO_LINEAR_ADVANCED_ID][1] / COUNT);
//			fprintf(fp, "\nAverage ANDnodes: %.0f", resultsAOsearch[scopeSize][AO_LINEAR_ADVANCED_ID][2] / COUNT);
//			fprintf(fp, "\nAverage ORnodes: %.0f", resultsAOsearch[scopeSize][AO_LINEAR_ADVANCED_ID][3] / COUNT);
//			fprintf(fp, "\nNumber of deadends: %.0f", resultsAOsearch[scopeSize][AO_LINEAR_ADVANCED_ID][5] / COUNT);
//			fprintf(fp, "\nNumber of zero deadends: %.0f", resultsAOsearch[scopeSize][AO_LINEAR_ADVANCED_ID][6] / COUNT);
//			fprintf(fp, "\nAverage Probability of Query: %e", resultsAOsearch[scopeSize][AO_LINEAR_ADVANCED_ID][4] / COUNT);
//		}
//
//		if(Algorithms & AO_LINEAR_BE)
//		{
//			fprintf(fp, "\n\n AO LINEAR - BE");
//			fprintf(fp, "\nAverage CPUtime: %f", resultsAOsearch[scopeSize][AO_LINEAR_BE_ID][0] / COUNT);
//			fprintf(fp, "\nAverage Nnodes: %.0f", resultsAOsearch[scopeSize][AO_LINEAR_BE_ID][1] / COUNT);
//			fprintf(fp, "\nAverage ANDnodes: %.0f", resultsAOsearch[scopeSize][AO_LINEAR_BE_ID][2] / COUNT);
//			fprintf(fp, "\nAverage ORnodes: %.0f", resultsAOsearch[scopeSize][AO_LINEAR_BE_ID][3] / COUNT);
//			fprintf(fp, "\nNumber of deadends: %.0f", resultsAOsearch[scopeSize][AO_LINEAR_BE_ID][5] / COUNT);
//			fprintf(fp, "\nNumber of zero deadends: %.0f", resultsAOsearch[scopeSize][AO_LINEAR_BE_ID][6] / COUNT);
//			fprintf(fp, "\nAverage Probability of Query: %e", resultsAOsearch[scopeSize][AO_LINEAR_BE_ID][4] / COUNT);
//		}
//
//		if(Algorithms & AO_ADVANCED_ADVANCED)
//		{
//			fprintf(fp, "\n\n AO ADVANCED - ADVANCED");
//			fprintf(fp, "\nAverage CPUtime: %f", resultsAOsearch[scopeSize][AO_ADVANCED_ADVANCED_ID][0] / COUNT);
//			fprintf(fp, "\nAverage Nnodes: %.0f", resultsAOsearch[scopeSize][AO_ADVANCED_ADVANCED_ID][1] / COUNT);
//			fprintf(fp, "\nAverage ANDnodes: %.0f", resultsAOsearch[scopeSize][AO_ADVANCED_ADVANCED_ID][2] / COUNT);
//			fprintf(fp, "\nAverage ORnodes: %.0f", resultsAOsearch[scopeSize][AO_ADVANCED_ADVANCED_ID][3] / COUNT);
//			fprintf(fp, "\nNumber of deadends: %.0f", resultsAOsearch[scopeSize][AO_ADVANCED_ADVANCED_ID][5] / COUNT);
//			fprintf(fp, "\nNumber of zero deadends: %.0f", resultsAOsearch[scopeSize][AO_ADVANCED_ADVANCED_ID][6] / COUNT);
//			fprintf(fp, "\nAverage Probability of Query: %e", resultsAOsearch[scopeSize][AO_ADVANCED_ADVANCED_ID][4] / COUNT);
//		}
//
//
//		if(Algorithms & AO_ADVANCED_BE)
//		{
//			fprintf(fp, "\n\n AO ADVANCED - BE ");
//			fprintf(fp, "\nAverage CPUtime: %f", resultsAOsearch[scopeSize][AO_ADVANCED_BE_ID][0] / COUNT);
//			fprintf(fp, "\nAverage Nnodes: %.0f", resultsAOsearch[scopeSize][AO_ADVANCED_BE_ID][1] / COUNT);
//			fprintf(fp, "\nAverage ANDnodes: %.0f", resultsAOsearch[scopeSize][AO_ADVANCED_BE_ID][2] / COUNT);
//			fprintf(fp, "\nAverage ORnodes: %.0f", resultsAOsearch[scopeSize][AO_ADVANCED_BE_ID][3] / COUNT);
//			fprintf(fp, "\nNumber of deadends: %.0f", resultsAOsearch[scopeSize][AO_ADVANCED_BE_ID][5] / COUNT);
//			fprintf(fp, "\nNumber of zero deadends: %.0f", resultsAOsearch[scopeSize][AO_ADVANCED_BE_ID][6] / COUNT);
//			fprintf(fp, "\nAverage Probability of Query: %e", resultsAOsearch[scopeSize][AO_ADVANCED_BE_ID][4] / COUNT);
//		}
//
//		// AND/OR PRUNING_CONSTRAINTS_ONLY
//		if(Algorithms & AO_CONSTRAINTS)
//		{
//			fprintf(fp, "\n\n CHECKING CONSTRAINTS ONLY");
//			fprintf(fp, "\nAverage CPUtime: %f", resultsAOsearch[scopeSize][PRUNING_CONSTRAINTS_ONLY][0] / COUNT);
//			fprintf(fp, "\nAverage Nnodes: %.0f", resultsAOsearch[scopeSize][PRUNING_CONSTRAINTS_ONLY][1] / COUNT);
//			fprintf(fp, "\nAverage ANDnodes: %.0f", resultsAOsearch[scopeSize][PRUNING_CONSTRAINTS_ONLY][2] / COUNT);
//			fprintf(fp, "\nAverage ORnodes: %.0f", resultsAOsearch[scopeSize][PRUNING_CONSTRAINTS_ONLY][3] / COUNT);
//			fprintf(fp, "\nNumber of deadends: %.0f", resultsAOsearch[scopeSize][PRUNING_CONSTRAINTS_ONLY][5] / COUNT);
//			fprintf(fp, "\nNumber of zero deadends: %.0f", resultsAOsearch[scopeSize][PRUNING_CONSTRAINTS_ONLY][6] / COUNT);
//			fprintf(fp, "\nAverage Probability of Query: %e", resultsAOsearch[scopeSize][PRUNING_CONSTRAINTS_ONLY][4] / COUNT);
//		}
//
//		// AND/OR PRUNING_FORWARD_CHECKING
//		if(Algorithms & AO_FC)
//		{
//			fprintf(fp, "\n\n FORWARD CHECKING");
//			fprintf(fp, "\nAverage CPUtime: %f", resultsAOsearch[scopeSize][PRUNING_FORWARD_CHECKING][0] / COUNT);
//			fprintf(fp, "\nAverage Nnodes: %.0f", resultsAOsearch[scopeSize][PRUNING_FORWARD_CHECKING][1] / COUNT);
//			fprintf(fp, "\nAverage ANDnodes: %.0f", resultsAOsearch[scopeSize][PRUNING_FORWARD_CHECKING][2] / COUNT);
//			fprintf(fp, "\nAverage ORnodes: %.0f", resultsAOsearch[scopeSize][PRUNING_FORWARD_CHECKING][3] / COUNT);
//			fprintf(fp, "\nNumber of deadends: %.0f", resultsAOsearch[scopeSize][PRUNING_FORWARD_CHECKING][5] / COUNT);
//			fprintf(fp, "\nNumber of zero deadends: %.0f", resultsAOsearch[scopeSize][PRUNING_FORWARD_CHECKING][6] / COUNT);
//			fprintf(fp, "\nAverage Probability of Query: %e", resultsAOsearch[scopeSize][PRUNING_FORWARD_CHECKING][4] / COUNT);
//		}
//
//		// AND/OR PRUNING_FC_PROJECTION
//		if(Algorithms & AO_FC_PROJ)
//		{
//			fprintf(fp, "\n\n FC + PROJECTION");
//			fprintf(fp, "\nAverage CPUtime: %f", resultsAOsearch[scopeSize][PRUNING_FC_PROJECTION][0] / COUNT);
//			fprintf(fp, "\nAverage Nnodes: %.0f", resultsAOsearch[scopeSize][PRUNING_FC_PROJECTION][1] / COUNT);
//			fprintf(fp, "\nAverage ANDnodes: %.0f", resultsAOsearch[scopeSize][PRUNING_FC_PROJECTION][2] / COUNT);
//			fprintf(fp, "\nAverage ORnodes: %.0f", resultsAOsearch[scopeSize][PRUNING_FC_PROJECTION][3] / COUNT);
//			fprintf(fp, "\nNumber of deadends: %.0f", resultsAOsearch[scopeSize][PRUNING_FC_PROJECTION][5] / COUNT);
//			fprintf(fp, "\nNumber of zero deadends: %.0f", resultsAOsearch[scopeSize][PRUNING_FC_PROJECTION][6] / COUNT);
//			fprintf(fp, "\nAverage Probability of Query: %e", resultsAOsearch[scopeSize][PRUNING_FC_PROJECTION][4] / COUNT);
//		}
//	}
//
//
//
//	// ==========================================ALGORITHMS AUXILIARY =====================================
//	if( (Algorithms & AO_AUX) || (Algorithms & BE_AUX) )
//	{
//		fprintf(fp, "\n\n -----------------------------------");
//		fprintf(fp, "\n\n AUXILIARY NETWORK RESULTS:");
//		fprintf(fp, "\n Average induced width: %d", problemMeasure[1][0] / COUNT);
//		fprintf(fp, "\n Average tree height: %d", problemMeasure[1][1] / COUNT);
//
//		// AND OR SEARCH ON AUXILIARY NETWORK
//		if(Algorithms & AO_AUX)
//		{
//			fprintf(fp, "\n\n AND/OR SEARCH ON AUXILIARY:");
//			fprintf(fp, "\nAverage CPUtime: %f", resultsAOsearch[scopeSize][AUXILIARY][0] / COUNT);
//			fprintf(fp, "\nAverage Nnodes: %.0f", resultsAOsearch[scopeSize][AUXILIARY][1] / COUNT);
//			fprintf(fp, "\nAverage ANDnodes: %.0f", resultsAOsearch[scopeSize][AUXILIARY][2] / COUNT);
//			fprintf(fp, "\nAverage ORnodes: %.0f", resultsAOsearch[scopeSize][AUXILIARY][3] / COUNT);
//			fprintf(fp, "\nNumber of zero deadends: %.0f", resultsAOsearch[scopeSize][AUXILIARY][6] / COUNT);
//			fprintf(fp, "\nAverage Probability of Query: %e", resultsAOsearch[scopeSize][AUXILIARY][4] / COUNT);
//		}
//
//		// BUCKET ELIMINATION ON AUXILIARY NETWORK
//		if(Algorithms & BE_AUX)
//		{
//			fprintf(fp, "\n\n BUCKET ELIMINATION ON AUXILIARY:");
//			if(BE_Count > 0)
//			{
//				fprintf(fp, "\nAverage CPUtime: %f", resultsAOsearch[0][BE_AUXILIARY][0] / BE_Count);
//				fprintf(fp, "\nAverage Probability of Query: %e", resultsAOsearch[0][BE_AUXILIARY][4] / BE_Count);
//			}
//			else
//				fprintf(fp, "\n induced width too big");
//		}
//	}
//
//
//	// ==========================================ALGORITHMS OR =====================================
//	if( (Algorithms & OR_CONSTRAINTS) || (Algorithms & OR_FC) || (Algorithms & OR_FC_PROJ))
//	{
//		fprintf(fp, "\n\n -----------------------------------");
//		fprintf(fp, "\n\n OR SPACE RESULTS");
//	}
//
//
//	// OR NO PRUNING
//	if(Algorithms & OR_NO_PRUNING)
//	{
//		fprintf(fp, "\n\n OR NO PRUNING");
//		fprintf(fp, "\nAverage CPUtime: %f", resultsAOsearch[scopeSize][OR_PRUNING_NO][0] / COUNT);
//		fprintf(fp, "\nAverage Nnodes: %.0f", resultsAOsearch[scopeSize][OR_PRUNING_NO][1] / COUNT);
//		fprintf(fp, "\nAverage ANDnodes: %.0f", resultsAOsearch[scopeSize][OR_PRUNING_NO][2] / COUNT);
//		fprintf(fp, "\nAverage ORnodes: %.0f", resultsAOsearch[scopeSize][OR_PRUNING_NO][3] / COUNT);
//		fprintf(fp, "\nNumber of deadends: %.0f", resultsAOsearch[scopeSize][OR_PRUNING_NO][5] / COUNT);
//		fprintf(fp, "\nNumber of zero deadends: %.0f", resultsAOsearch[scopeSize][OR_PRUNING_NO][6] / COUNT);
//		fprintf(fp, "\nAverage Probability of Query: %e", resultsAOsearch[scopeSize][OR_PRUNING_NO][4] / COUNT);
//	}
//
//	// OR PRUNING_CONSTRAINTS_ONLY
//	if(Algorithms & OR_CONSTRAINTS)
//	{
//		fprintf(fp, "\n\n OR CONSTRAINTS ONLY");
//		fprintf(fp, "\nAverage CPUtime: %f", resultsAOsearch[scopeSize][OR_PRUNING_CONSTRAINTS_ONLY][0] / COUNT);
//		fprintf(fp, "\nAverage Nnodes: %.0f", resultsAOsearch[scopeSize][OR_PRUNING_CONSTRAINTS_ONLY][1] / COUNT);
//		fprintf(fp, "\nAverage ANDnodes: %.0f", resultsAOsearch[scopeSize][OR_PRUNING_CONSTRAINTS_ONLY][2] / COUNT);
//		fprintf(fp, "\nAverage ORnodes: %.0f", resultsAOsearch[scopeSize][OR_PRUNING_CONSTRAINTS_ONLY][3] / COUNT);
//		fprintf(fp, "\nNumber of deadends: %.0f", resultsAOsearch[scopeSize][OR_PRUNING_CONSTRAINTS_ONLY][5] / COUNT);
//		fprintf(fp, "\nNumber of zero deadends: %.0f", resultsAOsearch[scopeSize][OR_PRUNING_CONSTRAINTS_ONLY][6] / COUNT);
//		fprintf(fp, "\nAverage Probability of Query: %e", resultsAOsearch[scopeSize][OR_PRUNING_CONSTRAINTS_ONLY][4] / COUNT);
//	}
//
//	// OR PRUNING_FORWARD_CHECKING
//	if(Algorithms & OR_FC)
//	{
//		fprintf(fp, "\n\n FORWARD CHECKING");
//		fprintf(fp, "\nAverage CPUtime: %f", resultsAOsearch[scopeSize][OR_PRUNING_FORWARD_CHECKING][0] / COUNT);
//		fprintf(fp, "\nAverage Nnodes: %.0f", resultsAOsearch[scopeSize][OR_PRUNING_FORWARD_CHECKING][1] / COUNT);
//		fprintf(fp, "\nAverage ANDnodes: %.0f", resultsAOsearch[scopeSize][OR_PRUNING_FORWARD_CHECKING][2] / COUNT);
//		fprintf(fp, "\nAverage ORnodes: %.0f", resultsAOsearch[scopeSize][OR_PRUNING_FORWARD_CHECKING][3] / COUNT);
//		fprintf(fp, "\nNumber of deadends: %.0f", resultsAOsearch[scopeSize][OR_PRUNING_FORWARD_CHECKING][5] / COUNT);
//		fprintf(fp, "\nNumber of zero deadends: %.0f", resultsAOsearch[scopeSize][OR_PRUNING_FORWARD_CHECKING][6] / COUNT);
//		fprintf(fp, "\nAverage Probability of Query: %e", resultsAOsearch[scopeSize][OR_PRUNING_FORWARD_CHECKING][4] / COUNT);
//	}
//
//	// OR PRUNING_FC_PROJECTION
//	if(Algorithms & OR_FC_PROJ)
//	{
//		fprintf(fp, "\n\n FC + PROJECTION");
//		fprintf(fp, "\nAverage CPUtime: %f", resultsAOsearch[scopeSize][OR_PRUNING_FC_PROJECTION][0] / COUNT);
//		fprintf(fp, "\nAverage Nnodes: %.0f", resultsAOsearch[scopeSize][OR_PRUNING_FC_PROJECTION][1] / COUNT);
//		fprintf(fp, "\nAverage ANDnodes: %.0f", resultsAOsearch[scopeSize][OR_PRUNING_FC_PROJECTION][2] / COUNT);
//		fprintf(fp, "\nAverage ORnodes: %.0f", resultsAOsearch[scopeSize][OR_PRUNING_FC_PROJECTION][3] / COUNT);
//		fprintf(fp, "\nNumber of deadends: %.0f", resultsAOsearch[scopeSize][OR_PRUNING_FC_PROJECTION][5] / COUNT);
//		fprintf(fp, "\nNumber of zero deadends: %.0f", resultsAOsearch[scopeSize][OR_PRUNING_FC_PROJECTION][6] / COUNT);
//		fprintf(fp, "\nAverage Probability of Query: %e", resultsAOsearch[scopeSize][OR_PRUNING_FC_PROJECTION][4] / COUNT);
//	}
//
//
//	fclose(fp);
//
//
//////================================================================================================================================================
//////=========== print the summary file =============================================================================================================
//////================================================================================================================================================
//
//	char summaryOutputFile[256] = {0} ;
//	strcpy(summaryOutputFile, "summary_") ;
//	strcat(summaryOutputFile, OutputFile);
//	FILE* fpsum = fopen(summaryOutputFile, "w");
//
////	fprintf(fpsum, "\n\n================================================================================");
////	fprintf(fpsum, "Total results");
//	fprintf(fpsum, " N = %d", N);
//	fprintf(fpsum, "\n K = %d", K);
//	fprintf(fpsum, "\n C = %d", C);
//	fprintf(fpsum, "\n P_bayes = %d", P_bayes);
//	fprintf(fpsum, "\n N_Const = %d", N_Const);
//	fprintf(fpsum, "\n P_Const = %d", P_Const);
//	fprintf(fpsum, "\n tightness = %f", tightness);
//	fprintf(fpsum, "\n Number of Instances = %d", COUNT);
//	fprintf(fpsum, "\n\n Algorithms run:");
//
//	if(Algorithms & AO_NO_PRUNING)					fprintf(fpsum, "\n    AND/OR NO PRUNING");
//	if(Algorithms & AO_CONSTRAINTS)					fprintf(fpsum, "\n    AND/OR CONSTRAINTS ONLY");
//	if(Algorithms & AO_FC)							fprintf(fpsum, "\n    AND/OR FC");
//	if(Algorithms & AO_FC_PROJ)						fprintf(fpsum, "\n    AND/OR RFC");
//	if(Algorithms & AO_AUX)							fprintf(fpsum, "\n    AUXILIARY AND/OR");
//	if(Algorithms & BE_AUX)							fprintf(fpsum, "\n    AUXILIARY BUCKET ELIMINATION");
//	if(Algorithms & OR_NO_PRUNING)					fprintf(fpsum, "\n    OR NO PRUNING");
//	if(Algorithms & OR_CONSTRAINTS)					fprintf(fpsum, "\n    OR CONSTRAINTS ONLY");
//	if(Algorithms & OR_FC)							fprintf(fpsum, "\n    OR FC");
//	if(Algorithms & OR_FC_PROJ)						fprintf(fpsum, "\n    OR RFC");
//	if(Algorithms & AO_LINEAR_ADVANCED)				fprintf(fpsum, "\n    AO LINEAR - ADVANCED");
//	if(Algorithms & AO_LINEAR_BE)					fprintf(fpsum, "\n    AO LINEAR - BE");
//	if(Algorithms & AO_ADVANCED_ADVANCED)			fprintf(fpsum, "\n    AO ADVANCED - ADVANCED");
//	if(Algorithms & AO_ADVANCED_BE)					fprintf(fpsum, "\n    AO ADVANCED - BE");
//
//	fprintf(fpsum, "\n -----------------------------------");
//
////	fprintf(fpsum, "\n\n MIXED NETWORK RESULTS:");
//	fprintf(fpsum, "\n Average induced width: %.0f", (double) problemMeasure[mixed][0] / COUNT);
//	fprintf(fpsum, "\n Average tree height: %.0f", (double) problemMeasure[mixed][1] / COUNT);
//
////	fprintf(fpsum, "\n\n AUXILIARY NETWORK RESULTS:");
////	fprintf(fpsum, "\n AO SEARCH Average induced width: %d", problemMeasure[auxiliary][0] / COUNT);
////	fprintf(fpsum, "\n AO SEARCH Average tree height: %d", problemMeasure[auxiliary][1] / COUNT);
//
////	fprintf(fpsum, "\n\n BUKCET ELIM Average induced width: %d", problemMeasure[aux_BE][0] / COUNT);
////	fprintf(fpsum, "\n BUCKET ELIM Average tree height: %d", problemMeasure[aux_BE][1] / COUNT);
//
//
//
//	// w-cutset depth
//	fprintf(fpsum, "\n\n\n=========================================================================================");
//	fprintf(fpsum, "\n\t\t\t\tAverage w-cutset depth:\n");
//	//	for (scopeSize = iStart; scopeSize <= iEnd; scopeSize += ((iStep < (iEnd-scopeSize)) || (scopeSize == iEnd)) ? iStep : (iEnd-scopeSize))
//	for (scopeSize = iEnd; scopeSize >= iStart; scopeSize -= ((iStep < (scopeSize-iStart)) || (scopeSize == iStart)) ? iStep : (scopeSize-iStart) )
//	{
//		if(Algorithms & AO_ADVANCED_BE)
//		{
//			fprintf(fpsum, "\n  w = %d \t depth = %3.0f %3.0f", scopeSize,	resultsAOsearch[scopeSize][AO_LINEAR_ADVANCED_ID][7] / COUNT, resultsSecondTree[scopeSize][AO_LINEAR_ADVANCED_ID][7] / COUNT);
//			fprintf(fpsum, "     ||     f = %3.0f %3.0f", resultsAOsearch[scopeSize][AO_LINEAR_ADVANCED_ID][7] / COUNT + scopeSize, resultsSecondTree[scopeSize][AO_LINEAR_ADVANCED_ID][7] / COUNT + scopeSize);
//		}
//	}
//
//
//	// CPU time
//	fprintf(fpsum, "\n\n\n=========================================================================================");
//	fprintf(fpsum, "\n\t\t\t\tAverage CPUtime (seconds):");
//	if((Algorithms & BE_AUX) && (BE_Count > 0))
//			fprintf(fpsum, "\n AUXILIARY BE \t\t\t%15.5f", resultsAOsearch[0][BE_AUXILIARY][0] / BE_Count);
//		else
//			if (Algorithms && BE_AUX) fprintf(fpsum, "\n AUXILIARY BE: w* too big");
//
//
////	for (scopeSize = iStart; scopeSize <= iEnd; scopeSize += ((iStep < (iEnd-scopeSize)) || (scopeSize == iEnd)) ? iStep : (iEnd-scopeSize))
//	for (scopeSize = iEnd; scopeSize >= iStart; scopeSize -= ((iStep < (scopeSize-iStart)) || (scopeSize == iStart)) ? iStep : (scopeSize-iStart) )
//
//	//scopeSize = iStart;
//	//while(scopeSize <= Width+1)
//	{
//	//	scopeSize += iStep ;
//		fprintf(fpsum, "\n\n--- (i = %d) ----------", scopeSize);
//		if(Algorithms & AO_NO_PRUNING)			fprintf(fpsum, "\n AO(i) - OLD \t\t\t\t%15.5f \t%15.5f", resultsAOsearch[scopeSize][PRUNING_NO][0] / COUNT, resultsSecondTree[scopeSize][PRUNING_NO][0] / COUNT);
//		if(Algorithms & AO_CONSTRAINTS)			fprintf(fpsum, "\n AND/OR CONSTRAINTS ONLY \t%15.5f",	resultsAOsearch[scopeSize][PRUNING_CONSTRAINTS_ONLY][0] / COUNT);
//		if(Algorithms & AO_FC)					fprintf(fpsum, "\n AND/OR FC \t\t\t%15.5f",				resultsAOsearch[scopeSize][PRUNING_FORWARD_CHECKING][0] / COUNT);
//		if(Algorithms & AO_FC_PROJ)				fprintf(fpsum, "\n AND/OR RFC \t\t\t%15.5f",			resultsAOsearch[scopeSize][PRUNING_FC_PROJECTION][0] / COUNT);
//		if(Algorithms & AO_AUX)					fprintf(fpsum, "\n AUXILIARY AND/OR \t\t%15.5f",		resultsAOsearch[scopeSize][AUXILIARY][0] / COUNT);
//		if(Algorithms & OR_NO_PRUNING)			fprintf(fpsum, "\n OR NO PRUNING \t\t\t%15.5f",			resultsAOsearch[scopeSize][OR_PRUNING_NO][0] / COUNT);
//		if(Algorithms & OR_CONSTRAINTS)			fprintf(fpsum, "\n OR CONSTRAINTS \t\t%15.5f",			resultsAOsearch[scopeSize][OR_PRUNING_CONSTRAINTS_ONLY][0] / COUNT);
//		if(Algorithms & OR_FC)					fprintf(fpsum, "\n OR FC \t\t\t\t%15.5f",				resultsAOsearch[scopeSize][OR_PRUNING_FORWARD_CHECKING][0] / COUNT);
//		if(Algorithms & OR_FC_PROJ)				fprintf(fpsum, "\n OR RFC \t\t\t%15.5f",				resultsAOsearch[scopeSize][OR_PRUNING_FC_PROJECTION][0] / COUNT);
//		if(Algorithms & AO_LINEAR_ADVANCED)	fprintf(fpsum, "\n AO LINEAR - ADVANCED  \t\t%15.5f \t%15.5f",	resultsAOsearch[scopeSize][AO_LINEAR_ADVANCED_ID][0] / COUNT, resultsSecondTree[scopeSize][AO_LINEAR_ADVANCED_ID][0] / COUNT);
//		if(Algorithms & AO_LINEAR_BE)	fprintf(fpsum, "\n AO LINEAR - BE \t\t%15.5f \t%15.5f",	resultsAOsearch[scopeSize][AO_LINEAR_BE_ID][0] / COUNT, resultsSecondTree[scopeSize][AO_LINEAR_BE_ID][0] / COUNT);
//		if(Algorithms & AO_ADVANCED_ADVANCED)		fprintf(fpsum, "\n AO ADVANCED - ADVANCED \t%15.5f \t%15.5f",			resultsAOsearch[scopeSize][AO_ADVANCED_ADVANCED_ID][0] / COUNT, resultsSecondTree[scopeSize][AO_ADVANCED_ADVANCED_ID][0] / COUNT);
//		if(Algorithms & AO_ADVANCED_BE) fprintf(fpsum, "\n AO ADVANCED - BE \t\t%15.5f \t%15.5f",			resultsAOsearch[scopeSize][AO_ADVANCED_BE_ID][0] / COUNT, resultsSecondTree[scopeSize][AO_ADVANCED_BE_ID][0] / COUNT);
//
//	}
//
//	// Average number of expansions:
//	fprintf(fpsum, "\n\n\n=========================================================================================");
//	fprintf(fpsum, "\n\t\t\t\tAverage number of expansions:");
//	if((Algorithms & BE_AUX) && (BE_Count > 0))
//		fprintf(fpsum, "\n AUXILIARY BE \t\t\t%15.0f", resultsAOsearch[0][BE_AUXILIARY][1] / BE_Count);
//	else fprintf(fpsum, "\n AUXILIARY BE: w* too big");
////	for (scopeSize = iStart; scopeSize <= iEnd; scopeSize += ((iStep < (iEnd-scopeSize)) || (scopeSize == iEnd)) ? iStep : (iEnd-scopeSize))
//	for (scopeSize = iEnd; scopeSize >= iStart; scopeSize -= ((iStep < (scopeSize-iStart)) || (scopeSize == iStart)) ? iStep : (scopeSize-iStart) )
//	//scopeSize = (Width/2);
//	//while(scopeSize <= Width+1)
//	{
//	//	scopeSize += iStep ;
//		fprintf(fpsum, "\n\n--- (i = %d) ----------", scopeSize);
//		if(Algorithms & AO_NO_PRUNING)	fprintf(fpsum, "\n AO(i) - OLD \t\t\t\t%15.0f \t%15.0f", resultsAOsearch[scopeSize][PRUNING_NO][1] / COUNT, resultsSecondTree[scopeSize][PRUNING_NO][1] / COUNT);
//		if(Algorithms & AO_CONSTRAINTS)	fprintf(fpsum, "\n AND/OR CONSTRAINTS ONLY \t%15.0f", resultsAOsearch[scopeSize][PRUNING_CONSTRAINTS_ONLY][1] / COUNT);
//		if(Algorithms & AO_FC)			fprintf(fpsum, "\n AND/OR FC \t\t\t%15.0f", resultsAOsearch[scopeSize][PRUNING_FORWARD_CHECKING][1] / COUNT);
//		if(Algorithms & AO_FC_PROJ)		fprintf(fpsum, "\n AND/OR RFC \t\t\t%15.0f", resultsAOsearch[scopeSize][PRUNING_FC_PROJECTION][1] / COUNT);
//		if(Algorithms & AO_AUX)			fprintf(fpsum, "\n AUXILIARY AND/OR \t\t%15.0f", resultsAOsearch[scopeSize][AUXILIARY][1] / COUNT);
//		if(Algorithms & OR_NO_PRUNING)	fprintf(fpsum, "\n OR NO PRUNING \t\t\t%15.0f", resultsAOsearch[scopeSize][OR_PRUNING_NO][1] / COUNT);
//		if(Algorithms & OR_CONSTRAINTS)	fprintf(fpsum, "\n OR CONSTRAINTS \t\t%15.0f", resultsAOsearch[scopeSize][OR_PRUNING_CONSTRAINTS_ONLY][1] / COUNT);
//		if(Algorithms & OR_FC)			fprintf(fpsum, "\n OR FC \t\t\t\t%15.0f", resultsAOsearch[scopeSize][OR_PRUNING_FORWARD_CHECKING][1] / COUNT);
//		if(Algorithms & OR_FC_PROJ)		fprintf(fpsum, "\n OR RFC \t\t\t%15.0f", resultsAOsearch[scopeSize][OR_PRUNING_FC_PROJECTION][1] / COUNT);
//		if(Algorithms & AO_LINEAR_ADVANCED)	fprintf(fpsum, "\n AO LINEAR - ADVANCED  \t\t%15.0f \t%15.0f",	resultsAOsearch[scopeSize][AO_LINEAR_ADVANCED_ID][1] / COUNT, resultsSecondTree[scopeSize][AO_LINEAR_ADVANCED_ID][1] / COUNT);
//		if(Algorithms & AO_LINEAR_BE)	fprintf(fpsum, "\n AO LINEAR - BE \t\t%15.0f \t%15.0f",	resultsAOsearch[scopeSize][AO_LINEAR_BE_ID][1] / COUNT, resultsSecondTree[scopeSize][AO_LINEAR_BE_ID][1] / COUNT);
//		if(Algorithms & AO_ADVANCED_ADVANCED)		fprintf(fpsum, "\n AO ADVANCED - ADVANCED \t%15.0f \t%15.0f", resultsAOsearch[scopeSize][AO_ADVANCED_ADVANCED_ID][1] / COUNT, resultsSecondTree[scopeSize][AO_ADVANCED_ADVANCED_ID][1] / COUNT);
//		if(Algorithms & AO_ADVANCED_BE) fprintf(fpsum, "\n AO ADVANCED - BE \t\t%15.0f \t%15.0f",			resultsAOsearch[scopeSize][AO_ADVANCED_BE_ID][1] / COUNT, resultsSecondTree[scopeSize][AO_ADVANCED_BE_ID][1] / COUNT);
//
//	}
//
//
//	// Average number of AND nodes
//	fprintf(fpsum, "\n\n\n=========================================================================================");
//	fprintf(fpsum, "\n\t\t\t\tAverage number of AND nodes:");
//	if((Algorithms & BE_AUX) && (BE_Count > 0))
//		fprintf(fpsum, "\n AUXILIARY BE \t\t\t%15.0f", resultsAOsearch[0][BE_AUXILIARY][2] / BE_Count);
//	else fprintf(fpsum, "\n AUXILIARY BE: w* too big");
////	for (scopeSize = iStart; scopeSize <= iEnd; scopeSize += ((iStep < (iEnd-scopeSize)) || (scopeSize == iEnd)) ? iStep : (iEnd-scopeSize))
//	for (scopeSize = iEnd; scopeSize >= iStart; scopeSize -= ((iStep < (scopeSize-iStart)) || (scopeSize == iStart)) ? iStep : (scopeSize-iStart) )
//	//scopeSize = (Width/2);
//	//while(scopeSize <= Width+1)
//	{
//	//	scopeSize += iStep ;
//		fprintf(fpsum, "\n\n--- (i = %d) ----------", scopeSize);
//		if(Algorithms & AO_NO_PRUNING)	fprintf(fpsum, "\n AO(i) - OLD\t\t\t\t%15.0f  \t%15.0f", resultsAOsearch[scopeSize][PRUNING_NO][2] / COUNT, resultsSecondTree[scopeSize][PRUNING_NO][2] / COUNT);
//		if(Algorithms & AO_CONSTRAINTS)	fprintf(fpsum, "\n AND/OR CONSTRAINTS ONLY \t%15.0f", resultsAOsearch[scopeSize][PRUNING_CONSTRAINTS_ONLY][2] / COUNT);
//		if(Algorithms & AO_FC)			fprintf(fpsum, "\n AND/OR FC \t\t\t%15.0f", resultsAOsearch[scopeSize][PRUNING_FORWARD_CHECKING][2] / COUNT);
//		if(Algorithms & AO_FC_PROJ)		fprintf(fpsum, "\n AND/OR RFC \t\t\t%15.0f", resultsAOsearch[scopeSize][PRUNING_FC_PROJECTION][2] / COUNT);
//		if(Algorithms & AO_AUX)			fprintf(fpsum, "\n AUXILIARY AND/OR \t\t%15.0f", resultsAOsearch[scopeSize][AUXILIARY][2] / COUNT);
//		if(Algorithms & OR_NO_PRUNING)	fprintf(fpsum, "\n OR NO PRUNING \t\t\t%15.0f", resultsAOsearch[scopeSize][OR_PRUNING_NO][2] / COUNT);
//		if(Algorithms & OR_CONSTRAINTS)	fprintf(fpsum, "\n OR CONSTRAINTS \t\t%15.0f", resultsAOsearch[scopeSize][OR_PRUNING_CONSTRAINTS_ONLY][2] / COUNT);
//		if(Algorithms & OR_FC)			fprintf(fpsum, "\n OR FC \t\t\t\t%15.0f", resultsAOsearch[scopeSize][OR_PRUNING_FORWARD_CHECKING][2] / COUNT);
//		if(Algorithms & OR_FC_PROJ)		fprintf(fpsum, "\n OR RFC \t\t\t%15.0f", resultsAOsearch[scopeSize][OR_PRUNING_FC_PROJECTION][2] / COUNT);
//		if(Algorithms & AO_LINEAR_ADVANCED)	fprintf(fpsum, "\n AO LINEAR - ADVANCED  \t\t%15.0f \t%15.0f",	resultsAOsearch[scopeSize][AO_LINEAR_ADVANCED_ID][2] / COUNT, resultsSecondTree[scopeSize][AO_LINEAR_ADVANCED_ID][2] / COUNT);
//		if(Algorithms & AO_LINEAR_BE)	fprintf(fpsum, "\n AO LINEAR - BE \t\t%15.0f \t%15.0f",	resultsAOsearch[scopeSize][AO_LINEAR_BE_ID][2] / COUNT, resultsSecondTree[scopeSize][AO_LINEAR_BE_ID][2] / COUNT);
//		if(Algorithms & AO_ADVANCED_ADVANCED)		fprintf(fpsum, "\n AO ADVANCED - ADVANCED \t%15.0f \t%15.0f", resultsAOsearch[scopeSize][AO_ADVANCED_ADVANCED_ID][2] / COUNT, resultsSecondTree[scopeSize][AO_ADVANCED_ADVANCED_ID][2] / COUNT);
//		if(Algorithms & AO_ADVANCED_BE) fprintf(fpsum, "\n AO ADVANCED - BE \t\t%15.0f \t%15.0f",			resultsAOsearch[scopeSize][AO_ADVANCED_BE_ID][2] / COUNT, resultsSecondTree[scopeSize][AO_ADVANCED_BE_ID][2] / COUNT);
//
//	}
//
//
//	// Average number of OR nodes
//	fprintf(fpsum, "\n\n\n=========================================================================================");
//	fprintf(fpsum, "\n\n\n\t\t\t\tAverage number of OR nodes:");
//	if((Algorithms & BE_AUX) && (BE_Count > 0))
//		fprintf(fpsum, "\n AUXILIARY BE \t\t\t%15.0f", resultsAOsearch[0][BE_AUXILIARY][3] / BE_Count);
//	else fprintf(fpsum, "\n AUXILIARY BE: w* too big");
////	for (scopeSize = iStart; scopeSize <= iEnd; scopeSize += ((iStep < (iEnd-scopeSize)) || (scopeSize == iEnd)) ? iStep : (iEnd-scopeSize))
//	for (scopeSize = iEnd; scopeSize >= iStart; scopeSize -= ((iStep < (scopeSize-iStart)) || (scopeSize == iStart)) ? iStep : (scopeSize-iStart) )
//	//scopeSize = (Width/2);
//	//while(scopeSize <= Width+1)
//	{
//	//	scopeSize += iStep ;
//		fprintf(fpsum, "\n\n--- (i = %d) ----------", scopeSize);
//		if(Algorithms & AO_NO_PRUNING)	fprintf(fpsum, "\n AO(i) - OLD\t\t\t\t%15.0f \t%15.0f", resultsAOsearch[scopeSize][PRUNING_NO][3] / COUNT, resultsSecondTree[scopeSize][PRUNING_NO][3] / COUNT);
//		if(Algorithms & AO_CONSTRAINTS)	fprintf(fpsum, "\n AND/OR CONSTRAINTS ONLY \t%15.0f", resultsAOsearch[scopeSize][PRUNING_CONSTRAINTS_ONLY][3] / COUNT);
//		if(Algorithms & AO_FC)			fprintf(fpsum, "\n AND/OR FC \t\t\t%15.0f", resultsAOsearch[scopeSize][PRUNING_FORWARD_CHECKING][3] / COUNT);
//		if(Algorithms & AO_FC_PROJ)		fprintf(fpsum, "\n AND/OR RFC \t\t\t%15.0f", resultsAOsearch[scopeSize][PRUNING_FC_PROJECTION][3] / COUNT);
//		if(Algorithms & AO_AUX)			fprintf(fpsum, "\n AUXILIARY AND/OR \t\t%15.0f", resultsAOsearch[scopeSize][AUXILIARY][3] / COUNT);
//		if(Algorithms & OR_NO_PRUNING)	fprintf(fpsum, "\n OR NO PRUNING \t\t\t%15.0f", resultsAOsearch[scopeSize][OR_PRUNING_NO][3] / COUNT);
//		if(Algorithms & OR_CONSTRAINTS)	fprintf(fpsum, "\n OR CONSTRAINTS \t\t%15.0f", resultsAOsearch[scopeSize][OR_PRUNING_CONSTRAINTS_ONLY][3] / COUNT);
//		if(Algorithms & OR_FC)			fprintf(fpsum, "\n OR FC \t\t\t\t%15.0f", resultsAOsearch[scopeSize][OR_PRUNING_FORWARD_CHECKING][3] / COUNT);
//		if(Algorithms & OR_FC_PROJ)		fprintf(fpsum, "\n OR RFC \t\t\t%15.0f", resultsAOsearch[scopeSize][OR_PRUNING_FC_PROJECTION][3] / COUNT);
//		if(Algorithms & AO_LINEAR_ADVANCED)	fprintf(fpsum, "\n AO LINEAR - ADVANCED  \t\t%15.0f \t%15.0f",	resultsAOsearch[scopeSize][AO_LINEAR_ADVANCED_ID][3] / COUNT, resultsSecondTree[scopeSize][AO_LINEAR_ADVANCED_ID][3] / COUNT);
//		if(Algorithms & AO_LINEAR_BE)	fprintf(fpsum, "\n AO LINEAR - BE \t\t%15.0f \t%15.0f",	resultsAOsearch[scopeSize][AO_LINEAR_BE_ID][3] / COUNT, resultsSecondTree[scopeSize][AO_LINEAR_BE_ID][3] / COUNT);
//		if(Algorithms & AO_ADVANCED_ADVANCED)		fprintf(fpsum, "\n AO ADVANCED - ADVANCED \t%15.0f \t%15.0f", resultsAOsearch[scopeSize][AO_ADVANCED_ADVANCED_ID][3] / COUNT, resultsSecondTree[scopeSize][AO_ADVANCED_ADVANCED_ID][3] / COUNT);
//		if(Algorithms & AO_ADVANCED_BE) fprintf(fpsum, "\n AO ADVANCED - BE \t\t%15.0f \t%15.0f",			resultsAOsearch[scopeSize][AO_ADVANCED_BE_ID][3] / COUNT, resultsSecondTree[scopeSize][AO_ADVANCED_BE_ID][3] / COUNT);
//
//	}
//
//
//	// Average number of deadends
//	fprintf(fpsum, "\n\n\n=========================================================================================");
//	fprintf(fpsum, "\n\n\n\t\t\t\tAverage number of deadends:");
//	if((Algorithms & BE_AUX) && (BE_Count > 0))
//		fprintf(fpsum, "\n AUXILIARY BE \t\t\t%15.0f", resultsAOsearch[0][BE_AUXILIARY][5] / BE_Count);
//	else fprintf(fpsum, "\n AUXILIARY BE: w* too big");
////	for (scopeSize = iStart; scopeSize <= iEnd; scopeSize += ((iStep < (iEnd-scopeSize)) || (scopeSize == iEnd)) ? iStep : (iEnd-scopeSize))
//	for (scopeSize = iEnd; scopeSize >= iStart; scopeSize -= ((iStep < (scopeSize-iStart)) || (scopeSize == iStart)) ? iStep : (scopeSize-iStart) )
//	//scopeSize = (Width/2);
//	//while(scopeSize <= Width+1)
//	{
//	//	scopeSize += iStep ;
//		fprintf(fpsum, "\n\n--- (i = %d) ----------", scopeSize);
//		if(Algorithms & AO_NO_PRUNING)	fprintf(fpsum, "\n AO(i) - OLD\t\t\t\t%15.0f \t%15.0f", resultsAOsearch[scopeSize][PRUNING_NO][5] / COUNT, resultsSecondTree[scopeSize][PRUNING_NO][5] / COUNT);
//		if(Algorithms & AO_CONSTRAINTS)	fprintf(fpsum, "\n AND/OR CONSTRAINTS ONLY \t%15.0f", resultsAOsearch[scopeSize][PRUNING_CONSTRAINTS_ONLY][5] / COUNT);
//		if(Algorithms & AO_FC)			fprintf(fpsum, "\n AND/OR FC \t\t\t%15.0f", resultsAOsearch[scopeSize][PRUNING_FORWARD_CHECKING][5] / COUNT);
//		if(Algorithms & AO_FC_PROJ)		fprintf(fpsum, "\n AND/OR RFC \t\t\t%15.0f", resultsAOsearch[scopeSize][PRUNING_FC_PROJECTION][5] / COUNT);
//		if(Algorithms & AO_AUX)			fprintf(fpsum, "\n AUXILIARY AND/OR \t\t%15.0f", resultsAOsearch[scopeSize][AUXILIARY][5] / COUNT);
//		if(Algorithms & OR_NO_PRUNING)	fprintf(fpsum, "\n OR NO PRUNING \t\t\t%15.0f", resultsAOsearch[scopeSize][OR_PRUNING_NO][5] / COUNT);
//		if(Algorithms & OR_CONSTRAINTS)	fprintf(fpsum, "\n OR CONSTRAINTS \t\t%15.0f", resultsAOsearch[scopeSize][OR_PRUNING_CONSTRAINTS_ONLY][5] / COUNT);
//		if(Algorithms & OR_FC)			fprintf(fpsum, "\n OR FC \t\t\t\t%15.0f", resultsAOsearch[scopeSize][OR_PRUNING_FORWARD_CHECKING][5] / COUNT);
//		if(Algorithms & OR_FC_PROJ)		fprintf(fpsum, "\n OR RFC \t\t\t%15.0f", resultsAOsearch[scopeSize][OR_PRUNING_FC_PROJECTION][5] / COUNT);
//		if(Algorithms & AO_LINEAR_ADVANCED)	fprintf(fpsum, "\n AO LINEAR - ADVANCED  \t\t%15.0f \t%15.0f",	resultsAOsearch[scopeSize][AO_LINEAR_ADVANCED_ID][5] / COUNT, resultsSecondTree[scopeSize][AO_LINEAR_ADVANCED_ID][5] / COUNT);
//		if(Algorithms & AO_LINEAR_BE)	fprintf(fpsum, "\n AO LINEAR - BE \t\t%15.0f \t%15.0f",	resultsAOsearch[scopeSize][AO_LINEAR_BE_ID][5] / COUNT, resultsSecondTree[scopeSize][AO_LINEAR_BE_ID][5] / COUNT);
//		if(Algorithms & AO_ADVANCED_ADVANCED)		fprintf(fpsum, "\n AO ADVANCED - ADVANCED \t%15.0f \t%15.0f", resultsAOsearch[scopeSize][AO_ADVANCED_ADVANCED_ID][5] / COUNT, resultsSecondTree[scopeSize][AO_ADVANCED_ADVANCED_ID][5] / COUNT);
//		if(Algorithms & AO_ADVANCED_BE) fprintf(fpsum, "\n AO ADVANCED - BE \t\t%15.0f \t%15.0f",			resultsAOsearch[scopeSize][AO_ADVANCED_BE_ID][5] / COUNT, resultsSecondTree[scopeSize][AO_ADVANCED_BE_ID][5] / COUNT);
//
//	}
//
//
//	// Average number of zerodeadends
//	fprintf(fpsum, "\n\n\n=========================================================================================");
//	fprintf(fpsum, "\n\n\n\t\t\t\tAverage number of zero deadends:");
//	if((Algorithms & BE_AUX) && (BE_Count > 0))
//		fprintf(fpsum, "\n AUXILIARY BE \t\t\t%15.0f", resultsAOsearch[0][BE_AUXILIARY][6] / BE_Count);
//	else fprintf(fpsum, "\n AUXILIARY BE: w* too big");
////	for (scopeSize = iStart; scopeSize <= iEnd; scopeSize += ((iStep < (iEnd-scopeSize)) || (scopeSize == iEnd)) ? iStep : (iEnd-scopeSize))
//	for (scopeSize = iEnd; scopeSize >= iStart; scopeSize -= ((iStep < (scopeSize-iStart)) || (scopeSize == iStart)) ? iStep : (scopeSize-iStart) )
//	//scopeSize = (Width/2);
//	//while(scopeSize <= Width+1)
//	{
//	//	scopeSize += iStep ;
//		fprintf(fpsum, "\n\n--- (i = %d) ----------", scopeSize);
//		if(Algorithms & AO_NO_PRUNING)	fprintf(fpsum, "\n AO(i) - OLD\t\t\t\t%15.0f \t%15.0f", resultsAOsearch[scopeSize][PRUNING_NO][6] / COUNT, resultsSecondTree[scopeSize][PRUNING_NO][6] / COUNT);
//		if(Algorithms & AO_CONSTRAINTS)	fprintf(fpsum, "\n AND/OR CONSTRAINTS ONLY \t%15.0f", resultsAOsearch[scopeSize][PRUNING_CONSTRAINTS_ONLY][6] / COUNT);
//		if(Algorithms & AO_FC)			fprintf(fpsum, "\n AND/OR FC \t\t\t%15.0f", resultsAOsearch[scopeSize][PRUNING_FORWARD_CHECKING][6] / COUNT);
//		if(Algorithms & AO_FC_PROJ)		fprintf(fpsum, "\n AND/OR RFC \t\t\t%15.0f", resultsAOsearch[scopeSize][PRUNING_FC_PROJECTION][6] / COUNT);
//		if(Algorithms & AO_AUX)			fprintf(fpsum, "\n AUXILIARY AND/OR \t\t%15.0f", resultsAOsearch[scopeSize][AUXILIARY][6] / COUNT);
//		if(Algorithms & OR_NO_PRUNING)	fprintf(fpsum, "\n OR NO PRUNING \t\t\t%15.0f", resultsAOsearch[scopeSize][OR_PRUNING_NO][6] / COUNT);
//		if(Algorithms & OR_CONSTRAINTS)	fprintf(fpsum, "\n OR CONSTRAINTS \t\t%15.0f", resultsAOsearch[scopeSize][OR_PRUNING_CONSTRAINTS_ONLY][6] / COUNT);
//		if(Algorithms & OR_FC)			fprintf(fpsum, "\n OR FC \t\t\t\t%15.0f", resultsAOsearch[scopeSize][OR_PRUNING_FORWARD_CHECKING][6] / COUNT);
//		if(Algorithms & OR_FC_PROJ)		fprintf(fpsum, "\n OR RFC \t\t\t%15.0f", resultsAOsearch[scopeSize][OR_PRUNING_FC_PROJECTION][6] / COUNT);
//		if(Algorithms & AO_LINEAR_ADVANCED)	fprintf(fpsum, "\n AO LINEAR - ADVANCED  \t\t%15.0f \t%15.0f",	resultsAOsearch[scopeSize][AO_LINEAR_ADVANCED_ID][6] / COUNT, resultsSecondTree[scopeSize][AO_LINEAR_ADVANCED_ID][6] / COUNT);
//		if(Algorithms & AO_LINEAR_BE)	fprintf(fpsum, "\n AO LINEAR - BE \t\t%15.0f \t%15.0f",	resultsAOsearch[scopeSize][AO_LINEAR_BE_ID][6] / COUNT, resultsSecondTree[scopeSize][AO_LINEAR_BE_ID][6] / COUNT);
//		if(Algorithms & AO_ADVANCED_ADVANCED)		fprintf(fpsum, "\n AO ADVANCED - ADVANCED \t%15.0f \t%15.0f", resultsAOsearch[scopeSize][AO_ADVANCED_ADVANCED_ID][6] / COUNT, resultsSecondTree[scopeSize][AO_ADVANCED_ADVANCED_ID][6] / COUNT);
//		if(Algorithms & AO_ADVANCED_BE) fprintf(fpsum, "\n AO ADVANCED - BE \t\t%15.0f \t%15.0f",			resultsAOsearch[scopeSize][AO_ADVANCED_BE_ID][6] / COUNT, resultsSecondTree[scopeSize][AO_ADVANCED_BE_ID][6] / COUNT);
//	}
//
//	// Average probability of query
//	fprintf(fpsum, "\n\n\n=========================================================================================");
//	fprintf(fpsum, "\n\n\n\t\t\t\tAverage probability of query:");
//	if((Algorithms & BE_AUX) && (BE_Count > 0))
//		fprintf(fpsum, "\n AUXILIARY BE \t\t\t%15e", resultsAOsearch[0][BE_AUXILIARY][4] / BE_Count);
//	else fprintf(fpsum, "\n AUXILIARY BE: w* too big");
////	for (scopeSize = iStart; scopeSize <= iEnd; scopeSize += ((iStep < (iEnd-scopeSize)) || (scopeSize == iEnd)) ? iStep : (iEnd-scopeSize))
//	for (scopeSize = iEnd; scopeSize >= iStart; scopeSize -= ((iStep < (scopeSize-iStart)) || (scopeSize == iStart)) ? iStep : (scopeSize-iStart) )
//	//scopeSize = (Width/2);
//	//while(scopeSize <= Width+1)
//	{
//	//	scopeSize += iStep ;
//		fprintf(fpsum, "\n\n--- (i = %d) ----------", scopeSize);
//		if(Algorithms & AO_NO_PRUNING)	fprintf(fpsum, "\n AO(i) - OLD\t\t\t\t%15.0f \t%15.0f", resultsAOsearch[scopeSize][PRUNING_NO][4] / COUNT, resultsSecondTree[scopeSize][PRUNING_NO][4] / COUNT);
//		if(Algorithms & AO_CONSTRAINTS)	fprintf(fpsum, "\n AND/OR CONSTRAINTS ONLY \t%15e", resultsAOsearch[scopeSize][PRUNING_CONSTRAINTS_ONLY][4] / COUNT);
//		if(Algorithms & AO_FC)			fprintf(fpsum, "\n AND/OR FC \t\t\t%15e", resultsAOsearch[scopeSize][PRUNING_FORWARD_CHECKING][4] / COUNT);
//		if(Algorithms & AO_FC_PROJ)		fprintf(fpsum, "\n AND/OR RFC \t\t\t%15e", resultsAOsearch[scopeSize][PRUNING_FC_PROJECTION][4] / COUNT);
//		if(Algorithms & AO_AUX)			fprintf(fpsum, "\n AUXILIARY AND/OR \t\t%15e", resultsAOsearch[scopeSize][AUXILIARY][4] / COUNT);
//		if(Algorithms & OR_NO_PRUNING)	fprintf(fpsum, "\n OR NO PRUNING \t\t\t%15e", resultsAOsearch[scopeSize][OR_PRUNING_NO][4] / COUNT);
//		if(Algorithms & OR_CONSTRAINTS)	fprintf(fpsum, "\n OR CONSTRAINTS \t\t%15e", resultsAOsearch[scopeSize][OR_PRUNING_CONSTRAINTS_ONLY][4] / COUNT);
//		if(Algorithms & OR_FC)			fprintf(fpsum, "\n OR FC \t\t\t\t%15e", resultsAOsearch[scopeSize][OR_PRUNING_FORWARD_CHECKING][4] / COUNT);
//		if(Algorithms & OR_FC_PROJ)		fprintf(fpsum, "\n OR RFC \t\t\t%15e", resultsAOsearch[scopeSize][OR_PRUNING_FC_PROJECTION][4] / COUNT);
//		if(Algorithms & AO_LINEAR_ADVANCED)	fprintf(fpsum, "\n AO LINEAR - ADVANCED  \t\t%15.0f \t%15.0f",	resultsAOsearch[scopeSize][AO_LINEAR_ADVANCED_ID][4] / COUNT, resultsSecondTree[scopeSize][AO_LINEAR_ADVANCED_ID][4] / COUNT);
//		if(Algorithms & AO_LINEAR_BE)	fprintf(fpsum, "\n AO LINEAR - BE \t\t%15.0f \t%15.0f",	resultsAOsearch[scopeSize][AO_LINEAR_BE_ID][4] / COUNT, resultsSecondTree[scopeSize][AO_LINEAR_BE_ID][4] / COUNT);
//		if(Algorithms & AO_ADVANCED_ADVANCED)		fprintf(fpsum, "\n AO ADVANCED - ADVANCED \t%15.0f \t%15.0f", resultsAOsearch[scopeSize][AO_ADVANCED_ADVANCED_ID][4] / COUNT, resultsSecondTree[scopeSize][AO_ADVANCED_ADVANCED_ID][4] / COUNT);
//		if(Algorithms & AO_ADVANCED_BE) fprintf(fpsum, "\n AO ADVANCED - BE \t\t%15.0f \t%15.0f",			resultsAOsearch[scopeSize][AO_ADVANCED_BE_ID][4] / COUNT, resultsSecondTree[scopeSize][AO_ADVANCED_BE_ID][4] / COUNT);
//	}
//
//	fclose(fpsum);
//
//	return 1;
//}
//// end run instances Robert








// test function for and-or cutset
//int runTestAndOrCutset(int argc, char* argv[])
//{
//	CACHE_AT_AND_NODES		= true;
//	CACHE_AT_OR_NODES		= true;
//	CUTSET_CACHE			= false;
//	BRUTE_FORCE_W_CUTSET	= false; // to explore the cutset with linear space
//	SWITCH_TO_BE			= false; // this is to switch to bucket elimination when space permits
//
//	int		N = 5;					//number of variables
//	int		C = 3;					//number of CPTs
//	int		R = 2;					//number of root nodes
//	int		P_bayes = 2;			//number of parents for bayes
//	int		K = 2;					//number of values per variable
//
//	int		N_Const = 0;			//number of constraints
//	int		P_Const = 0;			//number of variables per constraint
//	double	tightness = 0.10;		//tightness of constraints
//	int		Algorithms = 4;
//
//	int iStart =0;
//	int iEnd = 30;
//	int iStep = 1;
//
//	int COUNT = 10;
//
//	char inputFile[256];
//	strcpy(inputFile, "");
//	//=======================================================================
//	double avgWidth = 0;
//	double avgHeight = 0;
//	double totalAvgWidth = 0;
//	double totalAvgHeight = 0;
//
//	int BE_Count = 0;
//
//	if (argc <= 1)
//		return 0;
//
//	// begin initialization
//	initExperiment(N, C,  R,  P_bayes,  K,  N_Const,  P_Const,  tightness,
//					 Algorithms,  iStart,  iEnd,  iStep,  COUNT, argc, argv, inputFile);
//
//
//
//	bool probCreatedFlag = false;
//	CProblemMixed* prob;
//
////*1111111111111111111111111111111
//	// create random problem
//	while (!probCreatedFlag)
//	{
//		prob = new CProblemMixed(GRAPH_RANDOM_FIXED, N, K, C, N_Const, P_Const, tightness);
//		prob->setParamParents(P_bayes);
//		prob->setParamConnectivity(UNKNOWN);
//		prob->setParamWidth(UNKNOWN);
//		prob->setParamHeight(UNKNOWN);
//		prob->setParamSigma(UNKNOWN);
//		prob->create();
//
////		prob->extractDeterminism();
//
//		// create problem
//		if (prob->preprocess2(VO_MINFILL, true))
//			probCreatedFlag = true;
//		else
//			if (prob) delete prob;
//
//
//	}
//
////11111111111111111111111111111111
////*/
//
//
////---------------------------------------
//
///*2222222222222222222222222222222222
//	// read CPCS file
//	prob = new CProblemMixed();
//
//	// CPCS
////	prob->load(INP_SIMPLE, "C:\\Robert\\Research\\Benchmarks\\CPCS\\cpcs54.simple");
//
//	// GENETIC LINKAGE
////	prob->load2(INP_SIMPLE, "C:\\Robert\\Research\\Benchmarks\\GENETIC_LINKAGE\\BayesXml\\EA\\simple\\comp01_fileEA4.simple");
//
//
////	prob->extractDeterminism();
//
//	if (prob->preprocess(VO_MINFILL, true))
//		probCreatedFlag = true;
//	else
//	{
//		if (prob) delete prob;
//	}
////222222222222222222222222222222222222
//*/
//
//
//
//
//
////=================================================================================================//
///*
//
//
//	// create a tree decomposition
//	prob->makeTreeDecomposition(VO_MINFILL);
//	CTreeDecomposition* treeDecomp = prob->getTreeDecomposition();
//	treeDecomp->computeSeparatorWeights();
//
//	prob->getTreeDecomposition()->print();
//
//	prob->makeAndOrCutsetTree();
//
//	CLegalTree* w_cutsetTree = prob->getTreeCutset();
//	int* w_values = w_cutsetTree->getWCutsetHeight();
//
//	for (int i = prob->getN(); i>0 ; i--)
//		cout << "\n w=" << i << " - cutset has depth: " << w_values[i-1] + 1 << " ---- f(" << i << ") = " << w_values[i-1] + i + 1 ;
//
////	prob->getTreeCutset()->print();
////	prob->getTreeDecomposition()->print();
//
//
//
//	// try some belief updating
//
//	prob->preprocessCutsetTree();
//
//
////	prob->print();
//
////=================================================================================================//
//*/
//
//	double tm;
//	long exps = 0;
//	long andN = 0, orN = 0;
//
//	int scope = 10;
//
////	fprintf(fp, "\n\n AND/OR FC");
//	((CProblemMixed*) prob)->setAdvancedCacheFlag(scope);
//
//	// for Brute Force W cutset;
//	((CProblemMixed*) prob)->wCutsetInit();
//
//	prob->cacheInit(scope);
////	prob->print();
//
//	prob->execAO(-1, tm, exps, andN, orN, false, PRUNING_NO, scope);
//
//	prob->cacheDestroy();
//
//	if(prob) delete prob;
//
//	return 1;
//}
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//



void initCompetitionUAI(int& elimType, int& iBound, int& usePruning,  bool& bucketElimination, double& expMax, char inputFile[], int argc, char* argv[])
{
	int i ;
	for (i = 1 ; i < argc ; i++)
	{
		char *arg = argv[i] ;
		//printf("\n[%s]", arg) ;

		// read in switches
		if (arg[0] == '-')
		{
			switch (arg[1])
			{
				case 'i':
				case 'I':
					// format is -iBound=15
					{
						char *iPar ;
						iPar = strstr(arg, "Bound") ;
						if (iPar)
						{
							iPar += 6 ;
							iBound = atoi(iPar) ;
							break ;
						}

					}
				case 'e':
				case 'E':
					// format is -elim=0
					{
						char *iPar ;
						iPar = strstr(arg, "elim") ;
						if (iPar)
						{
							iPar += 5 ;
							elimType = atoi(iPar) ;
							break ;
						}

					}

				case 'x':
				case 'X':
					// format is -xp=700
					{
						char *iPar ;
						iPar = strstr(arg, "xp") ;
						if (iPar)
						{
							iPar += 3 ;
							expMax = atof(iPar) ;
							break ;
						}

					}

				case 'f':
				case 'F':
					// format is -file=file_name
					{
						char *iPar;
						iPar = strstr(arg, "file") ;
						if (iPar)
						{
							iPar += 5 ;
							strcpy(inputFile, iPar);
							break ;
						}
					}

				case 'u':
				case 'U':
					// format is -usePruning=3
					{
						char *iPar;
						iPar = strstr(arg, "usePruning") ;
						if (iPar)
						{
							iPar += 11 ;
							usePruning = atoi(iPar);
							break ;
						}
					}

				case 'b':
				case 'B':
					// format is -be=0 -be=1
					{
						char *iPar;
						iPar = strstr(arg, "be") ;
						if (iPar)
						{
							iPar += 3 ;
							int temp = atoi(iPar);
							if(temp == 0)
								bucketElimination = false;
							else
								bucketElimination = true;
							break ;
						}
					}


			}
		}
	}
}





int runCompetitionUAI06(int argc, char* argv[])
{
	CACHE_AT_AND_NODES		= false;
	CACHE_AT_OR_NODES		= true;
	CUTSET_CACHE			= true;
	SWITCH_TO_BE			= false;

	double logProb = 0.0;

	int iBound = 10;
	char inputFile[256];
	char outputFile[256];
	strcpy(outputFile, "longResults.txt");
	char networkName[256];
	int usePruning = 0;
	int elimType = VO_MINFILL;
	double expMax = 700.0;

	initCompetitionUAI(elimType, iBound, usePruning, SWITCH_TO_BE, expMax, inputFile,argc,argv);
	cout << "iBound = " << iBound << endl;
	cout << "usePruning = " << usePruning << endl;
	if (SWITCH_TO_BE)
		cout << "SWITCH TO BE = true" << endl;
	else
		cout << "SWITCH TO BE = false" << endl;
	cout << "file = " << inputFile << endl;



	CProblemMixed* prob = new CProblemMixed();
	double totalTime_start, totalTime_end;
	double readTime = 0.0;
	double preprocessTime = 0.0;
	double searchTime = 0.0;
	double t_start, t_end;

	t_start = cpuTime();
	totalTime_start = t_start;

	cout << " Starting to read file ... " << endl;
	prob->parser_simple2(inputFile);
	cout << " Done reading file. " << endl;
	t_end = cpuTime();
	readTime += t_end - t_start;

	FILE* fp = fopen(outputFile, "a");
	if (NULL == fp) return 0;

	char *iPar;
	iPar = strstr(inputFile, "BN") ;
	strcpy(networkName, iPar);
    char * pch;
    pch = strstr (networkName,".simple");
    strncpy (pch,"",7);

	fprintf(fp, "\n\n ====================================================================================================");
	fprintf(fp, "\n\n Network name: %s", networkName);
	fprintf(fp, "\n N = %d", prob->getN());
	fprintf(fp, "\n max K = %d", prob->getK());

	//prob->findDeterministicVariables();
	//char graphViz[256];
	//strcpy(graphViz, networkName);
	//strcat(graphViz, "Viz.txt");
	//FILE* fff = fopen(graphViz, "w");
	//fprintf (fff, "digraph G {");
	//fprintf(fff, "\n subgraph G1 { \n node [style=filled color=red];");
	//for (int i=0; i < prob->getN(); i++)
	//	if (prob->isDeterministic(i))
	//		fprintf(fff, "\n %d;",i);
	//fprintf(fff, "\n }");

	//function_map::iterator itf = prob->functions().begin();
	//for( ; itf != prob->functions().end() ; ++itf )
	//{
	//	CFunction* f = (*itf).second;
	//	int argc = f->getArgc();
	//	int* argv = f->getArgv();
	//	for (int i = 0; i < argc -1 ; i++ )
	//	{
	//		fprintf(fff, "\n   %d -> %d;", argv[i], argv[argc-1]);
	//	}
	//}
	//fprintf (fff, "\n}");
	//fclose(fff);

//	return 1;

	t_start = cpuTime();
		prob->markBarren_removeEvidence_makeConnectedComponents();
		logProb += log(globalCostFromEvidence);
	t_end = cpuTime();
	preprocessTime += t_end - t_start;


	fprintf(fp, "\n Relevant N = %d", prob->getN());
	fprintf(fp, "\n Number of components = %d", independentComponentProblems.size());
	fflush(fp);

	// compute scaling factor
	double maxN = 0.0;
	vector<CProblemMixed*>::iterator it1 = independentComponentProblems.begin();
	for ( ; it1 != independentComponentProblems.end() ; ++it1)
	{
		CProblem* prob = (*it1);
		int n = prob->getN();
		if ( (double) n > maxN )
			maxN = (double) n;
	}

	double maxScalingFactor;
	if(expMax / maxN > 700)
		maxScalingFactor = exp(700.0);
	else
		maxScalingFactor = exp( expMax / maxN );

	vector<CProblemMixed*>::iterator it = independentComponentProblems.begin();
	int count = 0;
	for ( ; it != independentComponentProblems.end() ; ++it)
	{
		prob = (*it);
		count ++;

		scalingFactor = exp ( log(maxScalingFactor) * ( (double) maxN / (double)prob->getN() ) );
		cout << "\n scaling factor = " << scalingFactor;

		fprintf(fp, "\n ---------------------------------");
		fprintf(fp, "\n Component %d (of %s)", count, networkName);
		fprintf(fp, "\n N = %d", prob->getN());
		fprintf(fp, "\n max K = %d", prob->getK());

		t_start = cpuTime();
			bool connected = prob->preprocessAdaptiveCaching(elimType, true);
			cout << " Induced width: " << prob->getInducedWidth() << endl;
			fprintf(fp, "\n w* = %d", prob->getInducedWidth());

			//continue;
			if(usePruning != 0)
				prob->extractDeterminism();

			//prob->createShannonTrees();

			//	fprintf(fp, "\n\n AND/OR FC");
			((CProblemMixed*) prob)->setAdvancedCacheFlag(iBound);

			// for Brute Force W cutset - needed for SWITCH_TO_BE also;
			// initialize the OR adaptive contexts
			((CProblemMixed*) prob)->wCutsetInit();

			prob->cacheInit(iBound);
			prob->initConstraintPropagation();
		t_end = cpuTime();
		preprocessTime += t_end - t_start;
		fprintf(fp, "\n Preprocess time: %f", t_end - t_start );

		double tm;
		long exps = 0;
		long andN = 0, orN = 0;


		///////////////////////////////////////////////////////////////////////////////////////////
		//FILE* fff;
		//prob->findDeterministicVariables();

		//char graphViz[256];
		//strcpy(graphViz, networkName);
		//strcat(graphViz,"_comp");
		//char buffer[256];
		//strcat(graphViz,itoa(count,buffer,10));
		//strcat(graphViz, "Viz.txt");
		//fff = fopen(graphViz, "w");
		//fprintf(fff, "digraph G {");
		//fprintf(fff, "\n subgraph G1 { \n node [style=filled color=red];");
		//for (int i=0; i < prob->getN(); i++)
		//	if (prob->isDeterministic(i))
		//		fprintf(fff, "\n %d;",i);
		//fprintf(fff, "\n }");
		//function_map::iterator itf = prob->functions().begin();
		//for( ; itf != prob->functions().end() ; ++itf )
		//{
		//	CFunction* f = (*itf).second;
		//	int argc = f->getArgc();
		//	int* argv = f->getArgv();
		//	for (int i = 0; i < argc -1 ; i++ )
		//	{
		//		fprintf(fff, "\n   %d -> %d;", argv[i], argv[argc-1]);
		//	}
		//}
		//fprintf (fff, "\n}");
		//fclose(fff);
		/////////////////////////////////////////////////////////////////////////////////////////


		percentageDoneVariable = prob->findVariableToComputePercentageDone(topChainVariables);
		t_start = cpuTime();
			prob->execAO(-1, tm, exps, andN, orN, false, usePruning, iBound);
		t_end = cpuTime();
		searchTime += t_end - t_start;
		fprintf(fp, "\n Solving time: %f", t_end - t_start );

		t_start = cpuTime();
			// cut out the scaling factor
			logProb += log(prob->getProbabilityOfQuery()) - expMax;
			prob->cacheDestroy();
			delete prob;
		t_end = cpuTime();
		preprocessTime += t_end - t_start;
	}
	independentComponentProblems.clear();
	topChainVariables.clear();

	totalTime_end = cpuTime();

	cout << endl;
	cout << "Log probability of evidence: " << logProb << endl;

	cout << "TIME: " << endl;
	cout << "   - Read file time:       " << readTime  << endl;
	cout << "   - Preprocessing time:   " << preprocessTime  << endl;
	cout << "   - Solving time:         " << searchTime  << endl;
	cout << "                         --------------------------------" << endl;
	cout << "   - TOTAL                 " << totalTime_end - totalTime_start << endl;

	fprintf(fp, "\n ---------------------------------" );
	fprintf(fp, "\n Final results (%s):", networkName);
	fprintf(fp, "\n Log probability of evidence: %f", logProb );
	fprintf(fp, "\n TIME:" );
	fprintf(fp, "\n    - Read file time:       %f", readTime );
	fprintf(fp, "\n    - Preprocessing time:   %f", preprocessTime );
	fprintf(fp, "\n    - Solving time:         %f", searchTime );
	fprintf(fp, "\n                         --------------------------------");
	fprintf(fp, "\n    - TOTAL                 %f", totalTime_end - totalTime_start );
	fprintf(fp, "\n    - ALL BUT READ IN       %f", totalTime_end - (totalTime_start + readTime) );
	fclose(fp);

	char outputFileShort[256];
	strcpy(outputFileShort, "shortResults.txt");

	fp = fopen(outputFileShort, "a");
	if (NULL == fp) return 0;

	fprintf(fp, "\n%s", networkName );
	fprintf(fp, "\t%f", logProb );
	fprintf(fp, "\t%f", totalTime_end - totalTime_start  );
	fprintf(fp, "\t%f", totalTime_end - (totalTime_start + readTime)  );

	fprintf(fp, "\t-elim=%d", elimType);
	fprintf(fp, " -iBound=%d", iBound);
	fprintf(fp, " -usePruning=%d", usePruning);
	fprintf(fp, " -be=%d", SWITCH_TO_BE);
	fprintf(fp, " -xp=%f", expMax);
	fprintf(fp, " -file=\"%s\"", inputFile);
	fclose(fp);

	return 1;
}




int main(int argc, char *argv[])
{
	//	int x = runExperiment(argc, argv);

//	runTestAndOrCutset(argc, argv);

	int x = runCompetitionUAI06(argc,argv);
//		getchar();

	return 1;
}
