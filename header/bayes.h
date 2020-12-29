// bayes.h -- Bayesian Networks.

/*
 * Copyright (C) 2004 Radu Marinescu
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

/*
 * NOTE: This is an interal header file.
 * You should not attempt to use it directly.
 */

#pragma warning (disable : 4786)


#ifndef _MPELIB_BAYES_H
#define _MPELIB_BAYES_H

#include "ProblemMixed.h"
#include "problem.h"
#include "graph.h"
#include "legaltree.h"
#include "bucketstruct.h"
#include "aonode.h"

extern long deadEnds;   // to count the deadends due to the constraints
extern long zeroDeadEnds; // to count the deadends due to 0's in CPT's


////////////////////////////////////////////////////////////////////////
// CProblemBayes class implementation.

class CProblemBayes : public CProblem
{
//protected:
//	bool m_connected;				// Flag indicating if underlying graph is connected.
//	int  m_type;					// Bayes Network type (random, random fixed, grid, coding)
//	
//	CGraph* m_graph;				// Underlying moral graph.
//	CLegalTree* m_tree;				// Rooted tree arrangement.
//	CBucketStruct* m_buckets;		// Bucket structure associated with the problem instance.
//
//	int m_paramConnectivity;		// Random network controlled by connectivity [Darwiche99].
//	int m_paramParents;				// Random network with fixed number of parents.
//	int m_paramHeight;				// Grid network - height.
//	int m_paramWidth;				// Grid network - width.
//	double m_paramSigma;			// Coding network - channel noise.
//
//	FILE* m_outputFile;
//
//	double m_probabilityOfQuery;	// The probability of the constraint part in mixed networks
//	double m_bestSolutionCost;		// Best solution cost.
//	int* m_bestSolution;			// Best solution assignement.
//
//	int* m_orderingDFS;				// Ordering of variables (DFS ordering of legal tree).
//	int* m_positionDFS;				// Position of variables (DFS ordering of legal tree).
//
//	// Search tree level.
//	SLINFO* m_levelInfo;			// Domain manager.
//	variable_v* m_descendants;
//
//protected:
//	bool init();					// Initialize the network.
//public:
//	bool initialize(function_map& functions_bayes, function_map& functions_constraint, CProblemMixed* prob, int algType);				// Initialize an auxiliary network from a mixed one
//	CLegalTree* getTree() { return m_tree; };
//	void setTree(CLegalTree* tree) { m_tree = tree; };
//	CGraph* getGraph() { return m_graph; };
//	double getProbabilityOfQuery() { return m_probabilityOfQuery; };
//
//
//public:
//	bool create();					// Create a network instance.
//	bool preprocess(int voType, bool withTree);
//	void print();
//
//	void createCustom();
//
//	int getTreeHeight();
//	int getInducedWidth();
//	variable_v& getTreeDescendants(int var);
//
//	virtual int* getOrdering(bool dfsOrd = false);
//	virtual int* getPosition(bool dfsOrd = false);
//
//	function_v& getBucket(int var, bool dfsOrd = true);
//	double getBestSolutionCost()		{ return m_bestSolutionCost; };
//
//public:
//	void setParamParents(int p)			{ m_paramParents = p; };
//	void setParamConnectivity(int p)	{ m_paramConnectivity = p; };
//	void setParamWidth(int p)			{ m_paramWidth = p; };
//	void setParamHeight(int p)			{ m_paramHeight = p; };
//	void setParamSigma(double p)		{ m_paramSigma = p; };
//
//	bool assertEvidence(int ne);
//
//protected:
//	int selectValue(int var, int& val, bool prune = true);		// Value selection (search).
//	int selectVariable(int svType);			// Variable selection (search).
//
//	bool initLevelInfo();					// Init domain list stack.
//	void destroyLevelInfo();				// Destroy domain list stack.
//	void pushLevelInfo();					// Push domain stack (increase size).
//	void popLevelInfo();					// Pop domain stack (decrease size).
//
//	void backupDomains();					// Backup domains.
//	void restoreDomains();					// Rollback domains.
//
//	void saveCurrentSolution();				// Save current solution as best solution.
//	double utility(int var = -1);			// Computes the utility of the current (partial) assignment.
//
//	void evalStaticHeuristic(int var);		// Evaluate pre-computed heuristic for current variable.
//	void evalDynamicHeuristic(int var, int ibound);
//
//	bool precompute(int ibound, double& cpuPreprocess);
//
//	// AND-OR search tree node expansion.
//	int expand(CAONode* node, stack<CAONode*>& succ, bool heuristic, int ibound);
//	bool forwardPrune(CAONode* node, double estimate);
//	void removeNode(list<CAONode*>& l, CAONode* node);
//public:
//	void test();
//
//
//	// Branch and Bound with static (precomputed) MB heuristics.
//	int execBBMBs(int ibound, long tmLimit, double& tmCpuSearch,
//		long& backtracks, long& expansions, bool silent = false);
//
//	// Branch and Bound with dynamic MBE heuristics.
//	int execBBMBd(int ibound, long tmLimit, double& tmCpuSearch,
//		long& backtracks, long& expansions, bool silent = false);
//
//	// Backtracking search (full tree expansion).
//	int execBT(long tmLimit, double& tmCpuSearch, long& backtracks, 
//		long& expansions, bool silent = false);
//
//	// AND-OR search (full tree expansion).
//	int execAO(long tmLimit, double& tmCpuSearch, long& expansions, 
//		long& nodesAND, long& nodesOR, bool silent = false);
//
//	// AND-OR search with MBE pruning (Alpha-Beta pruning style of the tree).
//	int execAOMB(int ibound, long tmLimit, double& tmCpuSearch, 
//		long& expansions, long& nodesAND, long& nodesOR, bool silent = false);
//
//	// Bucket Elimination on the auxiliary network
//	int execBE_AUX(double& timeCpuBe, double& result);
//	
//public:
//	CProblemBayes();
//	CProblemBayes(int t, int n, int k, int p);
//
//	virtual ~CProblemBayes();
};

#endif	// _MPELIB_BAYES_H

// Local Variables:
// mode: C++
// End:
