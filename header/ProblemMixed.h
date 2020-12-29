// ProblemMixed.h: interface for the CProblemMixed class -- Mixed Networks.

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

/*
 * NOTE: This is an internal header file.
 * You should not attempt to use it directly.
 */

#pragma warning (disable : 4786)

#ifndef _MPELIB_MIXED_H
#define _MPELIB_MIXED_H

#include "problem.h"
#include "graph.h"
#include "legaltree.h"
#include "bucketstruct.h"
#include "aonode.h"
#include "constraintnode.h"
#include "ConstraintTable.h"
#include "CacheTable.h"
#include "TreeDecomposition.h"
#include "bitvector.h"
#include <string>

extern long deadEnds;   // to count the deadends due to the constraints
extern long zeroDeadEnds; // to count the deadends due to 0's in CPT's


////////////////////////////////////////////////////////////////////////
// CProblemMixed class implementation.

class CProblemMixed : public CProblem
{
protected:
	
	int m_Constraints;				// Number of constraints (for the constraint part of the mixed net)

protected:
    // The constraints in this problem.  We store them as pointers
    // we want to be able to have variables and constraints of
    // different types. We also prefer hash tables to lists and
    // vectors since they allow O(1) search & removals w/o any
    // need for address readjustment.

	function_map m_functionsConstraint;		// Constraints of the mixed problem.

	CFunctionList** m_adjacenciesConstraint;	// The adjacency matrix for the constraint problem

protected:
	bool m_connected;				// Flag indicating if underlying moral mixed graph is connected.
    bool m_bayesConnected;			// Flag indicating if underlying moral bayes graph is connected.
	bool m_constraintConnected;				// Flag indicating if underlying constraint graph is connected.

	double m_tightness;				// tightness of the constraints (% of allowed tuples)

	CBitVector* m_domainsPropagation;  // Domains of variables for constraint propagation (0 = invalid, 1 = valid).


	int  m_type;					// Mixed Network type (random, random fixed, grid, coding)
	
	CGraph* m_graph;				// Underlying moral graph (mixed)
	CGraphHash m_graph2;
	CGraph* m_constraintGraph;		// Constraint graph (robert)
	CGraph* m_bayesGraph;			// Mixed graph (robert)
	CLegalTree* m_tree;				// Rooted tree arrangement.
	CBucketStruct* m_buckets;		// Bucket structure associated with the problem instance.

	int m_paramConnectivity;		// Random network controlled by connectivity [Darwiche99].
	int m_paramParents;				// Random network with fixed number of parents.
	int m_paramConstraintParents;	// Random mixed net with fixed number of variables per constraint
	int m_paramHeight;				// Grid network - height.
	int m_paramWidth;				// Grid network - width.
	double m_paramSigma;			// Coding network - channel noise.

	FILE* m_outputFile;

	double m_probabilityOfQuery;	// The probability of the constraint part in mixed networks
	double m_numberOfSolutions;		// The number of solutions of the constraint network
	double m_bestSolutionCost;		// Best solution cost.
	int* m_bestSolution;			// Best solution assignement.

	int* m_orderingDFS;				// Ordering of variables (DFS ordering of legal tree).
	int* m_positionDFS;				// Position of variables (DFS ordering of legal tree).

	int* m_ordering_backup;			// Backup Ordering of variables 
	int* m_position_backup;			// Backup Position of variables

	int* m_ordering_cutsetDFS;		// Ordering of variables for cutset tree (DFS ordering of m_tree_cutset).
	int* m_position_cutsetDFS;		// Position of variables for cutset tree (DFS ordering of m_tree_cutset).

	bool* m_isDeterministic;


	// Search tree level.
	SLINFO* m_levelInfo;			// Domain manager.
	variable_v* m_descendants;
	variable_v*	m_parentSepSet;		// parent separator set for caching in search (AND caching)
	variable_v* m_parentSet;		// parent set for caching (OR caching)
	variable_v* m_advancedCacheFlagAND; // for each variable it keeps a vector of variables for which to purge the cache
	variable_v* m_advancedCacheFlagOR; // for each variable it keeps a vector of variables for which to purge the cache

	CCacheTable** m_cacheFull;			// cache for search
	CCacheTable** m_cacheSep;			// cache for search
	CCacheTable** m_cacheFullCounting;	// cache for search counting
	CCacheTable** m_cacheSepCounting;	// cache for search counting

	CLegalTree* m_tree_cutset;
	CLegalTree* m_treeBackup;			// this does not need to be deleted

	int* m_wCutset;						// The max w for which (variable belongs to minimal height w-cutset)
	int* m_wCutsetHeight;				// Vector of heights of AO-w-cutsets; m_height_w_cutset[0] = m_height;

	CTreeDecomposition* m_treeDecomposition;	// a tree decomposition of the problem

	vector<int> m_superlinkOrder;
	vector<int> m_superlinkCutset;
	map<string, int> m_mapString2Int;
	map<int, string> m_mapInt2String;

	// for each variable, keep a vector of pointers to constraints that have it in scope
	vector<function_v> m_variablesInFunctions;

	int m_r; // max scope size

	bool** m_validValues;  // adjacency table for valid variable values (used in constraint propagation)
	int* m_validValuesCount; // number of valid values

public:
	function_map& functionsConstraint()		{ return m_functionsConstraint; };
	CLegalTree* getTree() { return m_tree; };
	void setTree(CLegalTree* tree) { m_tree = tree; };
	void setTreeCutset(CLegalTree* tree) { m_tree_cutset = tree; };

	bool isDescendant(int ancestor, int var);

	// Add a new function to the repository.
	virtual void addConstraintFunction(CConstraintTable* c, bool w_id = true);
	// Get a function from the repository.
	CConstraintTable* getConstraintFunctionByID(int id);
	double getProbabilityOfQuery() { return m_probabilityOfQuery; };
	double getNumberOfSolutions() { return m_numberOfSolutions; };
	void parentsetCreate();
	void parentsetDestroy();
	void parentsetInit();
	void parentsetInit(int* substitution);
	void cacheInit(int scopeBound);
	void cacheDestroy(); 
	void setAdvancedCacheFlag(int scope);
	void setHighestVarInScope();

	//tree decomposition functions
	void makeTreeDecomposition(int type);		// type of ordering
	void makeAndOrCutsetTree(int scopeSize);	// the legal tree will be annotated for w-cutset
	void makeGWCutsetTree(int scopeSize);		// Bozhena's GWC (regular cutset); but organized as AND/OR already
												// the legal tree will be annotated for w-cutset
	void makeCutsetOr(int scopeSize);			// traverse the scopeSize cutset and make it a chain

	void makeTreeDecompositionClusters(int type, treeDecompositionCluster_v &clusters);			// computes the clusters given the type of the ordering
	CTreeDecomposition* makeTreeDecompositionClustersAssemble(treeDecompositionCluster_v &clusters);			// assembles the clusters in a tree decomposition, computes separators

	CLegalTree* getTreeCutset()		{ return m_tree_cutset; };
	CTreeDecomposition* getTreeDecomposition() { return m_treeDecomposition; };



	bool load(int itype, char* filename);  // from radu, to read cpcs files
	bool load2(int itype, char* filename);  // from radu, to read genetic files
	bool parser_simple2(char *filename);

	void createGridProblem(int N);			// creates a grid problem (Shuki Bruck); no adjacent 1's, how many solutions?

	void extractDeterminism();				// extracts flat CPTs and adds them to m_functionsConstraint

	void createShannonTrees();				// creates Shannon Trees for each constraint; needs determined ordering

	void setWCutsetHeight();				// sets m_wCutsetHeight;

	void removeEvidence(bool* barrenVariables = NULL);

	void reset();
	void makeConnectedComponents();
	void createFromComponnet(int N, variable_v vars, function_map functions);
	void markBarrenVariables(bool*& barrenVariables);
	
	void markBarren_removeEvidence_makeConnectedComponents();
	int findVariableToComputePercentageDone(variable_v& vars);
	double updatePercentageDone();
	void createProblemFromComponent(variable_v& vars, function_map& functions, int* domains, int* visited, int componentNumber );
	bool findDeterministicVariables();
	void fixOrderingDeterministic(variable_s* parents);
	void fixOrderingDeterministic2(variable_s* parents, int* substitution);
	bool isDeterministic(int var)		{ return m_isDeterministic[var]; };  
	void findTopologicalOrdering(variable_v& order);
	void computeFractionDone(int currentVar);

protected:
	bool init();					// Initialize the network.

public:
	bool isEligibleToBeChecked(CFunction* fun); //Checks if the scope of the function is included in the current assignment
	bool consistent_old(int var, int val); // Checks the consistency of the current assignment ending in (var) for and/or search in mixed networks
	bool forwardChecking_old(int var, int val); // Checks the consistency of the current assignment ending in var by FC for and/or search in mixed networks
	bool forwardCheckingByProjection_old(int var, int val); // Checks the consistency by looking at all constraints, projecting over the current assignment
	bool forwardCheckingByProjection(int var, int val); // Checks the consistency by looking at all constraints, projecting over the current assignment

	bool consistent(int var, int val); // Checks the consistency of the current assignment ending in (var) for and/or search in mixed networks
	bool forwardChecking(int var, int val); // Checks the consistency of the current assignment ending in var by FC for and/or search in mixed networks
	bool forwardChecking(CAONode* node); // Checks the consistency of the current assignment ending in var by FC for and/or search in mixed networks
	bool forwardCheckingWithMarking(CAONode* node); // Checks the consistency of the current assignment ending in var by FC for and/or search in mixed networks
	bool supportIn(CFunction* fun, int var, CAONode* node); // Checks if there is support in fun for unassigned variable var

	bool initConstraintPropagation();  // Initialize constraint propagation data structures

	bool create();					// Create a network instance.
	bool preprocess(int voType, bool withTree);
	bool preprocess2(int voType, bool withTree);
	bool preprocessAdaptiveCaching(int voType, bool withTree);
	bool preprocessCutsetTree();

	void wCutsetInit();

	void print();

	void createCustom();

	int getNConstraints()				{ return m_Constraints; };
	int getTreeHeight();
	int getInducedWidth();
	variable_v& getTreeDescendants(int var);

	virtual int* getOrdering(bool dfsOrd = false);
	virtual int* getPosition(bool dfsOrd = false);

	void setOrdering(int* ordering);
	void setPosition(int* position);

	function_v& getBucket(int var, bool dfsOrd = true);
	double getBestSolutionCost()		{ return m_bestSolutionCost; };

	int* getWCutset()					{ return m_wCutset; };
	int* getWCutsetHeight()				{ return m_wCutsetHeight; };

public:
	void setParamParents(int p)			{ m_paramParents = p; };
	void setParamConnectivity(int p)	{ m_paramConnectivity = p; };
	void setParamWidth(int p)			{ m_paramWidth = p; };
	void setParamHeight(int p)			{ m_paramHeight = p; };
	void setParamSigma(double p)		{ m_paramSigma = p; };

	bool assertEvidence(int ne);

	void backupOrdering();
	void backupPosition();

	void restoreOrdering();
	void restorePosition();



protected:
	int selectValue(int var, int& val, bool prune = true);		// Value selection (search).
	int selectVariable(int svType);			// Variable selection (search).

	bool initLevelInfo();					// Init domain list stack.
	void destroyLevelInfo();				// Destroy domain list stack.
	void pushLevelInfo();					// Push domain stack (increase size).
	void popLevelInfo();					// Pop domain stack (decrease size).

	void backupDomains();					// Backup domains.
	void restoreDomains();					// Rollback domains.

	void saveCurrentSolution();				// Save current solution as best solution.
	double utility(int var = -1);			// Computes the utility of the current (partial) assignment.

	void evalStaticHeuristic(int var);		// Evaluate pre-computed heuristic for current variable.
	void evalDynamicHeuristic(int var, int ibound);

	bool precompute(int ibound, double& cpuPreprocess);

	// AND-OR search tree node expansion.
	int expand(CAONode* node, stack<CAONode*>& succ, bool heuristic, int ibound, int usePruning = 0);
	bool forwardPrune(CAONode* node, double estimate);
	void removeNode(list<CAONode*>& l, CAONode* node);

public:
	void test();


	// Branch and Bound with static (precomputed) MB heuristics.
	//int execBBMBs(int ibound, long tmLimit, double& tmCpuSearch,
	//	long& backtracks, long& expansions, bool silent = false);

	// Branch and Bound with dynamic MBE heuristics.
	//int execBBMBd(int ibound, long tmLimit, double& tmCpuSearch,
	//	long& backtracks, long& expansions, bool silent = false);

	// Backtracking search (full tree expansion).
	int execBT(long tmLimit, double& tmCpuSearch, long& backtracks, 
		long& expansions, bool silent = false);

	// AND-OR search (full tree expansion).
	
	int execAO(long tmLimit, double& tmCpuSearch, long& expansions, 
		long& nodesAND, long& nodesOR, bool silent = false, int usePruning = 0, int cacheSize = 0);

	// AND-OR search with MBE pruning (Alpha-Beta pruning style of the tree).
	//int execAOMB(int ibound, long tmLimit, double& tmCpuSearch, 
	//	long& expansions, long& nodesAND, long& nodesOR, bool silent = false);
	
public:
	CProblemMixed();
	CProblemMixed(int t, int n, int k, int p, int nConstraints, int sizeConstraint, double tightness);
	CProblemMixed(int n, int k);

	virtual ~CProblemMixed();
};

#endif	// _MPELIB_MIXED_H

// Local Variables:
// mode: C++
// End:
