// ProblemMixed.cpp -- Mixed Networks.

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




#include "ProblemMixed.h"
#include "cpt.h"
#include "function.h"
#include "CacheTable.h"
#include "defs.h"
#include <memory.h>
#include "TreeDecomposition.h"
#include "TreeDecompositionCluster.h"
#include "TreeDecompositionSeparator.h"
#include "timers.h"



double c_expand_AND_start;
double c_expand_AND_end;
double c_expand_OR_start;
double c_expand_OR_end;
double tmExpandAND = 0.0;
double tmExpandOR = 0.0;

extern bool CACHE_AT_AND_NODES;
extern bool CACHE_AT_OR_NODES;
extern bool CUTSET_CACHE;
extern bool COUNTING;
extern bool BRUTE_FORCE_W_CUTSET;
extern bool SWITCH_TO_BE;
extern bool g_useSuperlinkOrder;
extern vector<CProblemMixed*>	independentComponentProblems; // to keep the connected components of a problem
extern variable_v topChainVariables;
extern int percentageDoneVariable;
double percentageDone = 0.0;

double globalCostFromEvidence = 1.0;
extern double scalingFactor;

/////////////////////////////////////////////////////////////////////////
// CProblemMixed class implementation.
CProblemMixed::CProblemMixed() : CProblem()
{
	m_N = 0;
	m_Constraints = 0;
	m_K = 0;

	m_type = UNKNOWN;
	m_connected = false;

	m_graph = NULL;
	m_tree = NULL;
	m_tree_cutset = NULL;
	m_buckets = NULL;

	m_paramParents		= UNKNOWN;
	m_paramConnectivity = UNKNOWN;
	m_paramHeight		= UNKNOWN;
	m_paramWidth		= UNKNOWN;
	m_paramSigma		= UNKNOWN;

	m_bestSolutionCost = 0;
	m_bestSolution = NULL;

	m_outputFile = NULL;

	m_levelInfo = NULL;

	m_ordering_backup = NULL;
	m_position_backup = NULL;

	m_orderingDFS = NULL;
	m_positionDFS = NULL;

	m_ordering_cutsetDFS = NULL;
	m_position_cutsetDFS = NULL;

	m_isDeterministic = NULL;


	m_descendants = NULL;

	m_parentSepSet = NULL; 
	m_parentSet = NULL; 
	m_cacheFull = NULL;
	m_cacheSep = NULL;
	m_cacheFullCounting = NULL;
	m_cacheSepCounting = NULL;
	m_advancedCacheFlagAND = NULL;
	m_advancedCacheFlagOR = NULL;

	m_wCutset = NULL;
	m_wCutsetHeight = NULL;

	
	m_treeDecomposition = NULL;
}

CProblemMixed::CProblemMixed(int n, int k) : CProblem(n, k, 0)
{
	m_Constraints = 0;

	m_type = UNKNOWN;
	m_connected = false;

	m_graph = NULL;
	m_tree = NULL;
	m_tree_cutset = NULL;
	m_buckets = NULL;

	m_paramParents		= UNKNOWN;
	m_paramConnectivity = UNKNOWN;
	m_paramHeight		= UNKNOWN;
	m_paramWidth		= UNKNOWN;
	m_paramSigma		= UNKNOWN;

	m_bestSolutionCost = 0;
	m_bestSolution = NULL;

	m_outputFile = NULL;

	m_levelInfo = NULL;

	m_ordering_backup = NULL;
	m_position_backup = NULL;

	m_orderingDFS = NULL;
	m_positionDFS = NULL;

	m_ordering_cutsetDFS = NULL;
	m_position_cutsetDFS = NULL;

	m_isDeterministic = NULL;


	m_descendants = NULL;

	m_parentSepSet = NULL; 
	m_parentSet = NULL; 
	m_cacheFull = NULL;
	m_cacheSep = NULL;
	m_cacheFullCounting = NULL;
	m_cacheSepCounting = NULL;
	m_advancedCacheFlagAND = NULL;
	m_advancedCacheFlagOR = NULL;

	m_wCutset = NULL;
	m_wCutsetHeight = NULL;

	
	m_treeDecomposition = NULL;
}


CProblemMixed::CProblemMixed(int t, int n, int k, int p, int nConstraints, int sizeConstraint, double tightness) : 
				CProblem(n, k, p), m_type(t), m_Constraints(nConstraints), 
				m_paramConstraintParents(sizeConstraint), m_tightness(tightness)
{
	m_connected = false;

	m_graph = NULL;
	m_tree = NULL;
	m_tree_cutset = NULL;
	m_buckets = NULL;
	m_bayesGraph = NULL;
	m_constraintGraph = NULL;

	m_paramParents		= UNKNOWN;
	m_paramConnectivity = UNKNOWN;
	m_paramHeight		= UNKNOWN;
	m_paramWidth		= UNKNOWN;
	m_paramSigma		= UNKNOWN;

	m_bestSolutionCost = 0;
	m_bestSolution = NULL;

	m_outputFile = NULL;

	m_levelInfo = NULL;

	m_ordering_backup = NULL;
	m_position_backup = NULL;

	m_orderingDFS = NULL;
	m_positionDFS = NULL;

	m_ordering_cutsetDFS = NULL;
	m_position_cutsetDFS = NULL;

	m_isDeterministic = NULL;

	m_descendants = NULL;

	m_parentSepSet = NULL; 
	m_parentSet = NULL; 
	m_cacheFull = NULL;
	m_cacheSep = NULL;
	m_cacheFullCounting = NULL;
	m_cacheSepCounting = NULL;
	m_advancedCacheFlagAND = NULL;
	m_advancedCacheFlagOR = NULL;

	m_wCutset = NULL;
	m_wCutsetHeight = NULL;
	
	m_treeDecomposition = NULL;
}

// Destroy the problem instance.
CProblemMixed::~CProblemMixed()
{
	if (m_graph) delete m_graph;

	if (m_bayesGraph) delete m_bayesGraph;

	if (m_constraintGraph) delete m_constraintGraph;

	if (m_tree) delete m_tree;

	if (m_tree_cutset) delete m_tree_cutset;

	if (m_buckets) delete m_buckets;

	if (m_bestSolution)	delete[] m_bestSolution;

	if (m_ordering_backup) delete[] m_ordering_backup;

	if (m_position_backup) delete[] m_position_backup;

	if (m_orderingDFS) delete[] m_orderingDFS;

	if (m_positionDFS) delete[] m_positionDFS;

	if (m_ordering_cutsetDFS) delete[] m_ordering_cutsetDFS;

	if (m_position_cutsetDFS) delete[] m_position_cutsetDFS;

	if (m_isDeterministic) delete[] m_isDeterministic;

	if (m_descendants)
	{
		for (int i = 0; i < m_N; ++i)
			m_descendants[i].clear();
		delete[] m_descendants;
	}

//#ifdef CACHE_AT_OR_NODES
	if (m_parentSet)
	{
		for (int i = 0; i < m_N; ++i)
			m_parentSet[i].clear();
		delete[] m_parentSet;
	}

	if (m_cacheFull) delete[] m_cacheFull;

	if (m_cacheFullCounting) delete[] m_cacheFullCounting;
//#endif CACHE_AT_OR_NODES


//#ifdef CACHE_AT_AND_NODES
	if (m_parentSepSet)
	{
		for (int i = 0; i < m_N; ++i)
			m_parentSepSet[i].clear();
		delete[] m_parentSepSet;
	}

	if (m_cacheSep)	delete[] m_cacheSep;

	if (m_cacheSepCounting) delete[] m_cacheSepCounting;
//#endif CACHE_AT_AND_NODES

	if (m_advancedCacheFlagAND)	delete[] m_advancedCacheFlagAND;

	if (m_advancedCacheFlagOR) delete[] m_advancedCacheFlagOR;

	function_map::iterator it = m_functionsConstraint.begin();
	for ( ; it != m_functionsConstraint.end(); ++it)
		delete (*it).second;
	m_functionsConstraint.clear();

	if (m_treeDecomposition) delete m_treeDecomposition;	

	if (m_wCutset) delete[] m_wCutset;

//	if (m_wCutsetHeight) delete[] m_wCutsetHeight;
}

void CProblemMixed::backupOrdering()
{
	if (m_ordering_backup)
		delete[] m_ordering_backup;

	m_ordering_backup = new int[m_N];
	assert(m_ordering_backup != NULL);
	memcpy(m_ordering_backup, m_ordering, m_N * sizeof(int));
}

void CProblemMixed::backupPosition()
{
	if (m_position_backup)
		delete[] m_position_backup;

	m_position_backup = new int[m_N];
	assert(m_position_backup != NULL);
	memcpy(m_position_backup, m_position, m_N * sizeof(int));
}

void CProblemMixed::restoreOrdering()
{
	if (m_ordering)
		delete[] m_ordering;

	m_ordering = new int[m_N];
	assert(m_ordering != NULL);
	memcpy(m_ordering, m_ordering_backup, m_N * sizeof(int));
}

void CProblemMixed::restorePosition()
{
	if (m_position)
		delete[] m_position;

	m_position = new int[m_N];
	assert(m_position != NULL);
	memcpy(m_position, m_position_backup, m_N * sizeof(int));
}

CConstraintTable* CProblemMixed::getConstraintFunctionByID(int id)
{
	function_map::iterator it = m_functionsConstraint.find(id);

	return (it != m_functionsConstraint.end()) ?
		((CConstraintTable*)(*it).second) : NULL;
}

void CProblemMixed::addConstraintFunction(CConstraintTable* c, bool w_id)
{
	unsigned long id = (w_id) ? ++m_Constraints : c->getID();
	if (w_id) c->setID(id);

	m_functionsConstraint.insert(make_pair(id, c));
}


int CProblemMixed::getTreeHeight()
{
	if (m_tree)
		return m_tree->height();

	return 0;
}

variable_v& CProblemMixed::getTreeDescendants(int var)
{
	return m_descendants[var];
}

bool CProblemMixed::isDescendant(int ancestor, int var)
{
	variable_v::iterator it = m_descendants[ancestor].begin();
	for (; it != m_descendants[ancestor].end(); ++it)
	{
		int v = (*it);

		if (v == var) return true;
	}
	return false;
}



int CProblemMixed::getInducedWidth()
{
	if (m_graph)
		return m_graph->getWidth();

	return 0;
}

int* CProblemMixed::getOrdering(bool dfsOrd)
{
	return ((!dfsOrd) ? m_ordering : m_orderingDFS);
}

int* CProblemMixed::getPosition(bool dfsOrd)
{
	return ((!dfsOrd) ? m_position : m_positionDFS);
}

void CProblemMixed::setOrdering(int* ordering)
{
	memcpy(m_ordering, ordering, sizeof(int) * m_N);
}

void CProblemMixed::setPosition(int* position)
{
	memcpy(m_position, position, sizeof(int) * m_N);
}

function_v& CProblemMixed::getBucket(int var, bool dfsOrd)
{
	// Safety checks.
	assert(m_buckets != NULL);

	int pos = m_position[var];
	CBucket* bkt = m_buckets->getBucketAt(pos);
	assert(bkt != NULL);

	return bkt->functions();
}

bool CProblemMixed::init()
{
	// Safety checks.
	assert(m_evidence == NULL);
	assert(m_assignment == NULL);
	assert(m_adjacencies == NULL);
	assert(m_graph == NULL);

	m_assignment = new int[m_N];
	for(int i=0; i<m_N; i++)
		m_assignment[i] = -1;	

	m_evidence = new bool[m_N];
	for(int i=0; i<m_N; i++)
		m_evidence[i] = false;	

	m_backupAssignment = new int[m_N];
	for(int i=0; i<m_N; i++)
		m_backupAssignment[i] = -1;	
	
	m_backupEvidence = new bool[m_N];
	for(int i=0; i<m_N; i++)
		m_backupEvidence[i] = false;	

	m_domains = new int[m_N];
	for (int k = 0; k < m_N; ++k) m_domains[k] = m_K;


	// Initialize the moral mixed graph.
	m_graph = new CGraph(m_N);

	// Initialize the constraint graph.
	m_constraintGraph = new CGraph(m_N);

	// Initialize the bayes graph.
	m_bayesGraph = new CGraph(m_N);


	return true;
}

// Initialize constraint propagation data structures
bool CProblemMixed::initConstraintPropagation()
{
	// initialize the domains of variables for propagation
	m_domainsPropagation = new CBitVector[m_N];
	for (int i=0; i<m_N; i++)
		m_domainsPropagation[i] = CBitVector(m_domains[i], true);

	m_validValues = new bool*[m_N];
	m_validValuesCount = new int[m_N];
	for (int i=0; i<m_N; i++)
	{
		m_validValuesCount[i] = m_domains[i];
		m_validValues[i] = new bool[m_domains[i]];
		for(int j=0; j<m_domains[i]; j++)
			m_validValues[i][j] = true;
	}

	// make lists of functions that each variable participates in
	//function_v tempListFunc = function_v(NULL);
	//for (int i=0; i<m_N; i++)
	//	m_variablesInFunctions.push_back(tempListFunc);
	m_variablesInFunctions.resize(m_N);

	function_map::iterator itf = m_functionsConstraint.begin();
	for ( ; itf != m_functionsConstraint.end(); ++itf )
	{
		CFunction* f = (*itf).second;
		int argc = f->getArgc();
		int* argv = f->getArgv();

		for (int k = 0; k < argc; k++)
		{
			int var = argv[k];
			m_variablesInFunctions[var].push_back(f);
		}
	}

	//itf = m_functions.begin();
	//for ( ; itf != m_functions.end(); ++itf )
	//{
	//	CFunction* f = (*itf).second;
	//	int argc = f->getArgc();
	//	int* argv = f->getArgv();
	//	cout << "(";
	//	for (int k = 0; k < argc; k++)
	//	{
	//		cout << argv[k] << ",";
	//	}
	//	cout << ")" << endl;
	//}


	//for (int i = 0; i < m_N; i++ )
	//{
	//	cout << endl;
	//	cout << "variable: " << i << endl;
	//	function_v::iterator it = m_variablesInFunctions[i].begin();
	//	for ( ; it != m_variablesInFunctions[i].end() ; ++it)
	//	{
	//		CFunction* f = (*it);
	//		int argc = f->getArgc();
	//		int* argv = f->getArgv();
	//		cout << "(";
	//		for (int k = 0; k < argc; k++)
	//		{
	//			cout << argv[k] << ",";
	//		}
	//		cout << ")";
	//	}
	//}


	return true;
}

// This function asserts likely evidence for the Bayesian network.
bool CProblemMixed::assertEvidence(int ne)
{
	// Safety checks.
	assert(m_N > 0);

	if (0 == ne)
		return false;

	// Reset previous evidence.
	for(int i=0; i<m_N; i++)
	{
		m_evidence[i] = false;
		m_assignment[i] = -1;
	}

	// Declare variables.
	int numVarWithValues = 0, i;

	// Generate a value for every variable.
	while (numVarWithValues < m_N) 
	{
		function_map::iterator it = m_functions.begin();
		for (; it != m_functions.end(); ++it) 
		{
			CFunction* fun = (*it).second;
			CProbabilityTable* pt = dynamic_cast<CProbabilityTable*>(fun);
			assert(fun != NULL);

			int child = pt->getChild();
			if (m_assignment[child] >= 0) continue ;
			
			// Check if all parents have a value.
			int n = pt->getNumberOfParents() ;
			for (i = 0 ; i < n ; i++) 
			{
				int p = pt->getParent(i);
				if (!isEvidence(p)) break;
			}

			// If all parents have a value, set a value for the child.
			if (i < n) continue;

			// Compute the total sum of the PT for this variable given fixed parent values.
			double x = 0.0;
			for (i = 0; i < m_domains[child]; ++i) 
			{
				setEvidence(child, i);

				double f;
				pt->getCurrentValue(&f);
				
				x += f;
			}

			double y = x * RandUniformDouble();
			x = 0.0;
			for (i = 0; i < m_domains[child]; ++i) 
			{
				setEvidence(child, i);

				double f;
				pt->getCurrentValue(&f);
				
				x += f ;
				
				if (x >= y) 
				{
					setEvidence(child, i);

					break ;
				}
			}

			++numVarWithValues;
		}
	}

	// Pick randomly evidence variables.
	variable_s evidence;
	for (int ei = 0 ; ei < ne ; ei++) 
	{
		while (1) 
		{
			int evar = RandUniform(m_N);

			// Check if already assigned.
			if (evidence.find(evar) == evidence.end())
			{
				evidence.insert(evar);
				break;
			}
		}
	}

	// Set evidence variables.
	for (int evar = 0 ; evar < m_N ; ++evar) 
	{
		if (evidence.find(evar) == evidence.end())
		{
			resetValue(evar);
			resetEvidence(evar);
		}
	}

	// Free temporary buffers.
	evidence.clear();

	return true;
}

bool CProblemMixed::create()
{
	// Safety checks.
	assert(m_ordering == NULL);
	assert(m_position == NULL);

	// Init the instance.
	init();

	switch (m_type)
	{
	case GRAPH_RANDOM:
		{
			// Create a random bayesian structure.
			// For each child set the number of parents.
			int i, j;
			int* numParents = new int[m_N];
			for (i = 0; i < m_N; ++i) 
			{
				double r = RandUniformDouble();
				if (r <= .2)
					numParents[i] = 0;
				else if (r <= .3)
					numParents[i] = 1;
				else if (r <= .5)
					numParents[i] = 2;
				else if (r <= .75)
					numParents[i] = 3;
				else if (r <= .95)
					numParents[i] = 4;
				else
					numParents[i] = 5;
				
				numParents[i] = min(numParents[i], i);
				numParents[i] = min(numParents[i], m_paramConnectivity);
			}

			// Create graph structure.
			for (int child = 0; child < m_N; ++child) 
			{
				variable_s parents;
				for (j = 0; j < numParents[child];) 
				{
					int parent;
					double r = RandUniformDouble();
					if (child > m_paramConnectivity)
						parent = (int)(i - 1 - floor(m_paramConnectivity * r));
					else
						parent = (int) floor(child * r);

					if (parents.find(parent) == parents.end())
					{
						parents.insert(parent);
						++j;
					}
				}

				// Create a function with child and parents.
				int argc = 1 + numParents[child];
				int* argv = new int[argc];
				assert(argv != NULL);
		
				// Set the function's parents.
				int pos = 0;
				variable_s::iterator it = parents.begin();
				for (; it != parents.end(); ++it) argv[pos++] = (*it);

				// Set the function's child (last position).
				argv[pos] = child;

				// Add the function to the problem instance.
				int id = child;
				CProbabilityTable* cpt = new CProbabilityTable(id, CPT, argc, argv);
				cpt->setOwner(this);
				cpt->create(PT_UNIFORM);
				addFunction(cpt);

				parents.clear();
			}

			// Create the moral graph.
			m_graph->init(m_functions);
			m_connected = m_graph->isConnected();

			break;
		}
	case GRAPH_RANDOM_FIXED:
		{
			// Create a random bayesian structure 
			// with fixed number of parents per child.

			int N = m_N;
			int P = m_paramParents;
			int C = m_C;
			int nConstraints = m_Constraints;
			m_Constraints = 0; //will be incremented again by addConstraintFunction

			int* ordering = new int[N];
			int* position = new int[N];

			// Create a random ordering of the variables.
			int i;
			for (i = 0; i < N; ++i)
			{
				ordering[i] = i;
				position[i] = i;
			}

			// Randomly, switch pairs of variables.
			for (i = 0 ; i < N ; ++i) 
			{
				int j = RandUniform(N);
			
				// Switch variable i and j.
				int k = ordering[j];
				ordering[j] = ordering[i];
				ordering[i] = k;

				position[ordering[j]] = j;
				position[ordering[i]] = i;
			}

			bool* cpts = new bool[N];
			memset(cpts, false, N * sizeof(bool));

			//create C CPTs
			int count = 0;
			while (count < C)
			{
				// Pick the child.
				int child = ordering[RandUniform(C)];

				// Check if variable was visited
				if (cpts[child]) continue;
				cpts[child] = 1;

				// Pick P parents for the child.
				int numHigherVars = N - position[child] - 1;
				// Notice : number of parents cannot be larger than 'numHigherVars'.
				int numParents = P;
				if (numParents > numHigherVars) 
					numParents = numHigherVars;
				
				variable_s parents;
				for (i = 0; i < numParents;) 
				{
					int parent = ordering[position[child] + 1 + RandUniform(numHigherVars)];
					// Check that this parent is not the same as the child.
					if (child == parent) continue ;
					// Check other parents.
					if (parents.find(parent) == parents.end())
					{
						parents.insert(parent);
						++i;
					}
				}

				// Create a function with child and parents.
				int argc = 1 + numParents;
				int* argv = new int[argc];
				assert(argv != NULL);

				// Set the function's parents.
				int pos = 0;
				variable_s::iterator it = parents.begin();
				for (; it != parents.end(); ++it) argv[pos++] = (*it);

				// Set the function's child (last position).
				argv[pos] = child;

				// Add the function to the problem instance.
				int id = child;
				CProbabilityTable* cpt = new CProbabilityTable(id, CPT, argc, argv);
				cpt->setOwner(this);
				cpt->create(PT_UNIFORM);
				addFunction(cpt);

				parents.clear();

				++count;
			}

///*
			// Create priors.
			for (i = 0; i < N; ++i)
			{
				if (cpts[i]) continue;

				int argc = 1;
				int* argv = new int[argc];
				assert(argv != NULL);

				// Set the function's child (last position).
				argv[0] = i;
					
				// Add the function to the problem instance.
				int id = i;
				CProbabilityTable* cpt = new CProbabilityTable(id, CPT, argc, argv);
				cpt->setOwner(this);
				cpt->create(PT_UNIFORM);
				addFunction(cpt);

			}
			//end create C CPTs
//*/

			//now create nConstraints constraint tables
			count = 0;
			while (count < nConstraints)
			{
				variable_s parents;
				for (i = 0; i < m_paramConstraintParents;) 
				{
					int parent = ordering[RandUniform(N)];			
					// Check other parents.
					if (parents.find(parent) == parents.end())
					{
						parents.insert(parent);
						++i;
					}
				}
				
				// Create a constraint with parents as scope
				int argc = m_paramConstraintParents;
				int* argv = new int[argc];
				assert(argv != NULL);

				// Set the constraint's scope.
				int pos = 0;
				variable_s::iterator it = parents.begin();
				for (; it != parents.end(); ++it) argv[pos++] = (*it);

				// Add the function to the problem instance.
				int id = count;
				CConstraintTable* constraint = new CConstraintTable(id, CONSTRAINT, argc, argv);
				constraint->setOwner(this);
				constraint->create(m_tightness);
				addConstraintFunction(constraint);

				parents.clear();

				++count;


			}
			//end create nConstraints constraint tables


			// Create the moral mixed graph.
			m_graph->init(m_functions, m_functionsConstraint);
			m_connected = m_graph->isConnected();

			
			// Create the moral bayes graph.
			m_bayesGraph->init(m_functions);
			m_bayesConnected = m_bayesGraph->isConnected();

			// Create the constraint graph.
			m_constraintGraph->init(m_functionsConstraint);
			m_constraintConnected = m_constraintGraph->isConnected();

			delete[] ordering;
			delete[] position;
			delete[] cpts;

			break;
		}
	case GRAPH_GRID:
		{
			break;
		}
	case GRAPH_CODING:
		{
			break;
		}
	};

	return true;
}


void CProblemMixed::createCustom()
{
	// Create a problem.
	assert(m_N > 0);
	assert(m_K > 0);

	m_N = 7;
	m_K = 2;

	init();

	// Create CPTs.
	int argc;
	int* argv = NULL;
	CProbabilityTable* cpt = NULL;

	// CPT 0.
	argc = 1;
	argv = new int[argc];
	argv[0] = 0;

	cpt = new CProbabilityTable(0, CPT, argc, argv);
	cpt->setOwner(this);
	cpt->create(PT_UNIFORM);
	addFunction(cpt);

	// CPT 1.
	argc = 2;
	argv = new int[argc];
	argv[0] = 0;
	argv[1] = 1;

	cpt = new CProbabilityTable(1, CPT, argc, argv);
	cpt->setOwner(this);
	cpt->create(PT_UNIFORM);
	addFunction(cpt);

	// CPT 2.
	argc = 2;
	argv = new int[argc];
	argv[0] = 0;
	argv[1] = 2;

	cpt = new CProbabilityTable(2, CPT, argc, argv);
	cpt->setOwner(this);
	cpt->create(PT_UNIFORM);
	addFunction(cpt);

	// CPT 3.
	argc = 2;
	argv = new int[argc];
	argv[0] = 1;
	argv[1] = 3;

	cpt = new CProbabilityTable(3, CPT, argc, argv);
	cpt->setOwner(this);
	cpt->create(PT_UNIFORM);
	addFunction(cpt);

	// CPT 4.
	argc = 3;
	argv = new int[argc];
	argv[0] = 0;
	argv[1] = 1;
	argv[2] = 4;

	cpt = new CProbabilityTable(4, CPT, argc, argv);
	cpt->setOwner(this);
	cpt->create(PT_UNIFORM);
	addFunction(cpt);

	// CPT 5.
	argc = 2;
	argv = new int[argc];
	argv[0] = 2;
	argv[1] = 5;

	cpt = new CProbabilityTable(5, CPT, argc, argv);
	cpt->setOwner(this);
	cpt->create(PT_UNIFORM);
	addFunction(cpt);

	// CPT 6.
	argc = 2;
	argv = new int[argc];
	argv[0] = 2;
	argv[1] = 6;

	cpt = new CProbabilityTable(6, CPT, argc, argv);
	cpt->setOwner(this);
	cpt->create(PT_UNIFORM);
	addFunction(cpt);

	// Create the moral graph.
	m_graph->init(m_functions);
	m_connected = m_graph->isConnected();

}

// This function preprocesses the network's graph. It creates the
// variable ordering, computes the induced width and from the induced
// graph it creates the rooted tree arrangement [Bayardo96].
bool CProblemMixed::preprocess(int voType, bool withTree)
{
	// Safety checks.
	assert(m_graph != NULL);

	if (!m_connected)
	{
		printf("\n not connected");
		return false;	// not connected.
	}

	// Create variable ordering, induced graph.
	m_graph->order(voType, m_ordering, m_position);

	// Check if we have to create the rooted tree.
	if (withTree)
	{
		// Create the rooted tree arrangement from the induced graph.
		assert(m_tree == NULL);

		m_tree = new CLegalTree(m_N, m_ordering, m_position);
		m_tree->create(m_graph);
		m_tree->order(m_orderingDFS, m_positionDFS);

		m_descendants = new variable_v[m_N];

		assert(m_descendants != NULL);

		for (int var = 0; var < m_N; ++var)	
			m_tree->descendants(var, m_descendants[var]);

		memcpy(m_ordering, m_orderingDFS, m_N*sizeof(int));
		memcpy(m_position, m_positionDFS, m_N*sizeof(int));
	}

	
	// robert: Create the parent separator set
	// this includes earlier vars connected to vars below, union the current var
	// it does not include the earlier vars only connected to the current, but not to ones below
	
	parentsetCreate();
	
	parentsetInit();
		
	
	// Create the bucket structure.
	if (m_buckets)
	{
		delete m_buckets;
		m_buckets = NULL;
	}

	m_buckets = new CBucketStruct(this);

	assert(m_buckets != NULL);

	m_buckets->init(withTree);

	// Robert
	// Set the m_highestVarInScope for all the constraints
	// This is equivalent to making a bucket structure for the constraints;

	setHighestVarInScope();

	// Reset original functions.
	function_map& funs = m_functionsConstraint;
	function_map::iterator it = funs.begin();
	for (; it != funs.end(); ++it)
	{
		CFunction* fun = (*it).second;
		fun->setUsed(false);
	}

	return true;
}

// This function preprocesses the network's graph. It creates the
// variable ordering, computes the induced width and from the induced
// graph it creates the rooted tree arrangement [Bayardo96].
//bool CProblemMixed::preprocess2(int voType, bool withTree)
//{
//	// Safety checks.
//	assert(m_graph != NULL);
//
//	if (!m_connected)
//	{
//		printf("\n not connected");
//		return false;	// not connected.
//	}
//
//	// Assert as evidence all singleton variables.
//	// Reindex variables, scopes, moral graph.
//	// removeEvidence();
//	if (!m_connected)
//	{
//		cout << "\nAdjusted graph is not connected!";
//		cout << "\nAborting ..." << endl;
//		return false;	// not connected.
//	}
//
//	int i, var, pos;
//
//	vector<int> order;
//	int width;
//
//	if (g_useSuperlinkOrder)
//	{
//		// Set the new superlink order as default order.
//		assert(m_N == (int)m_superlinkOrder.size());
//
////*		// ---------------
//		// Move superlink cutset (search starts with cutset).
//		vector<int> temp;
//		for (i = 0; i < m_N; ++i)
//		{
//			var = m_superlinkOrder[i];
//			vector<int>::iterator it1 = find(m_superlinkCutset.begin(), m_superlinkCutset.end(), var);
//			if (it1 == m_superlinkCutset.end())
//			{
//				temp.push_back(var);
//			}
//		}
//
//		for (i = 0; i < (int)m_superlinkCutset.size(); ++i)
//			temp.push_back(m_superlinkCutset[i]);
//
//		m_superlinkOrder.clear();
//		copy(temp.begin(), temp.end(), back_inserter(m_superlinkOrder));
//		temp.clear();
//		
//		// ---------------
////*/
//
//		m_ordering = new int[m_N];
//		m_position = new int[m_N];
//		pos = m_N - 1;
//		for (i = 0; i < m_N; ++i)
//		{
//			var = m_superlinkOrder[i];
//			m_ordering[pos] = var;
//			m_position[var] = pos;
//			--pos;
//		}
//	}
//	else
//	{
//		// Create variable ordering, induced graph.
//
//		// Create the graph hash.
//		function_map::iterator itf = m_functions.begin();
//		for (; itf != m_functions.end(); ++itf)
//		{
//			CFunction* f = (*itf).second;
//			m_graph2.addClique(f->getArgc(), f->getArgv());
//		}
//
//		// Create the elimination order.
//		m_graph2.eliminate(voType, order, m_domains);
//		m_ordering = new int[m_N];
//		m_position = new int[m_N];
//		int pos = m_N - 1;
//		for (int i = 0; i < m_N; ++i)
//		{
//			int var = order[i];
//			m_ordering[pos] = var;
//			m_position[var] = pos;
//			--pos;
//		}
//	}
//
//	// Check if we have to create the rooted tree.
//	if (withTree)
//	{
//				// Create induced graph.
//		m_graph->backup();
//		m_graph->createInduced(m_ordering, m_position);
//
//		// Create pseudo-tree.
//		m_tree = new CLegalTree(m_N, m_ordering, m_position);
//		m_tree->create(m_graph);
//
//		// Create the DFS order from the pseudo-tree.
//		int *ordDFS = NULL, *posDFS = NULL;
//		m_tree->order(ordDFS, posDFS);
//
//		m_descendants = new variable_v[m_N];
//		assert(m_descendants != NULL);
//		for (int var = 0; var < m_N; ++var)	
//			m_tree->descendants(var, m_descendants[var]);
//
//		// Record the new ordering.
//		memcpy(m_ordering, ordDFS, m_N*sizeof(int));
//		memcpy(m_position, posDFS, m_N*sizeof(int));
//
//		// Create the induced graph, wrt new ordering.
//		m_graph->restore();
//		m_graph->createInduced(m_ordering, m_position);
//		
//		delete[] ordDFS;
//		delete[] posDFS;
//	}
//	else
//	{
//		m_graph2.induce(order, width);
//		m_graph->setWidth(width);
//	}
//
//	// robert: Create the parent separator set
//	// this includes earlier vars connected to vars below, union the current var
//	// it does not include the earlier vars only connected to the current, but not to ones below
//	
//	parentsetCreate();
//	
//	parentsetInit();
//	
//	// Create the bucket structure.
//	if (m_buckets)
//	{
//		delete m_buckets;
//		m_buckets = NULL;
//	}
//
//	m_buckets = new CBucketStruct(this);
//	assert(m_buckets != NULL);
//	m_buckets->init(false);
//
//	// Robert
//	// Set the m_highestVarInScope for all the constraints
//	// This is equivalent to making a bucket structure for the constraints;
//
//	setHighestVarInScope();
//
//	// Reset original functions.
//	function_map& funs = m_functionsConstraint;
//	function_map::iterator it_f = funs.begin();
//	for (; it_f != funs.end(); ++it_f)
//	{
//		CFunction* fun = (*it_f).second;
//		fun->setUsed(false);
//	}
//
//	order.clear();
//
//	return true;
//}
//

// This function preprocesses the network's graph. It creates the
// variable ordering, computes the induced width and from the induced
// graph it creates the rooted tree arrangement [Bayardo96].
// 
// This function simulates Geigers elimination plus conditioning algorithm
// When inducing the graph, we skip the conditioning variables
// then move conditioning set at the beginning
// then create the pseudo-tree in the new ordering, on the special-induced graph

void CProblemMixed::markBarren_removeEvidence_makeConnectedComponents()
{
	// Identify barren nodes (that don't have evidence descendants)
	// initially all are barren, then evidence is propagated up
	bool* barrenVariables = new bool[m_N];
	for (int i=0; i<m_N; i++)
		barrenVariables[i] = false;
	
	markBarrenVariables(barrenVariables);

	cout << "\n initial number of variables   : " << m_N;
	removeEvidence(barrenVariables);
	cout << "\n number of non-barren variables: " << m_N; 

	makeConnectedComponents();

	delete[] barrenVariables;
}

bool CProblemMixed::preprocessAdaptiveCaching(int voType, bool withTree)
{	
	if (m_N == 0)
	{
		cout << "\n probability of the evidence is: " << globalCostFromEvidence;
		getchar();
		exit(1);
	}

	// Create new moral graph.
	m_graph->init(m_functions);
	m_connected = m_graph->isConnected();

	// Safety checks.
	assert(m_graph != NULL);

	if (!m_connected)
	{
		printf("\n not connected");
		return false;	// not connected.
	}

	int i, var, pos;

	vector<int> order;
	int width;

	bool isDeterministic = false;
	// Find deterministic variables
	if (voType == VO_TOPOLOGICAL)
		isDeterministic = findDeterministicVariables(); 

	//int count = 0;
	//cout << endl << " This network" << (isDeterministic ? " has" : " does not have") << " determinism" << endl;
	//if (isDeterministic)
	//{
	//	cout << " The deterministic variables are:" << endl;
	//	for (int i=0 ; i<m_N ; i++)
	//	{
	//		if (deterministic[i])
	//		{
	//			cout << i << " ";
	//			count++;
	//		}
	//	}
	//}
	//cout << endl << " N = " << m_N << endl << "N_deterministic = " << count;
	
	// Create variable ordering, induced graph.

	variable_s* parents = new variable_s[m_N];
	// Create the graph hash.
	function_map::iterator itf = m_functions.begin();
	for (; itf != m_functions.end(); ++itf)
	{
		CFunction* f = (*itf).second;
		m_graph2.addClique(f->getArgc(), f->getArgv());

		if(isDeterministic)
		{
			int argc = f->getArgc();
			int* argv = f->getArgv();
			if ( m_isDeterministic[argv[argc-1]] )
				for (int i = 0; i < argc - 1 ; i++)
					parents[argv[argc-1]].insert(argv[i]);
		}
	}

	//for (int i=0; i<m_N ; i++)
	//{
	//	cout << endl << "parents[" << i << "]: ";
	//	variable_s::iterator it = parents[i].begin();
	//	for ( ; it != parents[i].end(); ++it )
	//		cout << (*it) << " ";
	//}	
	
	// Create the elimination order.
	
	//if(isDeterministic)
		//m_graph2.addAncestralEdges(parents);

	if(voType == VO_TOPOLOGICAL)
		findTopologicalOrdering(order);
	else
		m_graph2.eliminate(voType, order, m_domains, m_isDeterministic, parents);

	m_ordering = new int[m_N];
	m_position = new int[m_N];
	pos = m_N - 1;
	for (int i = 0; i < m_N; ++i)
	{
		int var = order[i];
		m_ordering[pos] = var;
		m_position[var] = pos;
		--pos;
	}

	//cout << "\n position of variables:\n";
	//for (int i=0; i<m_N; i++)
	//	cout << " [" << i << "]" << m_position[i];

	//cout << "\n ordering of variables:\n";
	//for (int i=0; i<m_N; i++)
	//	cout << "\n [" << i << "]" << m_ordering[i];

	// move deterministic variables immediately after parents are instantiated
	// deterministic variables will be substituted in contexts, if they are fully determined
	int* contextSubstitution = new int[m_N];
	for (int i=0; i< m_N; i++)
		contextSubstitution[i] = i;

	//if(isDeterministic)
	//	fixOrderingDeterministic2(parents, contextSubstitution);

	if(isDeterministic)
		for (int i=0; i<m_N; i++)
		{
			cout << "\n position[" << i << "]: " << m_ordering[i] << " parents are: ";
			if (m_isDeterministic[m_ordering[i]])
			{
				variable_s::iterator it = parents[m_ordering[i]].begin();
				for ( ; it != parents[m_ordering[i]].end(); it++)
					cout << (*it) << " [" << m_position[*it] << "] ";
			}
			else
				cout << "                     -------------------- non deterministic" ;
			cout << "  subst with -> " << contextSubstitution[m_ordering[i]];
			/*	if (i%100 == 99)
				break;*/
				//getchar();
		}
	
	

	// Check if we have to create the rooted tree.m
	if (withTree)
	{
		// Create induced graph.
		m_graph->backup();
		m_graph->createInduced(m_ordering, m_position);

		// Create pseudo-tree.
		m_tree = new CLegalTree(m_N, m_ordering, m_position);
		m_tree->create(m_graph);

		// Create the DFS order from the pseudo-tree.
		int *ordDFS = NULL, *posDFS = NULL;
		m_tree->order(ordDFS, posDFS);

		m_descendants = new variable_v[m_N];
		assert(m_descendants != NULL);
		for (int var = 0; var < m_N; ++var)	
			m_tree->descendants(var, m_descendants[var]);

		// Record the new ordering.
		memcpy(m_ordering, ordDFS, m_N*sizeof(int));
		memcpy(m_position, posDFS, m_N*sizeof(int));

		// Create the induced graph, wrt new ordering.
		m_graph->restore();
		m_graph->createInduced(m_ordering, m_position);

		//m_tree->print();
		
		delete[] ordDFS;
		delete[] posDFS;
	}
	else
	{
		m_graph2.induce(order, width);
		m_graph->setWidth(width);
	}

	// Robert: Create the parent separator set
	// this includes earlier vars connected to vars below, union the current var
	// it does not include the earlier vars only connected to the current, but not to ones below
	
	// for deterministic networks, find the top chain of deterministic variables and exclude from contexts
	/*variable_s topDeterministicVars;
	if(isDeterministic)
		for (int i=0; i<m_N; i++)
		{
			int v = m_ordering[i];
			if (m_isDeterministic[v])
				topDeterministicVars.insert(v);
			else
				break;
		}*/

	parentsetCreate();
	
	//if(isDeterministic)
	//	parentsetInit(contextSubstitution);
	//else
		parentsetInit();

	//topDeterministicVars.clear();

	//for (int i=0; i<m_N; i++)
	//{
	//	int var = m_ordering[i];
	//	if (m_isDeterministic[var]) continue;

	//	cout << "\n context[" << var << "] = ";

	//	variable_v::iterator it = m_parentSet[var].begin();
	//	for ( ; it != m_parentSet[var].end(); it++)
	//		cout << (*it) << " ";

	//	//if (i%100 ==99) break;
	//}
	
	// Create the bucket structure.
	if (m_buckets)
	{
		delete m_buckets;
		m_buckets = NULL;
	}

	m_buckets = new CBucketStruct(this);
	assert(m_buckets != NULL);
	m_buckets->init(false);

	// Robert
	// Set the m_highestVarInScope for all the constraints
	// This is equivalent to making a bucket structure for the constraints;

	setHighestVarInScope();

	 //Reset original functions.
	function_map& funs = m_functionsConstraint;
	function_map::iterator it_f = funs.begin();
	for (; it_f != funs.end(); ++it_f)
	{
		CFunction* fun = (*it_f).second;
		fun->setUsed(false);
	}

	
	for (int i=0; i<m_N; i++)
		parents[i].clear();
	delete[] parents; 

	order.clear();

	return true;
}


bool CProblemMixed::preprocessCutsetTree()
{
	// change the tree
	m_treeBackup = getTree();
	setTree(m_tree_cutset);
	m_tree_cutset = m_treeBackup;

//	m_tree->print();

	backupOrdering();
	backupPosition();

	m_tree->order(m_orderingDFS, m_positionDFS);
	memcpy(m_ordering, m_orderingDFS, m_N*sizeof(int));
	memcpy(m_position, m_positionDFS, m_N*sizeof(int));

	// clean m_graph (remove old induced edges) and recreate induced graph
	m_graph->cleanConnections();
	
	m_graph->width(m_ordering, m_position);

	if (m_descendants)
		delete[] m_descendants;
	m_descendants = new variable_v[m_N];
	assert(m_descendants != NULL);

	for (int var = 0; var < m_N; ++var)	
		m_tree->descendants(var, m_descendants[var]);



	// robert: Create the parent separator set
	// this includes earlier vars connected to vars below, union the current var
	// it does not include the earlier vars only connected to the current, but not to ones below
	parentsetCreate();

	parentsetInit();
		
	// Create the bucket structure.
	if (m_buckets)
	{
		delete m_buckets;
		m_buckets = NULL;
	}

	m_buckets = new CBucketStruct(this);
	assert(m_buckets != NULL);

	m_buckets->init(false);

	// Robert
	// Set the m_highestVarInScope for all the constraints
	// This is equivalent to making a bucket structure for the constraints;

	setHighestVarInScope();

	// Reset original functions.
	function_map& funs = m_functionsConstraint;
	function_map::iterator it = funs.begin();
	for (; it != funs.end(); ++it)
	{
		CFunction* fun = (*it).second;
		fun->setUsed(false);
	}

	return true;
}

void CProblemMixed::setHighestVarInScope()
{
	for (int i = m_N - 1; i >= 0; --i)
	{
		int var = m_ordering[i];
		
		// Look for all constraints that have 'var' "highest"
		// in their scope and set their m_highestVarInScope
		function_map& funs = m_functionsConstraint;
		function_map::iterator it = funs.begin();
		for (; it != funs.end(); ++it)
		{
			CFunction* fun = (*it).second;
			if (fun->isUsed()) continue;

			if (fun->isMemberOf(var))
			{
				((CConstraintTable*)fun)->set_highestVarInScope(var);
				fun->setUsed(true);
			}
		}
	}	
}

void CProblemMixed::parentsetCreate()
{
	//if (CACHE_AT_OR_NODES)
	{
		if (m_parentSet)
			delete[] m_parentSet;
		m_parentSet = new variable_v[m_N];
	}

	//if (CACHE_AT_AND_NODES)
	{
		if(m_parentSepSet)
			delete[] m_parentSepSet;
		m_parentSepSet = new variable_v[m_N];
	}
}

void CProblemMixed::parentsetDestroy()
{
	if (m_parentSepSet)
	{
		for (int i = 0; i < m_N; ++i)
			m_parentSepSet[i].clear();
		delete[] m_parentSepSet;
	}

	if (m_parentSet)
	{
		for (int i = 0; i < m_N; ++i)
			m_parentSet[i].clear();
		delete[] m_parentSet;
	}
}

void CProblemMixed::parentsetInit()
{
	int var; 


	//if (CACHE_AT_OR_NODES)
	{
		for(var = 0 ; var < m_N; var++)
		{
			for(int index_ancestor = 0; index_ancestor < m_position[var];  index_ancestor++)
			{	
				int ancestor = m_ordering[index_ancestor];
				if (ancestor == var) continue;
				if(!isDescendant(ancestor, var))
					continue;

				variable_v::iterator it = m_descendants[var].begin();
				for (; it != m_descendants[var].end(); ++it)
				{
					int v = (*it);
	//					if (v == var) continue;
					if ( m_graph->isConnected(ancestor, v) )
					{
						//add ancestor to parent set of 'var'
						m_parentSet[var].push_back(ancestor);
						break;
					}						
				}				
			}
		}
	}


	//if (CACHE_AT_AND_NODES)
	{
		for(var = 0 ; var < m_N; var++)
		{
			for(int index_ancestor = 0; index_ancestor < m_position[var];  index_ancestor++)
			{	
				int ancestor = m_ordering[index_ancestor];
				if (ancestor == var) continue;
				if(!isDescendant(ancestor, var))
					continue;

				variable_v::iterator it = m_descendants[var].begin();
				for (; it != m_descendants[var].end(); ++it)
				{
					int v = (*it);
					if (v == var) continue;
					if ( m_graph->isConnected(ancestor, v) )
					{
						//add ancestor to parent set of 'var'
						m_parentSepSet[var].push_back(ancestor);
						break;
					}						
				}				
			}
			// this is for separator set!
			m_parentSepSet[var].push_back(var);
		}
	}
}



void CProblemMixed::parentsetInit(int* substitution)
{
	// cache at and nodes no longer supported;

	int var; 

	//if (CACHE_AT_OR_NODES)
	{
		for(var = 0 ; var < m_N; var++)
		{
			for(int index_ancestor = 0; index_ancestor < m_position[var];  index_ancestor++)
			{	
				int ancestor = m_ordering[index_ancestor];
				if (ancestor == var) continue;
				if(!isDescendant(ancestor, var)) continue;
				if (substitution[ancestor] == -1) continue;

				variable_v::iterator it = m_descendants[var].begin();
				for (; it != m_descendants[var].end(); ++it)
				{
					int v = (*it);
	//					if (v == var) continue;
					if ( m_graph->isConnected(ancestor, v) )
					{
						//add ancestor to parent set of 'var' 
						// see if not already in, because we use substitutions
						int anc = substitution[ancestor];
						bool isIn = false;
						variable_v::iterator itcont = m_parentSet[var].begin();
						for ( ; itcont != m_parentSet[var].end(); itcont++ )						
							if (anc == (*itcont))
							{
								isIn = true;
								break;
							}

						if (!isIn)
							m_parentSet[var].push_back(anc);
						break;
					}						
				}				
			}
		}
	}


	//if (CACHE_AT_AND_NODES)
	{
		for(var = 0 ; var < m_N; var++)
		{
			for(int index_ancestor = 0; index_ancestor < m_position[var];  index_ancestor++)
			{	
				int ancestor = m_ordering[index_ancestor];
				if (ancestor == var) continue;
				if(!isDescendant(ancestor, var))
					continue;

				variable_v::iterator it = m_descendants[var].begin();
				for (; it != m_descendants[var].end(); ++it)
				{
					int v = (*it);
					if (v == var) continue;
					if ( m_graph->isConnected(ancestor, v) )
					{
						//add ancestor to parent set of 'var'
						m_parentSepSet[var].push_back(ancestor);
						break;
					}						
				}				
			}
			// this is for separator set!
			m_parentSepSet[var].push_back(var);
		}
	}
}


void CProblemMixed::cacheInit(int scopeBound)
{
	int var;

	if (CACHE_AT_AND_NODES)		
	{
		assert(m_parentSepSet != NULL);

		m_cacheSep = new CCacheTable*[m_N];
		assert(m_cacheSep != NULL);

		if (COUNTING)
		{
			m_cacheSepCounting = new CCacheTable*[m_N];
			assert(m_cacheSepCounting != NULL);
		}
		
		for(var=0; var < m_N; var++)
		{
			int argc = m_parentSepSet[var].size();		
			int *argv, *argvCounting;

			if(argc == 0)
			{
				argv = NULL;
				argvCounting = NULL;
			}
			else
			{
				if (CUTSET_CACHE)
					argc = min (argc, scopeBound);

				argv = new int[argc];
				argvCounting = new int[argc]; // for counting sol
			
				variable_v::reverse_iterator it = m_parentSepSet[var].rbegin();
				int index = argc;
				for ( ; index > 0 ;index--)					
				{
					int v = (*it);
					argv[index-1] = v;
					argvCounting[index-1] = v; // for counting sol
					it++;
				}
			}

			// brute force kind of code, to cancel the instructions above
			// and have linear space exploration of the cutset
			if (BRUTE_FORCE_W_CUTSET && (m_wCutset[var] > scopeBound ) )
			{
				argc = scopeBound + 1;;
				if (argv) delete[] argv;
				if (argvCounting) delete[] argvCounting;
				argv = new int[argc];
				argvCounting = new int[argc];
				for(int i=0; i<m_N; i++)
				{
					argv[i] = 0;
					argvCounting[i] = 0;
				}
			}

			// check for dead cache: it is enough to check that 
			// |context[var]| = |context[parent[var]]| + 1; if yes, then dead-cache
			int position = m_tree->getPosition()[var];
			CLegalTreeNode* varNode = m_tree->getNodes()[position];
			if (varNode->variable() != var)
			{
				cout << "problem with dead cache allocation";
				exit(1);
			}
			
			CLegalTreeNode* parentNode = varNode->parent();
			if (parentNode != NULL)
			{				
				int parent = parentNode->variable();
				if( m_parentSepSet[var].size() == (m_parentSepSet[parent].size() + 1) )
				{
					argc = scopeBound + 1;;
					if (argv) delete[] argv;
					argv = new int[argc];
					for(int i=0; i<m_N; i++)
						argv[i] = 0;					
				}
			}


/*			// check for dead cache (when context[var] includes context[parent[var]]
			int position = m_tree->getPosition()[var];
			CLegalTreeNode* varNode = m_tree->getNodes()[position];
			if (varNode->variable() != var)
				exit(1);

			if (varNode->parent() != NULL)
			{
				int parent = varNode->parent()->variable();
				variable_v childContext;
				variable_v parentContext;
				copy(m_parentSepSet[var].begin(), m_parentSepSet[var].end(), inserter(childContext, childContext.begin() ));
				copy(m_parentSepSet[parent].begin(), m_parentSepSet[parent].end(), inserter(parentContext, parentContext.begin() ));

				sort(childContext.begin(), childContext.end());
				sort(parentContext.begin(), parentContext.end());
				if ( includes( childContext.begin(), childContext.end(),
					parentContext.begin(), parentContext.end() ) )
				{
					argc = scopeBound + 1;;
					if (argv) delete[] argv;
					if (argvCounting) delete[] argvCounting;
					argv = new int[argc];
					argvCounting = new int[argc];
					memset(argv, 0, m_N * sizeof(int));
					memset(argvCounting, 0, m_N * sizeof(int));						
				}
				childContext.clear();
				parentContext.clear();
			}
*/			

			CCacheTable* tempCache = new CCacheTable(var, CACHE, argc, argv);
			tempCache->setOwner(this);
			m_cacheSep[var] = tempCache;
			if (argc == 0)
				m_cacheSep[var]->setConstant(true);

			// allocate space for the cache table - should not allocate for brute-force w-cutset
 			if( (argc <= scopeBound) )
				m_cacheSep[var]->init();	

			if (COUNTING)
			{
				CCacheTable* tempCacheCounting = new CCacheTable(var, CACHE, argc, argvCounting);
				tempCacheCounting->setOwner(this);
				m_cacheSepCounting[var] = tempCacheCounting;
				if (argc == 0)
					m_cacheSepCounting[var]->setConstant(true);

				// allocate space for the cache table
				if( (m_wCutset[var] <= scopeBound) && (argc <= scopeBound) )
					m_cacheSepCounting[var]->init();
			}
			else
				if(argvCounting) delete[] argvCounting;
		}
	} // end CACHE_AT_AND_NODES


	if (CACHE_AT_OR_NODES)
	{
		assert(m_parentSet != NULL);

		m_cacheFull = new CCacheTable*[m_N];
		assert(m_cacheFull != NULL);

		if (COUNTING)
		{
			m_cacheFullCounting = new CCacheTable*[m_N];
			assert(m_cacheFullCounting != NULL);
		}
		
		for(var=0; var < m_N; var++)
		{
			int argc = m_parentSet[var].size();		
			int *argv, *argvCounting;

			if(argc == 0)
			{
				argv = NULL;
				argvCounting = NULL;
			}
			else
			{
				if (CUTSET_CACHE)
					argc = min (argc, scopeBound);

				argv = new int[argc];
				argvCounting = new int[argc]; // for counting sol
			
				variable_v::reverse_iterator it = m_parentSet[var].rbegin();
				int index = argc;
				for ( ; index > 0 ;index--)					
				{
					int v = (*it);
					argv[index-1] = v;
					argvCounting[index-1] = v; // for counting sol
					it++;
				}			
			}
			
			// brute force kind of code, to cancel the instructions above
			if (BRUTE_FORCE_W_CUTSET && (m_wCutset[var] > scopeBound ) )
			{
				argc = max(scopeBound,m_N) + 1;;
				if (argv) delete[] argv;
				if (argvCounting) delete[] argvCounting;
				argv = new int[argc];
				argvCounting = new int[argc];
				for(int i=0; i<m_N; i++)
				{
					argv[i] = 0;
					argvCounting[i] = 0;
				}
			}

			if (m_isDeterministic && m_isDeterministic[var])
			{
				argc = max(scopeBound,m_N) + 1;;
				if (argv) delete[] argv;
				argv = new int[argc];
				for(int i=0; i<m_N; i++)
					argv[i] = 0;	
			}
			else
			{
				// check for dead cache: it is enough to check that 
				// |context[var]| = |context[parent[var]]| + 1; if yes, then dead-cache
				int position = m_tree->getPosition()[var];
				CLegalTreeNode* varNode = m_tree->getNodes()[position];
				if (varNode->variable() != var)
					exit(1);
				
				CLegalTreeNode* parentNode = varNode->parent();
				if (parentNode != NULL)
				{				
					int parent = parentNode->variable();
					if( m_parentSet[var].size() == (m_parentSet[parent].size() + 1) )
					{
						argc = max(scopeBound,m_N) + 1;;
						if (argv) delete[] argv;
						argv = new int[argc];
						for(int i=0; i<m_N; i++)
							argv[i] = 0;					
					}
				}
			}
			
			/*			// check for dead cache (when context[var] includes context[parent[var]]
			// Actually, it is enough to check that |context[var]| = |context[parent[var]]| + 1; if yes, => dead-cache
			int position = m_tree->getPosition()[var];
			CLegalTreeNode* varNode = m_tree->getNodes()[position];
			if (varNode->variable() != var)
				exit(1);

			if (varNode->parent() != NULL)
			{
				int parent = varNode->parent()->variable();
				variable_v childContext;
				variable_v parentContext;
				copy(m_parentSet[var].begin(), m_parentSet[var].end(), inserter(childContext, childContext.begin() ));
				copy(m_parentSet[parent].begin(), m_parentSet[parent].end(), inserter(parentContext, parentContext.begin() ));

				sort(childContext.begin(), childContext.end());
				sort(parentContext.begin(), parentContext.end());
				if ( includes( childContext.begin(), childContext.end(),
					parentContext.begin(), parentContext.end() ) )
				{					
					argc = scopeBound + 1;;
					if (argv) delete[] argv;
					if (argvCounting) delete[] argvCounting;
					argv = new int[argc];
					argvCounting = new int[argc];
					memset(argv, 0, m_N * sizeof(int));
					memset(argvCounting, 0, m_N * sizeof(int));
				}
				childContext.clear();
				parentContext.clear();
			}
*/		

			CCacheTable* tempCache = new CCacheTable(var, CACHE, argc, argv);
			tempCache->setOwner(this);
			m_cacheFull[var] = tempCache;
			if (argc == 0)
				m_cacheFull[var]->setConstant(true);

			// allocate space for the cache table
 			if( (argc <= scopeBound) )
				m_cacheFull[var]->init();	

			if (COUNTING)
			{
				CCacheTable* tempCacheCounting = new CCacheTable(var, CACHE, argc, argvCounting);
				tempCacheCounting->setOwner(this);
				m_cacheFullCounting[var] = tempCacheCounting;
				if (argc == 0)
					m_cacheFullCounting[var]->setConstant(true);

				// allocate space for the cache table
				if((m_wCutset[var] <= scopeBound) && (argc <= scopeBound))
					m_cacheFullCounting[var]->init();	
			}
			else
				if(argvCounting) delete[] argvCounting;
		}
	} //END CACHE_AT_OR_NODES

}

// this function is to clean the cache structure in between different scopeBound runs 
// then need to run cacheInit for the next scopeBound
void CProblemMixed::cacheDestroy()
{
	int i;

	if (CACHE_AT_OR_NODES)
	{
		if (m_cacheFull)
		{
			for (i=0; i<m_N; i++)
				delete m_cacheFull[i];
			delete[] m_cacheFull;
			m_cacheFull = NULL;
		}
		if (m_cacheFullCounting)
		{
			for (i=0; i<m_N; i++)
				delete m_cacheFullCounting[i];
			delete[] m_cacheFullCounting;
			m_cacheFullCounting = NULL;
		}
	}//end CACHE_AT_OR_NODES

	if (CACHE_AT_AND_NODES)
	{
		if (m_cacheSep)
		{
			for (i=0; i<m_N; i++)
				delete m_cacheSep[i];
			delete[] m_cacheSep;
			m_cacheSep = NULL;
		}
		if (m_cacheSepCounting)
		{
			for (i=0; i<m_N; i++)
				delete m_cacheSepCounting[i];
			delete[] m_cacheSepCounting;
			m_cacheSepCounting = NULL;
		}
	} //end CACHE_AT_AND_NODES
}

// for each variable X in the pseudo tree, maintain a vector of variable for which 
// the cache should be purged for every new assignment of X=x
void CProblemMixed::setAdvancedCacheFlag(int scope)
{
	if (m_advancedCacheFlagAND)
		delete[] m_advancedCacheFlagAND;

	if (m_advancedCacheFlagOR)
		delete[] m_advancedCacheFlagOR;

	m_advancedCacheFlagAND = new variable_v[m_N];
	m_advancedCacheFlagOR = new variable_v[m_N];

	for (int var = 0; var < m_N; var++)
	{
		if (CACHE_AT_AND_NODES)
		{
			assert(m_parentSepSet != NULL);

			if ( m_parentSepSet[var].size() > scope ) 
			{			
				variable_v::reverse_iterator it = m_parentSepSet[var].rbegin();
				it += (scope);
				int flagParent = (*it);
				m_advancedCacheFlagAND[flagParent].push_back(var);
			}
		}


		if (CACHE_AT_OR_NODES)
		{
			assert(m_parentSet != NULL);

			if ( m_parentSet[var].size() > scope )
			{			
				variable_v::reverse_iterator it = m_parentSet[var].rbegin();
				it += (scope);
				m_advancedCacheFlagOR[ (*it) ].push_back(var);
			}			
		}
	}
}


void CProblemMixed::print()
{
	int i;

//	for (i = 0; i < m_N; ++i)
//	{
//		if (isEvidence(i))
//		{
//			cout << "\n - evidence: " << i
//				<< " set to: " << getValue(i);
//		}
//	}
//
	// Ordering:
	cout << "\n Ordering: ";
	for (i = 0; i < m_N; ++i)
	{
		cout << m_ordering[i] << " ";
	}
	cout << endl;

//	m_buckets->print();
//
//	// Legal Tree.
//	m_tree->print();

//	// Descendants in legal tree.
//	cout << "\n Legal Tree Descendants: ";
//	for (i = 0; i < m_N; ++i)
//	{
//		cout << endl;
//		cout << "descendants of " << i <<": " ;
//		copy(m_descendants[i].begin(), m_descendants[i].end(),
//			ostream_iterator<int>(cout, " "));
//	}
//
//
	if (CACHE_AT_AND_NODES)
	{
		// Parent separator set
		cout << "\n Parent Separator Sets: ";

		int maxSepSize = -1;
		for (i = 0; i < m_N; ++i)
		{
//			cout << endl;
//			cout << "parentSepSet[" << i << "]: ";
//			cout << "size = " << m_parentSepSet[i].size() << " --- ";
//			copy(m_parentSepSet[i].begin(), m_parentSepSet[i].end(),
//					ostream_iterator<int>(cout, " "));
			int temp = m_parentSepSet[i].size();
			if (maxSepSize < temp )
				maxSepSize = temp;
//
//			cout << " ***===*** ";
//			cout << "advancedParSepSet[" << i << "]: ";
//			cout << "size = " << m_cacheSep[i]->getArgc() << " --- ";
//			for (int j = 0; j < m_cacheSep[i]->getArgc(); j++)
//				cout << m_cacheSep[i]->getArgv()[j] << " ";
//
//			// check for dead cache
//			CLegalTreeNode** nodes = m_tree->getNodes();
//			int parent;
//			if (nodes[m_position[m_cacheSep[i]->getArgv()[0]]]->parent() == NULL )
//				parent = -1;
//			else
//				parent = nodes[m_position[m_cacheSep[i]->getArgv()[0]]]->parent()->variable();
//			cout << "parent = " << parent;
//
//			for (int k = 0; k < m_N; k++)
//			{
//				variable_v::iterator it = find(m_advancedCacheFlagAND[k].begin(), m_advancedCacheFlagAND[k].end(), i);
//				if ( it != m_advancedCacheFlagAND[k].end() )
//					cout << " flag parent = " << k;			
//			}
//
//
//			cout << "\n      FlagAND:";
//				copy(m_advancedCacheFlagAND[i].begin(), m_advancedCacheFlagAND[i].end(),
//						ostream_iterator<int>(cout, " "));
//			cout << "\n      FlagOR :";
//				copy(m_advancedCacheFlagOR[i].begin(), m_advancedCacheFlagOR[i].end(),
//						ostream_iterator<int>(cout, " "));
//
//
		}

		cout << endl << "max separator size:" << maxSepSize << endl;
	} //end CACHE_AT_AND_NODES

//	if (CACHE_AT_OR_NODES)
//	{
//		// Parent set
//		cout << "\n Parent Sets: ";
//		for (i = 0; i < m_N; ++i)
//		{
//			cout << endl;
//			cout << "parentFullSet[" << i << "]: ";
//			cout << "size = " << m_parentSet[i].size() << " ---";
//			copy(m_parentSet[i].begin(), m_parentSet[i].end(),
//				ostream_iterator<int>(cout, " "));
//		}
//		cout << endl;
//	} // end CACHE_AT_OR_NODES
//
//    if (CACHE_AT_AND_NODES && CACHE_AT_OR_NODES)
//	{
//		for (i = 0; i < m_N; ++i)
//		{
//			cout << endl;
//			cout << "variable " << i << ": ";
//			cout << " full size = " << m_parentSet[i].size() << " --- ";
//			cout << " separator size = " << m_parentSepSet[i].size();		
//		}
//	} // end CACHE_AT_OR_NODES && CACHE_AT_AND_NODES
//
//
//	function_map::iterator it = m_functionsConstraint.begin();
//	for (; it != m_functionsConstraint.end(); ++it)
//	{		
//		CFunction* fun = (*it).second;
//		((CConstraintTable*) fun)->print(true, true);
//	}
//	cout << "this was problem mixed";
//
//	// print the wCutset
//	if (BRUTE_FORCE_W_CUTSET)
//	{
//		for (i = 0; i < m_N; ++i)
//		{
//			cout << endl;
//			cout << "w-cutset[" << i << "]: ";
//			cout << m_wCutset[i];			
//		}	
//	}

}

///////////////
bool CProblemMixed::precompute(int ibound, double& cpuPrecompute)
{
	// Safety checks.
	assert(m_buckets != NULL);

	// Declare variables.
	int k;
	double* bounds = NULL;

	// Run MBE with augmentation.
	double c_start = cpuTime();
	bool bOk = m_buckets->process(MBE_AUGMENT, ibound, k, bounds);
	double c_end = cpuTime();

	// Record pre-computing CPU time.
	cpuPrecompute = (c_end - c_start);

	delete[] bounds;

	return bOk;
}

void CProblemMixed::test()
{
	int k;
	double* bounds = NULL;

	m_buckets->process(MBE_SIMPLE, 2, k, bounds);

	cout << "\n MBE results: [";
	for (int i = 0; i < k; ++i)
	{
		cout << bounds[i] << " ";
	}
	cout << "]" << endl;

	// Free temporary buffers.
	delete[] bounds;
}

// This function initialize the domains.
bool CProblemMixed::initLevelInfo()
{
	// Destroy previous level manager.
	destroyLevelInfo();

	// Create searh level manager.
	m_levelInfo = new SLINFO;
	m_levelInfo->prev = NULL;
	m_levelInfo->level = -1;
	m_levelInfo->status = UNKNOWN;

	int N = m_N;
	int K = m_K;

	// Create domain list/cost.
	m_levelInfo->domainSize = new int[m_N];
	m_levelInfo->domainCost = new double[m_N * m_K];
	m_levelInfo->domainState = new bool[m_N * m_K];
	
	int addr;
	for (int var = 0; var < N; ++var) 
	{
		m_levelInfo->domainSize[var] = getStaticDomainSize(var);
		addr = var*K;
		for (int j = 0; j < m_levelInfo->domainSize[var]; ++j, ++addr) 
		{
			m_levelInfo->domainCost[addr] = 0.0;		// Value cost.
			if (isEvidence(var))
			{
				if (j == getValue(var))
					m_levelInfo->domainState[addr] = true;
				else
					m_levelInfo->domainState[addr] = false;
			}
			else
			{
				m_levelInfo->domainState[addr] = true;		// Value state.
			}
		}
	}

	return true;
}

// This function destroys the search level information.
void CProblemMixed::destroyLevelInfo()
{
	// Safety checks.
	if (m_levelInfo)
	{
		while (NULL != m_levelInfo)
		{
			SLINFO* temp = m_levelInfo;
			m_levelInfo = temp->prev;
			
			temp->destroy();
			delete temp;
		}

		m_levelInfo = NULL;
	}
}

// This function saves the current search level information
// onto domain stack.
void CProblemMixed::pushLevelInfo()
{
	// Safety checks.
	assert(m_levelInfo != NULL);

	// Create next level info structure.
	SLINFO* temp = new SLINFO;
	temp->prev = m_levelInfo;
	temp->level = m_levelInfo->level + 1;
	temp->status = UNKNOWN;
	temp->domainSize = new int[m_N];
	temp->domainCost = new double[m_N * m_K];
    temp->domainState = new bool[m_N * m_K];

	memcpy(temp->domainSize, m_levelInfo->domainSize, m_N * sizeof(int));
	memcpy(temp->domainCost, m_levelInfo->domainCost, m_N * m_K * sizeof(double));
	memcpy(temp->domainState, m_levelInfo->domainState, m_N * m_K * sizeof(bool));

	// Update pointer link.
	m_levelInfo = temp;
}

// This function restores the previous search level info.
void CProblemMixed::popLevelInfo()
{
	// Safety checks.
	assert(m_levelInfo != NULL);

	SLINFO* temp = m_levelInfo->prev;
	
	// Destroy current top of the stack.
	m_levelInfo->destroy();
	delete m_levelInfo;

	// Update pointer link.
	m_levelInfo = temp;
}

// This function saves the current assignment as the best one.
void CProblemMixed::saveCurrentSolution()
{
	// Safety checks.
	assert(m_bestSolution != NULL);
	assert(m_assignment != NULL);
	assert(m_N > 0);

	memcpy(m_bestSolution, m_assignment, m_N * sizeof(int));
}

// This function computes the cost of the current (partial) assignment.
// If the variable in the argument is -1 then we have a full assignment,
// otherwise we compute the utility of the partial assignment up to the
// variable.
double CProblemMixed::utility(int var)
{
	// Safety checks.
	assert(m_assignment != NULL);
	assert(m_N > 0);

	double util = 1.0;
	
	if (-1 == var)
	{
		// Evaluate a full assignment.
		function_map::iterator it = m_functions.begin();
		for (; it != m_functions.end(); ++it)
		{
			double p;
			CFunction* fun = (*it).second;

			fun->getCurrentValue(&p);

			util *= p;
		}
	}
	else
	{
		// Get current search level.
		int level = m_position[var];

		// Evaluate a partial assignment.
		for (int pos = 0; pos <= level; ++pos)
		{
			// Process original functions only.
			CBucket* bucket = m_buckets->getBucketAt(pos);
			function_v& funs = bucket->functions();
			function_v::iterator it = funs.begin();
			for (; it != funs.end(); ++it)
			{
				CFunction* fun = (*it);
				if (!fun->isOriginal())
					continue;	// skip non-original functions.

				double p;
				assert(-1 != m_assignment[var]);

				fun->getCurrentValue(&p);

				util *= p;
			}
		}
	}

	return util;
}

// This function evalues the heuristic value costs at current search level.
void CProblemMixed::evalStaticHeuristic(int var)
{
	// Safety checks.
	assert(m_ordering != NULL);
	assert(m_position != NULL);
	assert((var >= 0) && (var < m_N));

	// Get current search tree level.
	int level = m_levelInfo->level;
	int k = getStaticDomainSize(var);

	double* domainCost = m_levelInfo->domainCost;
	assert(domainCost != 0);

	// Current variable's address index.
	int baseAddr = var * m_K;
	
	// Heuristic's f,g,h components.
	double* g = new double[k];
	double* h = new double[k];
	double* f = new double[k];
	memset(h, 0, k * sizeof(double));
	memset(f, 0, k * sizeof(double));
	memset(g, 0, k * sizeof(double));

	// Get the current level's bucket.
	CBucket* bucket = m_buckets->getBucketAt(level);

	// Compute "h" component.
	bucket->evalStaticHeuristic(k, h);

	// Check for evidence and keep the value.
	if (isEvidence(var))
	{
		// Get evidence value.
		int i = getValue(var);

		// Update the address.
		int addr = baseAddr + i;
	
		// Compute "g" component.
		g[i] = utility(var);

		// Heuristic evaluation function.
		f[i] = g[i] * h[i];

		// Update current value's cost.
		domainCost[addr] = f[i];
	}
	else
	{
		// Backup current assignment.
		backupAssignment();

		// Compute "g" component.
		for (int i = 0; i < k; ++i)
		{
			// Update the address.
			int addr = baseAddr + i;

			setValue(var, i);
			g[i] = utility(var);

			// Heuristic evaluation function.
			f[i] = g[i] * h[i];

			// Update current value's cost.
			domainCost[addr] = f[i];
		}

		// Restore current assignment.
		restoreAssignment();
	}

	delete[] g;
	delete[] h;
	delete[] f;
}

void CProblemMixed::evalDynamicHeuristic(int var, int ibound)
{
	// Safety checks.
	assert(m_buckets != NULL);


	// Declare variables.
	int k;
	double* h = NULL;

	// Get current variable's level (position).
	int level = m_position[var];

	// Run MBE with current variable as first in the ordering.
	m_buckets->process(MBE_PARTIAL, ibound, k, h, level);

	// Update current value costs.
	double* domainCost = m_levelInfo->domainCost;
	assert(domainCost != 0);

	// Current variable's address index.
	int baseAddr = var * m_K;
	double* g = new double[k];
	double* f = new double[k];
	memset(f, 0, k * sizeof(double));
	memset(g, 0, k * sizeof(double));
	
	// Check for evidence and keep the value.
	if (isEvidence(var))
	{
		// Get evidence value.
		int i = getValue(var);

		// Update the address.
		int addr = baseAddr + i;
	
		// Compute "g" component. Need to refer to the previous variable.
		int prev = (level > 0) ? m_ordering[level - 1] : var;
		g[i] = (level > 0) ? utility(prev) : 1.0;

		// Heuristic evaluation function.
		f[i] = g[i] * h[i];

		// Update current value's cost.
		domainCost[addr] = f[i];
	}
	else
	{
		// Backup current assignment.
		backupAssignment();

		for (int i = 0; i < k; ++i)
		{
			// Update the address.
			int addr = baseAddr + i;

			// Compute "g" component. Need to refer to previous variable.
			int prev = (level > 0) ? m_ordering[level - 1] : var;
			g[i] = (level > 0) ? utility(prev) : 1.0;

			// Heuristic evaluation function.
			f[i] = g[i] * h[i];

			// Update current value's cost.
			domainCost[addr] = f[i];
		}

		// Restore current assignment.
		restoreAssignment();
	}

	delete[] h;
	delete[] g;
	delete[] f;
}

// This function selects the next variable for a static ordering.
int CProblemMixed::selectVariable(int svType)
{
	// Safety checks.
	assert(m_assignment != NULL);
	assert(m_ordering != NULL);
	assert(m_position != NULL);

	int var;

	switch (svType)
	{
	case STATIC:
		{
			// Get current search tree level.
			int level = m_levelInfo->level;

			var = m_ordering[level + 1];

			break;
		}
	case DYNAMIC:
		{
			break;
		}
	}

	return var;
}

// This function selects the next value for the current variable.
int CProblemMixed::selectValue(int var, int& val, bool prune)
{
	// Safety checks.
	assert(m_assignment != NULL);

	int result = SV_FAILURE;		// Result of the function.

	int* domainSize = m_levelInfo->domainSize;
	double* domainCost = m_levelInfo->domainCost;
	bool* domainState = m_levelInfo->domainState;
	
	assert(domainSize != NULL);
	assert(domainCost != NULL);
	assert(domainState != NULL);

	// Select next best value.
	int bestValue = -1;
	double bestCost = -1.0;
	int baseAddr = var * m_K;

	int K = getStaticDomainSize(var);
	for (int k = 0; k < K; ++k)
	{
		int addr = baseAddr + k;
		if (false == domainState[addr])
			continue;	// skip already visited values.

		int value = k;
		double cost = domainCost[addr];

		if (prune)
		{
			if (cost <= m_bestSolutionCost)
			{
				// Current assignment cannot be extended 
				// to a better solution. Disable the value.
				domainState[addr] = false;
				continue;	
			}
			else
			{
				// Current assignment is a valid one.
				if (cost > bestCost)
				{
					bestValue = value;
					bestCost = cost;
				}
			}
		}
		else
		{
			bestValue = value;
			break;
		}
	}

	// Next value selection.
	val = bestValue;

	if (-1 == val)
	{
		// No value found.
		result = SV_NOVALUE;
	}
	else
	{
		// Found a value. Mark it as used.
		domainState[baseAddr + val] = false;

		result = SV_SUCCESS;
	}

	return result;
}

// This function runs a Branch and Bound search algorithm. It is guided
// by a set of pre-computed heuristics for each variable-value pair. Assumes
// static variable ordering and dynamic value ordering.
//int CProblemMixed::execBBMBs(int ibound, long tmLimit, double& tmCpuSearch,
//	long& backtracks, long& expansions, bool silent)
//{
//	// Safety checks.
//	assert(m_N > 0);
//	assert(m_K > 0);
////robert	assert(m_C > 0);
//
//	// DFS variable stack.
//	stack<int> dfsStack;
//
//	// General variables.
//	int result;						// Output code.
//	double bound;					// Current lower bound.
//	clock_t c_start;				// Start timer.
//	clock_t c_end;					// End timer.
//	char tmpbuf[1024];				// Buffer.
//
//	// Intermediate results.
//	long _timeNextSlice = 1024;		// Next time point for intermediate time (every 1 sec).
//	long _timeLastCheck = 0 ;		// Time when the time-check was done last time.
//	time_t time_counter ;
//	unsigned long current_clock ;
//	double tmSearch, tmPrecompute;
//
//	// Set time zone from TZ environment variable.
//	//_tzset();
//
//	// Prologue.
//	sprintf(tmpbuf, "\nBBMB-static started ...");
//	cout << tmpbuf;
//	if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//	sprintf(tmpbuf, "\n   i-bound: %d", ibound);
//	cout << tmpbuf;
//	if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//	sprintf(tmpbuf, "\n   time limit: %d", tmLimit);
//	cout << tmpbuf;
//	if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//	
//	// Precompute heuristics.
//	if (!precompute(ibound, tmPrecompute))
//	{
//		cout << "\n   FAILURE: Pre-Processing phase failed.";
//		cout << "\nBBMB-static aborted.";
//		result = S_FAILURE;
//		goto done;
//	}
//
//	// Output pre-processing statistics.
//	sprintf(tmpbuf, "\n   PP CPU time: %7g", tmPrecompute);
//	cout << tmpbuf;
//	if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//	sprintf(tmpbuf, "\n   silent mode: %s", ((!silent) ? "true" : "false"));
//	cout << tmpbuf;
//	if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//	if (!silent)
//	{
//		cout << "\n   DETAILS";
//		if (m_outputFile) fprintf(m_outputFile, "\n   DETAILS");
//	}
//
//	// Start the timer (_clock)
//	c_start = clock();
//
//	// Start the local clock
//	startClock(&time_counter) ;
//
//	// Initialize search.
//	backtracks = 0;
//	expansions = 0;
//
//	// Init current lower bound.
//	m_bestSolutionCost = 0;
//	if (m_bestSolution)
//	{
//		delete[] m_bestSolution;
//		m_bestSolution = NULL;
//	}
//	m_bestSolution = new int[m_N];
//	memset(m_bestSolution, -1, m_N * sizeof(int));
//
//	// Init search tree level info.
//	initLevelInfo();
//
//	// Set first variable.
//	int var0; var0 = m_ordering[0];
//	
//	// Initialize the variable stack.
//	dfsStack.push(var0);
//
//	// Push search level for variable.
//	pushLevelInfo();
//	evalStaticHeuristic(var0);
//
//	result = S_FAILURE;
//
//	// Start search.
//	while (true)
//	{
//		// EXPAND step : 
//		//   1a) pick value for the current variable.
//		//   1b) if no values left, then 
//		//	      1b-1) roll back variable-stack.
//		//        1b-2) BACKTRACK.
//		//    2) pick next var (static order) and add to variable-stack (and make current).
//		//    3) EXPAND again
//
//		// BACKTRACK step :
//		//    1) restore domain of current variable.
//		//    2) EXPAND.
//
//expand:
//
//		// Check for silent run.
//		if (!silent)
//		{
//			current_clock = stopClock(0, &time_counter);
//			if (current_clock >= _timeNextSlice && 
//				_timeLastCheck < _timeNextSlice) 
//			{
//				// print intemediate results
//				double current_time = ((double) current_clock) / 1000.0 ;
//				
//				sprintf(tmpbuf, "\n    -int: time %7g, nodes %d, bkts %d, cost %7g",
//					current_time, expansions, backtracks, m_bestSolutionCost);
//
//				cout << tmpbuf;
//				if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//				
//				// Compute next time check bound as next 1024 multiple time point from the current one.
//				int _timeFraction = current_clock & 7 ;
//				_timeNextSlice = current_clock + 1024 - _timeFraction ;
//
//				// Save current time for next round.
//				_timeLastCheck = current_clock ;
//			}
//		}
//
//		// Check for time limit violation
//		if ((tmLimit > 0) && (clock() - c_start >= tmLimit))
//		{
//			result = S_TIMEOUT;		// timeout.
//
//			break;
//		}
//
//		// Node expansion.
//		int var = dfsStack.top();
//		
//		// Value selection.
//		int val;
//		switch (selectValue(var, val))
//		{
//		case SV_FAILURE:
//			{
//				result = S_FAILURE;
//				goto done;
//			}
//		case SV_NOVALUE:
//			{
//				// Remove current value for variable.
//				if (!isEvidence(var))
//					resetValue(var);
//
//				// Pop variable stack
//				dfsStack.pop();
//				++backtracks;
//
//				// Check for empty stack (done).
//				if (0 == dfsStack.size())
//				{
//					result = S_SUCCESS;
//					goto done;	// Complete search.
//				}
//
//				goto backtrack;
//
//				break;
//			}
//		case SV_SUCCESS:
//			{
//				++expansions;
//				
//				// Set value for top variable.
//				if (!isEvidence(var))
//					setValue(var, val);
//
//				// Remove value from top variable's domain.
//
//				// Check for solution.
//				if (m_N == dfsStack.size())
//				{
//					// Compute solution cost (lower bound).
//					bound = utility();
//					
//					if (bound > m_bestSolutionCost)
//					{
//						m_bestSolutionCost = bound;
//						saveCurrentSolution();
//					}
//					
//					goto expand;
//				}
//
//				break;
//			}
//		}
//
//		// Select next variable.
//		int nextVar;
//		nextVar = selectVariable(STATIC);
//
//		// Add a new variable.
//		dfsStack.push(nextVar);
//
//		// Push search level for variable.
//		pushLevelInfo();
//		evalStaticHeuristic(nextVar);
//
//		continue;
//
//backtrack:
//		
//		// Pop search level info from stack.
//		popLevelInfo();
//	}
//
//done:
//
//	// Done.
//	c_end = clock();
//	current_clock = stopClock(0, &time_counter) ;
//	double tm1 = (double)current_clock / 1000.0 ;
//	double tm2 = (double)(c_end - c_start) / 1000.0;
//	tmSearch = tm2;
//
//	tmCpuSearch = tmSearch;
//
//	// Output solution.
//	switch (result)
//	{
//	case S_SUCCESS:
//		{
//			sprintf(tmpbuf, "\n   out of time: false");
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n   optimal solution cost: %7g", m_bestSolutionCost);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n\n   CPU time (search): %7g", tmSearch);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n   CPU time (preproc.): %7g", tmPrecompute);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n   CPU time (total): %7g", tmSearch + tmPrecompute);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n   number of expansions: %d", expansions);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n   number of backtracks: %d", backtracks);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			cout << "\nEND.\n";
//			if (m_outputFile) fprintf(m_outputFile, "\nEND.\n");
//
//			break;
//		}
//	case S_TIMEOUT:
//		{
//			sprintf(tmpbuf, "\n   out of time: true");
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n   best solution cost: %7g", m_bestSolutionCost);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n\n   CPU time (search): %7g", tmSearch);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n   CPU time (preproc.): %7g", tmPrecompute);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n   CPU time (total): %7g", tmSearch + tmPrecompute);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n   number of expansions: %d", expansions);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n   number of backtracks: %d", backtracks);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			cout << "\nEND.\n";
//			if (m_outputFile) fprintf(m_outputFile, "\nEND.\n");
//
//			break;
//		}
//	};
//
//	// Clear memory.
//	destroyLevelInfo();
//	m_buckets->clean();
//	for (int v = 0; v < m_N; ++v)
//	{
//		if (isEvidence(v)) continue;
//
//		resetValue(v);
//	}
//
//	return result;
//}
//

// This function runs a Branch and Bound search algorithm. It is guided
// by a set of heuristics that are dynamically computed using MBE. Assumes
// static variable ordering and dynamic value ordering.
//int CProblemMixed::execBBMBd(int ibound, long tmLimit, double& tmCpuSearch,
//	long& backtracks, long& expansions, bool silent)
//{
//	// Safety checks.
//	assert(m_N > 0);
//	assert(m_K > 0);
////robert	assert(m_C > 0);
//
//	// DFS variable stack.
//	stack<int> dfsStack;
//
//	// General variables.
//	int result;						// Output code.
//	double bound;					// Current lower bound.
//	clock_t c_start;				// Start timer.
//	clock_t c_end;					// End timer.
//	char tmpbuf[1024];				// Buffer.
//
//	// Intermediate results.
//	long _timeNextSlice = 1024;		// Next time point for intermediate time (every 1 sec).
//	long _timeLastCheck = 0 ;		// Time when the time-check was done last time.
//	time_t time_counter ;
//	unsigned long current_clock ;
//	double tmSearch;
//
//	// Set time zone from TZ environment variable.
//	//_tzset();
//
//	// Prologue.
//	sprintf(tmpbuf, "\nBBMB-dynamic started ...");
//	cout << tmpbuf;
//	if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//	sprintf(tmpbuf, "\n   i-bound: %d", ibound);
//	cout << tmpbuf;
//	if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//	sprintf(tmpbuf, "\n   time limit: %d", tmLimit);
//	cout << tmpbuf;
//	if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//	
//	sprintf(tmpbuf, "\n   silent mode: %s", ((!silent) ? "true" : "false"));
//	cout << tmpbuf;
//	if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//	if (!silent)
//	{
//		cout << "\n   DETAILS";
//		if (m_outputFile) fprintf(m_outputFile, "\n   DETAILS");
//	}
//
//	// Start the timer (_clock)
//	c_start = clock();
//
//	// Start the local clock
//	startClock(&time_counter) ;
//
//	// Initialize search.
//	backtracks = 0;
//	expansions = 0;
//
//	// Init current lower bound.
//	m_bestSolutionCost = 0;
//	if (m_bestSolution)
//	{
//		delete[] m_bestSolution;
//		m_bestSolution = NULL;
//	}
//	m_bestSolution = new int[m_N];
//	memset(m_bestSolution, -1, m_N * sizeof(int));
//
//	// Init search tree level info.
//	initLevelInfo();
//
//	// Set first variable.
//	int var0; var0 = m_ordering[0];
//	
//	// Initialize the variable stack.
//	dfsStack.push(var0);
//
//	// Push search level for variable.
//	pushLevelInfo();
//	evalDynamicHeuristic(var0, ibound);
//
//	result = S_FAILURE;
//
//	// Start search.
//	while (true)
//	{
//		// EXPAND step : 
//		//   1a) pick value for the current variable.
//		//   1b) if no values left, then 
//		//	      1b-1) roll back variable-stack.
//		//        1b-2) BACKTRACK.
//		//    2) pick next var (static order) and add to variable-stack (and make current).
//		//    3) EXPAND again
//
//		// BACKTRACK step :
//		//    1) restore domain of current variable.
//		//    2) EXPAND.
//
//expand:
//
//		// Check for silent run.
//		if (!silent)
//		{
//			current_clock = stopClock(0, &time_counter);
//			if (current_clock >= _timeNextSlice && 
//				_timeLastCheck < _timeNextSlice) 
//			{
//				// print intemediate results
//				double current_time = ((double) current_clock) / 1000.0 ;
//				
//				sprintf(tmpbuf, "\n    -int: time %7g, nodes %d, bkts %d, cost %7g",
//					current_time, expansions, backtracks, m_bestSolutionCost);
//
//				cout << tmpbuf;
//				if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//				
//				// Compute next time check bound as next 1024 multiple time point from the current one.
//				int _timeFraction = current_clock & 7 ;
//				_timeNextSlice = current_clock + 1024 - _timeFraction ;
//
//				// Save current time for next round.
//				_timeLastCheck = current_clock ;
//			}
//		}
//
//		// Check for time limit violation
//		if ((tmLimit > 0) && (clock() - c_start >= tmLimit))
//		{
//			result = S_TIMEOUT;		// timeout.
//
//			break;
//		}
//
//		// Node expansion.
//		int var = dfsStack.top();
//
//		if (var == 0)
//			int b = 0;
//
//		// Value selection.
//		int val;
//		switch (selectValue(var, val))
//		{
//		case SV_FAILURE:
//			{
//				result = S_FAILURE;
//				goto done;
//			}
//		case SV_NOVALUE:
//			{
//				// Remove current value for variable.
//				if (!isEvidence(var))
//					resetValue(var);
//
//				// Pop variable stack
//				dfsStack.pop();
//				++backtracks;
//
//				// Check for empty stack (done).
//				if (0 == dfsStack.size())
//				{
//					result = S_SUCCESS;
//					goto done;	// Complete search.
//				}
//
//				goto backtrack;
//
//				break;
//			}
//		case SV_SUCCESS:
//			{
//				++expansions;
//				
//				// Set value for top variable.
//				if (!isEvidence(var))
//					setValue(var, val);
//
//				// Check for solution.
//				if (m_N == dfsStack.size())
//				{
//					// Compute solution cost (lower bound).
//					bound = utility();
//					
//					if (bound > m_bestSolutionCost)
//					{
//						m_bestSolutionCost = bound;
//						saveCurrentSolution();
//					}
//					
//					goto expand;
//				}
//
//				break;
//			}
//		}
//
//		// Select next variable.
//		int nextVar;
//		nextVar = selectVariable(STATIC);
//
//		// Add a new variable.
//		dfsStack.push(nextVar);
//
//		// Push search level for variable.
//		pushLevelInfo();
//		evalDynamicHeuristic(nextVar, ibound);
//
//		continue;
//
//backtrack:
//		
//		// Pop search level info from stack.
//		popLevelInfo();
//	}
//
//done:
//
//	// Done.
//	c_end = clock();
//	current_clock = stopClock(0, &time_counter) ;
//	double tm1 = (double)current_clock / 1000.0 ;
//	double tm2 = (double)(c_end - c_start) / 1000.0;
//	tmSearch = tm2;
//
//	tmCpuSearch = tmSearch;
//
//	// Output solution.
//	switch (result)
//	{
//	case S_SUCCESS:
//		{
//			sprintf(tmpbuf, "\n   out of time: false");
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n   optimal solution cost: %7g", m_bestSolutionCost);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n\n   CPU time (search): %7g", tmSearch);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n   number of expansions: %d", expansions);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n   number of backtracks: %d", backtracks);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			cout << "\nEND.\n";
//			if (m_outputFile) fprintf(m_outputFile, "\nEND.\n");
//
//			break;
//		}
//	case S_TIMEOUT:
//		{
//			sprintf(tmpbuf, "\n   out of time: true");
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n   best solution cost: %7g", m_bestSolutionCost);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n\n   CPU time (search): %7g", tmSearch);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n   number of expansions: %d", expansions);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n   number of backtracks: %d", backtracks);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			cout << "\nEND.\n";
//			if (m_outputFile) fprintf(m_outputFile, "\nEND.\n");
//
//			break;
//		}
//	};
//
//	// Clear memory.
//	destroyLevelInfo();
//	m_buckets->clean();
//	for (int v = 0; v < m_N; ++v)
//	{
//		if (isEvidence(v)) continue;
//
//		resetValue(v);
//	}
//
//	return result;
//}
//

// This function runs a Backtracking search algorithm. It has no heuristic
// guidance and assumes static variable ordering. It is used for comparison
// with the more advanced AND-OR search (without pruning) algorihm.
int CProblemMixed::execBT(long tmLimit, double& tmCpuSearch,
	long& backtracks, long& expansions, bool silent)
{
	// Safety checks.
	assert(m_N > 0);
	assert(m_K > 0);
//robert	assert(m_C > 0);

	// DFS variable stack.
	stack<int> dfsStack;

	// General variables.
	int result;						// Output code.
	double bound;					// Current lower bound.
	double c_start;				// Start timer.
	double c_end;					// End timer.
	char tmpbuf[1024];				// Buffer.

	// Intermediate results.
	long _timeNextSlice = 1024;		// Next time point for intermediate time (every 1 sec).
	long _timeLastCheck = 0 ;		// Time when the time-check was done last time.
	time_t time_counter ;
	unsigned long current_clock ;
	double tmSearch;

	// Set time zone from TZ environment variable.
	//_tzset();

	// Prologue.
	sprintf(tmpbuf, "\nBT started ...");
	cout << tmpbuf;
	if (m_outputFile) fprintf(m_outputFile, tmpbuf);

	sprintf(tmpbuf, "\n   time limit: %d", tmLimit);
	cout << tmpbuf;
	if (m_outputFile) fprintf(m_outputFile, tmpbuf);
	
	sprintf(tmpbuf, "\n   silent mode: %s", ((!silent) ? "true" : "false"));
	cout << tmpbuf;
	if (m_outputFile) fprintf(m_outputFile, tmpbuf);

	if (!silent)
	{
		cout << "\n   DETAILS";
		if (m_outputFile) fprintf(m_outputFile, "\n   DETAILS");
	}

	// Start the timer (_clock)
	c_start = cpuTime();

	//// Start the local clock
	//startClock(&time_counter) ;

	// Initialize search.
	backtracks = 0;
	expansions = 0;

	// Init current lower bound.
	m_bestSolutionCost = 0;
	if (m_bestSolution)
	{
		delete[] m_bestSolution;
		m_bestSolution = NULL;
	}
	m_bestSolution = new int[m_N];
	memset(m_bestSolution, -1, m_N * sizeof(int));

	// Init search tree level info.
	initLevelInfo();

	// Set first variable.
	int var0; var0 = m_ordering[0];
	
	// Initialize the variable stack.
	dfsStack.push(var0);

	// Push search level for variable.
	pushLevelInfo();

	result = S_FAILURE;

	// Start search.
	while (true)
	{
		// EXPAND step : 
		//   1a) pick value for the current variable.
		//   1b) if no values left, then 
		//	      1b-1) roll back variable-stack.
		//        1b-2) BACKTRACK.
		//    2) pick next var (static order) and add to variable-stack (and make current).
		//    3) EXPAND again

		// BACKTRACK step :
		//    1) restore domain of current variable.
		//    2) EXPAND.

expand:

		// Check for silent run.
		//if (!silent)
		//{
		//	current_clock = stopClock(0, &time_counter);
		//	if (current_clock >= _timeNextSlice && 
		//		_timeLastCheck < _timeNextSlice) 
		//	{
		//		// print intemediate results
		//		double current_time = ((double) current_clock) / 1000.0 ;
		//		
		//		sprintf(tmpbuf, "\n    -int: time %7g, nodes %d, bkts %d, cost %7g",
		//			current_time, expansions, backtracks, m_bestSolutionCost);

		//		cout << tmpbuf;
		//		if (m_outputFile) fprintf(m_outputFile, tmpbuf);
		//		
		//		// Compute next time check bound as next 1024 multiple time point from the current one.
		//		int _timeFraction = current_clock & 7 ;
		//		_timeNextSlice = current_clock + 1024 - _timeFraction ;

		//		// Save current time for next round.
		//		_timeLastCheck = current_clock ;
		//	}
		//}

		// Check for time limit violation
		if ((tmLimit > 0) && (cpuTime() - c_start >= tmLimit))
		{
			result = S_TIMEOUT;		// timeout.

			break;
		}

		// Node expansion.
		int var = dfsStack.top();
		
		// Value selection.
		int val;
		switch (selectValue(var, val, false))
		{
		case SV_FAILURE:
			{
				result = S_FAILURE;
				goto done;
			}
		case SV_NOVALUE:
			{
				// Remove current value for variable.
				if (!isEvidence(var))
					resetValue(var);

				// Pop variable stack
				dfsStack.pop();
				++backtracks;

				// Check for empty stack (done).
				if (0 == dfsStack.size())
				{
					result = S_SUCCESS;
					goto done;	// Complete search.
				}

				goto backtrack;

				break;
			}
		case SV_SUCCESS:
			{
				++expansions;
				
				// Set value for top variable.
				if (!isEvidence(var))
					setValue(var, val);

				// Check for solution.
				if (m_N == dfsStack.size())
				{
					// Compute solution cost (lower bound).
					bound = utility();
					
					if (bound > m_bestSolutionCost)
					{
						m_bestSolutionCost = bound;
						saveCurrentSolution();
					}
					
					goto expand;
				}

				break;
			}
		}

		// Select next variable.
		int nextVar;
		nextVar = selectVariable(STATIC);

		// Add a new variable.
		dfsStack.push(nextVar);

		// Push search level for variable.
		pushLevelInfo();

		continue;

backtrack:
		
		// Pop search level info from stack.
		popLevelInfo();
	}

done:

	// Done.
	c_end = cpuTime();
	// current_clock = stopClock(0, &time_counter) ;

	// double tm1 = (double)current_clock / 1000.0 ;
	double tm2 = (double)(c_end - c_start);
	tmSearch = tm2;

	tmCpuSearch = tmSearch;

	// Output solution.
	switch (result)
	{
	case S_SUCCESS:
		{
			sprintf(tmpbuf, "\n   out of time: false");
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			sprintf(tmpbuf, "\n   optimal solution cost: %7g", m_bestSolutionCost);
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			sprintf(tmpbuf, "\n\n   CPU time (search): %7g", tmSearch);
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			sprintf(tmpbuf, "\n   number of expansions: %d", expansions);
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			sprintf(tmpbuf, "\n   number of backtracks: %d", backtracks);
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			cout << "\nEND.\n";
			if (m_outputFile) fprintf(m_outputFile, "\nEND.\n");

			break;
		}
	case S_TIMEOUT:
		{
			sprintf(tmpbuf, "\n   out of time: true");
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			sprintf(tmpbuf, "\n   best solution cost: %7g", m_bestSolutionCost);
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			sprintf(tmpbuf, "\n\n   CPU time (search): %7g", tmSearch);
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			sprintf(tmpbuf, "\n   number of expansions: %d", expansions);
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			sprintf(tmpbuf, "\n   number of backtracks: %d", backtracks);
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			cout << "\nEND.\n";
			if (m_outputFile) fprintf(m_outputFile, "\nEND.\n");

			break;
		}
	};

	// Clear memory.
	destroyLevelInfo();

	return result;
}


void CProblemMixed::removeNode(list<CAONode*>& l, CAONode* node)
{
	list<CAONode*>::iterator it = l.begin();
	while (it != l.end())
	{
		CAONode* tmp = (*it);
		if (tmp == node)
		{
			l.erase(it);
			break;
		}

		++it;
	}
}

bool CProblemMixed::forwardPrune(CAONode* node, double estimate)
{
	// Safety checks.
	assert(node->type() == AND);

	CAONode* parentOR = node->parent();			// OR parent (consistent with the legal tree).
	assert(parentOR != NULL);

	if (parentOR->isUpdated())
	{
		double pseudoMax = estimate * node->getG();
		double gValueOR = parentOR->getG();

		if (pseudoMax <= gValueOR)
			return true;
		else
			return false;
	}

	return false;
}


// This function expands and AND-OR search tree node.
//int CProblemMixed::expand(CAONode* node, stack<CAONode*>& succ, bool heuristic, int ibound, int usePruning)
//{
//	// Safety checks.
//	assert(node != NULL);
//	assert(m_tree != NULL);
//
//	int result = E_FAILURE;
////	char tmpbuf[1024];
//
//	switch (node->type())
//	{
//	case AND:	// Expand an AND node (value).
//		{
//#ifdef DEBUG_TIME
//			c_expand_AND_start = clock();
//#endif
//
//			// Generate OR successors.
//			CAONode* parent = node->parent();
//			assert(parent != NULL);
//			assert(parent->type() == OR);
//
//			// Before expanding the node, must compute its g-value.
//			int var = parent->label();
//			int val = node->label();
//			
//			// Set current value.
//			setValue(var, val);
//
//			// purge cache for advanced caching
//			if (CUTSET_CACHE)
//			{
//				if (CACHE_AT_AND_NODES)
//				{
//					variable_v::iterator it = m_advancedCacheFlagAND[var].begin();
//					for ( ; it != m_advancedCacheFlagAND[var].end(); it++ )
//					{
//						m_cacheSep[(*it)]->purge();									
//						if (COUNTING)
//							m_cacheSepCounting[(*it)]->purge();
//					}
//				}
//				
//				if (CACHE_AT_OR_NODES)					
//				{
//					variable_v::iterator it = m_advancedCacheFlagOR[var].begin();
//					for ( ; it != m_advancedCacheFlagOR[var].end(); it++ )
//					{
//						m_cacheFull[(*it)]->purge();									
//						if (COUNTING) 
//							m_cacheFullCounting[(*it)]->purge();
//					}					
//				}
//			}
//
//			int pos = m_position[var];
//			CBucket* bucket = m_buckets->getBucketAt(pos);
//			//assert(bucket != NULL);
//			
//			// Multiply all original functions in the bucket.
//			// here g is 1 if the bucket is empty
//			
//			// robert: evaluate buckets for probabilities, initialize node with bucket result
//
//			double el = bucket->multiply();	
//			node->setL(el);			
//			node->setG(1.0);
//			
//			// for counting solutions initialize with 1.0
//			if(COUNTING)
//				node->setGCount(1.0);
//
//
////			sprintf(tmpbuf, "\nLOG:  -[AND] initial g-value is set to %7g", g);
////			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			// Check for zero l-value. Prune under it.
/////*
//			if (0.0 == node->getL())
//			{
//				result = S_SUCCESS;
//				zeroDeadEnds += 1;
//				break;
//			}	
////*/
//
//			// Get variable's children from the legal tree.
//			variable_v children;
//			m_tree->getChildren(var, children);
//
//			if (!heuristic)
//			{
////				sprintf(tmpbuf, "\nLOG:  -[AND] children are: ");
//				variable_v::iterator it = children.begin();
//				for (; it != children.end(); ++it)
//				{
//					int ch = (*it);
//					CAONode* child = new CAONode(OR, ch);
//					assert(child != NULL);
//
//					child->setParent(node);
//					node->addChild(child);
//
//					succ.push(child);
//
////					char tmp[10];
////					sprintf(tmp, "%d ", ch);
////					strcat(tmpbuf, tmp);
//				}
//			}
//			else
//			{
//				double pseudoG = 1.0;
////				sprintf(tmpbuf, "\nLOG:  -[AND] children are: ");
//				variable_v::iterator it = children.begin();
//				for (; it != children.end(); ++it)
//				{
//					int ch = (*it);
//					CAONode* child = new CAONode(OR, ch);
//					assert(child != NULL);
//
//					child->setParent(node);
//					node->addChild(child);
//
//					succ.push(child);
//
////					char tmp[10];
////					sprintf(tmp, "%d ", ch);
////					strcat(tmpbuf, tmp);
//
//					// Estimate value costs for the OR child.
//					double maxCost;
//					child->buildValueCosts(this, ibound, maxCost); 
//
//					pseudoG *= maxCost;
//				}
//
//				// Check if we need to prune the current AND node.
//				if (forwardPrune(node, pseudoG))
//				{
//					node->setG(-1.0);	// make sure that is not going to be maximized by OR parent.
//					node->clear();
//					while (!succ.empty())
//					{
//						CAONode* tmp = succ.top();
//						succ.pop();
//						delete tmp;
//					}
//
////					sprintf(tmpbuf, "\nLOG:  -[AND] children are: pruned");
//				}
//			}
//
//			result = E_SUCCESS;
//
////			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			// Free temporary buffers.
//			children.clear();
//
//#ifdef DEBUG_TIME
//			c_expand_AND_end = clock();
//			tmExpandAND += (double)(c_expand_AND_end - c_expand_AND_start) / 1000.0;
//#endif
//
//			break;
//		}
//	case OR:	// Expand an OR node (variable).
//		{
//#ifdef DEBUG_TIME
//			c_expand_OR_start = clock();
//#endif
//
//			// Check for heuristic computation (full AO tree expansion).
//			if (!heuristic)
//			{
//				// Generate AND successors.
//				int var = node->label();
//				vector<int> values;
//
//				// Check for evidence variable.
//				if (isEvidence(var))
//				{
//					int val = getValue(var);
//					values.push_back(val);
//				}
//				else
//				{
//					for (int val = 0; val < getStaticDomainSize(var); ++val)
//					{
//						values.push_back(val);
//					}
//				}
//
//				// Create AND successors for each value.
//
//				// Here we check different consistency levels:
//
//				vector<int>::iterator it = values.begin();
//				for (; it != values.end(); ++it)
//				{
//					int val = (*it);
//
//					// Check consistency of the current value.
//					// If no values have been found consistent,
//					// then g remains 0 as is set by default and 
//					// succ is empty, which will trigger propagation up
//
//					// if(!usePruning || forwardChecking(var, val))
//					
//					bool cons = false;
//					switch (usePruning)
//					{
//					case PRUNING_NO:
//						{
//							cons = true;
//							break;
//						}
//					case PRUNING_CONSTRAINTS_ONLY:
//						{
//							cons = consistent(var,val);
//							break;						
//						}
//					case PRUNING_FORWARD_CHECKING:
//						{
//							cons = (consistent(var, val) && forwardChecking(var, val));
//							break;
//						}
//					case PRUNING_FC_PROJECTION:
//						{
//							cons = (consistent(var, val) && forwardChecking(var, val) && forwardCheckingByProjection(var, val));
//							break;
//						}
//					case PRUNING_FULL_LOOK_AHEAD:	
//						break;
//					}
//
//
//					if (cons)
//					{
//						CAONode* child = new CAONode(AND, val);
//						assert(child != NULL);
//
//						// Set up links.
//						child->setParent(node);
//						node->addChild(child);
//
//						succ.push(child);
//					}
//					else 
//						deadEnds += 1;
//				}
//
//				result = E_SUCCESS;
//
//				// Free temporary buffers.
//				values.clear();
//			}
//
//			/*else
//			{
//				// Generate AND successors.
//				int var = node->label();
//				vector<VALUECOST*>& estimates = node->getValueCosts();				
//
//				// Create AND successors for each value, in decreasing order of its estimate.
//				vector<VALUECOST*>::iterator it = estimates.begin();
//				for (; it != estimates.end(); ++it)
//				{
//					int val = (*it)->value;
//					double cost = (*it)->cost;
//
//					CAONode* child = new CAONode(AND, val);
//					assert(child != NULL);
//
//					// Set up links.
//					child->setCost(cost);
//					child->setParent(node);
//					node->addChild(child);
//
//					succ.push(child);
//
////					sprintf(tmpbuf, "\nLOG:  -[OR] child %d - %7g", val, cost);
////					if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//				}
//
//				result = E_SUCCESS;
//			}
//			*/		
//#ifdef DEBUG_TIME
//			c_expand_OR_end = clock();
//			tmExpandOR += (double)(c_expand_OR_end - c_expand_OR_start) / 1000.0;
//#endif
//
//			break;
//		}
//	}
//
//	return result;
//}
//


// This function expands and AND-OR search tree node.
int CProblemMixed::expand(CAONode* node, stack<CAONode*>& succ, bool heuristic, int ibound, int usePruning)
{
	// Safety checks.
	assert(node != NULL);
	assert(m_tree != NULL);

	int result = E_FAILURE;
//	char tmpbuf[1024];

	switch (node->type())
	{
	case AND:	// Expand an AND node (value).
		{
#ifdef DEBUG_TIME
			c_expand_AND_start = cpuTime();
#endif

			// Generate OR successors.
			CAONode* parent = node->parent();
			assert(parent != NULL);
			assert(parent->type() == OR);

			// Before expanding the node, must compute its g-value.
			int var = parent->label();
			int val = node->label();
			
			// Set current value.
			setValue(var, val);

			// purge cache for advanced caching
			if (CUTSET_CACHE)
			{
				if (CACHE_AT_AND_NODES)
				{
					variable_v::iterator it = m_advancedCacheFlagAND[var].begin();
					for ( ; it != m_advancedCacheFlagAND[var].end(); it++ )
					{
						m_cacheSep[(*it)]->purge();									
						if (COUNTING)
							m_cacheSepCounting[(*it)]->purge();
					}
				}
				
				if (CACHE_AT_OR_NODES)					
				{
					variable_v::iterator it = m_advancedCacheFlagOR[var].begin();
					for ( ; it != m_advancedCacheFlagOR[var].end(); it++ )
					{
						m_cacheFull[(*it)]->purge();									
						if (COUNTING) 
							m_cacheFullCounting[(*it)]->purge();
					}					
				}
			}

			int pos = m_position[var];
			CBucket* bucket = m_buckets->getBucketAt(pos);
			//assert(bucket != NULL);
			
			// Multiply all original functions in the bucket.
			// here g is 1 if the bucket is empty
			
			// robert: evaluate buckets for probabilities, initialize node with bucket result

			double el = bucket->multiply();	
			node->setL(el);			
			node->setG(1.0);
			
			// for counting solutions initialize with 1.0
			if(COUNTING)
				node->setGCount(1.0);


//			sprintf(tmpbuf, "\nLOG:  -[AND] initial g-value is set to %7g", g);
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			// Check for zero l-value. Prune under it.
///*
			if (0.0 == node->getL())
			{
				result = S_SUCCESS;
				zeroDeadEnds += 1;
				deadEnds += 1;
    			break;
			}	

			bool cons = false;	
			switch (usePruning)
			{
			case PRUNING_NO:
				{
					cons = true;
					break;
				}
			case PRUNING_CONSTRAINTS_ONLY:
				{
					// not just constraints
					//cons = forwardCheckingWithMarking(node);
					cons = forwardChecking(var, val);
					break;						
				}
			case PRUNING_FORWARD_CHECKING:
				{
					cons = forwardCheckingByProjection(var, val);
					//cons = (forwardCheckingWithMarking(node) && forwardChecking(node));
					break;
				}
			case PRUNING_FC_PROJECTION:
				{					
					//cons = (forwardCheckingWithMarking(node) && forwardChecking(node) && forwardCheckingByProjection(var, val));
					cons = (forwardChecking(var, val) && forwardCheckingByProjection(var,val));
					break;
				}			
			case PRUNING_FULL_LOOK_AHEAD:	
				break;
			}



			//bool cons = forwardChecking(var,val);
			//bool consProj = forwardCheckingByProjection(var,val);
			//if (cons && !consProj) 
			//	cout << "1111" << endl;
			//if (!cons && consProj)
			//{
			//	cout << "    2222" << endl;
			//	cons = forwardChecking(var,val);
			//	consProj = forwardCheckingByProjection(var,val);
			//}
			//if (!cons && !consProj)
			//	cout << "        3333" << endl;
			//if (cons && consProj)
			//	cout << "            4444" << endl;

			if (!cons)
			{
				node->setG(0.0);
				result = S_SUCCESS;
				deadEnds += 1;
				break;
			}	
			
			//*/

			// Get variable's children from the legal tree.
			variable_v children;
			m_tree->getChildren(var, children);

//			sprintf(tmpbuf, "\nLOG:  -[AND] children are: ");
			variable_v::iterator it = children.begin();
			for (; it != children.end(); ++it)
			{
				int ch = (*it);
				CAONode* child = new CAONode(OR, ch);
				assert(child != NULL);

				child->setParent(node);
				node->addChild(child);

				succ.push(child);

//					char tmp[10];
//					sprintf(tmp, "%d ", ch);
//					strcat(tmpbuf, tmp);
			}

			result = E_SUCCESS;

//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			// Free temporary buffers.
			children.clear();

#ifdef DEBUG_TIME
			c_expand_AND_end = cpuTime();
			tmExpandAND += (double)(c_expand_AND_end - c_expand_AND_start);
#endif

			break;
		}
	case OR:	// Expand an OR node (variable).
		{
#ifdef DEBUG_TIME
			c_expand_OR_start = cpuTime();
#endif

			// Generate AND successors.
			int var = node->label();
			vector<int> values;

			// Check for evidence variable.
			if (isEvidence(var))
			{
				int val = getValue(var);
				values.push_back(val);
			}
			else
			{
				for (int val = 0; val < getStaticDomainSize(var); ++val)
				{
					// add only consistent values
					if (m_validValues[var][val])
						values.push_back(val);
				}
			}

			// Create AND successors for each value.

			// Here we check different consistency levels:

			vector<int>::iterator it = values.begin();
			for (; it != values.end(); ++it)
			{
				int val = (*it);

				CAONode* child = new CAONode(AND, val);
				assert(child != NULL);

				// Set up links.
				child->setParent(node);
				node->addChild(child);

				succ.push(child);
			}

			result = E_SUCCESS;

			// Free temporary buffers.
			values.clear();


	
#ifdef DEBUG_TIME
			c_expand_OR_end = cpuTime();
			tmExpandOR += (double)(c_expand_OR_end - c_expand_OR_start);
#endif

			break;
		}
	}

	return result;
}

 //This function executes AND-OR tree search (full tree expansion).
int CProblemMixed::execAO(long tmLimit, double& tmCpuSearch, long& expansions, long& nodesAND, long& nodesOR, 
						  bool silent, int usePruning, int cacheSize)
{
	// Safety checks.
	assert(m_N > 0);
	assert(m_K > 0);
//robert	assert(m_C > 0);
	assert(m_tree != NULL);

	// DFS stacks.
	stack<CAONode*> succ;
	stack<CAONode*> OPEN;
	list<CAONode*> CLOSED;
	
	// General variables.
	int result;						// Output code.
	double c_start;				// Start timer.
	double c_end;					// End timer.
	char tmpbuf[1024];				// Buffer.

	
#ifdef DEBUG_TIME
	double c_cache_start;			// variables to clock the caching time
	double c_cache_end;
	double c_expand_start;
	double c_expand_end;
	double c_propagate_start;
	double c_propagate_end;
	double c_true_expand_start;
	double c_true_expand_end;
	double c_false_expand_start;
	double c_false_expand_end;
	double tmCache = 0.0;
	double tmExpand = 0.0;
	double tmPropagate = 0.0;
	double tmTrueExpand = 0.0;
	double tmFalseExpand = 0.0;
	//double tm4 = 0.0;
#endif
	
	// Intermediate results.
	long _timeNextSlice = 1;		// Next time point for intermediate time (every 1 sec).
	long _timeLastCheck = 0 ;		// Time when the time-check was done last time.
	// time_t time_counter ;
	unsigned long current_clock ;
	double tmSearch;

	// Set time zone from TZ environment variable.
	//_tzset();

	// Prologue.
	sprintf(tmpbuf, "\n BEGIN AO SEARCH  ----------------------------------------------");
	cout << tmpbuf;
	if (m_outputFile) fprintf(m_outputFile, tmpbuf);

	sprintf(tmpbuf, "\n   time limit: %d", tmLimit);
	cout << tmpbuf;
	if (m_outputFile) fprintf(m_outputFile, tmpbuf);
	
	sprintf(tmpbuf, "\n   silent mode: %s", ((!silent) ? "true" : "false"));
	cout << tmpbuf;
	if (m_outputFile) fprintf(m_outputFile, tmpbuf);

	if (!silent)
	{
		cout << "\n   DETAILS";
		if (m_outputFile) fprintf(m_outputFile, "\n   DETAILS");
	}

	// Start the timer (_clock)
	c_start = cpuTime();

	// Start the local clock
	// startClock(&time_counter) ;

	// Initialize search.
	expansions = 0;
	nodesAND = 0;
	nodesOR = 0;
	deadEnds = 0;
	zeroDeadEnds = 0;
	m_probabilityOfQuery = -1;

	// Set first variable.
	int var0 = m_ordering[0];
	CAONode* root = new CAONode(OR, var0);
	assert(root != NULL);

	// Initialize the variable stack.
	OPEN.push(root);

	result = S_FAILURE;

	percentageDone = 0.0;
	// Start search.
	while (!OPEN.empty())
	{
		// EXPAND step : 
		//   1a) pick top of OPEN.
		//   1b) generate successors: 
		//	      1b-1) AND nodes - values of the variable.
		//        1b-2) OR nodes - child nodes of the variable in legal tree.

		// PROPAGATION step :
		//    1) propagate g-values to ancestors.
		//    2) EXPAND.

//expand:
#ifdef DEBUG_TIME
		c_expand_start = cpuTime();
#endif

		// Node expansion.
		CAONode* node = OPEN.top();
		OPEN.pop();
		CLOSED.push_front(node);

		// Check for silent run.
		if (!silent)
		{
			current_clock = cpuTime();
			if (current_clock >= _timeNextSlice && 
				_timeLastCheck < _timeNextSlice) 
			{
				// print intemediate results
				double current_time = (double) current_clock - c_start;
				
				sprintf(tmpbuf, "\n    -int: time %7g, nodes %d", current_time, expansions);

				cout << tmpbuf;
				
				// percentage done
				cout << "       done: " << percentageDone << " %";

				if (m_outputFile) fprintf(m_outputFile, tmpbuf);
				
				// Compute next time check bound as next 1024 multiple time point from the current one.
				// int _timeFraction = current_clock & 7 ;
				//_timeNextSlice = current_clock + 1024 - _timeFraction ;
				_timeNextSlice = current_clock + 3.0;

				// Save current time for next round.
				_timeLastCheck = current_clock ;
			}
		}

		// Check for time limit violation
		if ((tmLimit > 0) && (cpuTime() - c_start >= tmLimit))
		{
			result = S_TIMEOUT;		// timeout.

			break;
		}

#ifdef DEBUG_TIME
		c_false_expand_start = cpuTime();
#endif

	
#ifdef DEBUG_TIME
		c_false_expand_end = cpuTime();
		tmFalseExpand += (double)(c_false_expand_end - c_false_expand_start);
#endif

		bool inCache = false;
		// check cache

		if (CACHE_AT_OR_NODES)
		{
#ifdef DEBUG_TIME
			c_cache_start = cpuTime();
#endif	
			if ( (cacheSize > 0) && (node->type() == OR) && (node != root) )
			{
				int var = node->label();
				if(m_cacheFull[var]->getArgc() <= cacheSize)
				{
					double temp;
					double tempCounting;
					m_cacheFull[var]->getCurrentValue(&temp);
					if (temp >= 0)
					{
						inCache = true;
						node->setG(temp);
					}
										
					if (COUNTING)
					{
						m_cacheFullCounting[var]->getCurrentValue(&tempCounting);
					
						// Robert: attention, this may cause problems if inCache is set in 2 places!
						if (tempCounting >= 0)
						{
							//inCache = true;
							node->setGCount(tempCounting);
						}
					}
				}	
			}
		
#ifdef DEBUG_TIME
			c_cache_end = cpuTime();
			tmCache += (double)(c_cache_end - c_cache_start);
#endif
		} // end CACHE_AT_OR_NODES


 
		
		if (CACHE_AT_AND_NODES)
		{
#ifdef DEBUG_TIME
			c_cache_start = cpuTime();
#endif
			
			if ( (cacheSize > 0) && (node->type() == AND) )
			{	
				// Set current value temporarily
				int var = node->parent()->label();
				int val = node->label();			
				setValue(var, val);

				if(m_cacheSep[var]->getArgc() <= cacheSize)
				{
					double temp;
					double tempCounting;
					m_cacheSep[var]->getCurrentValue(&temp);
					if (temp >= 0)
					{
						inCache = true;
						node->setG(temp);
						// Robert - need to write the L value if we take G-value from cache
						int pos = m_position[var];
						CBucket* bucket = m_buckets->getBucketAt(pos);
						double el = bucket->multiply();	
						node->setL(el);	
					}
					
					if (COUNTING)
					{
						m_cacheSepCounting[var]->getCurrentValue(&tempCounting);					
						
						if (tempCounting >= 0)
						{
							inCache = true;
							node->setGCount(tempCounting);
							// Robert - need to write the L value if we take G-value from cache
						}
					}
				}				
			}
#ifdef DEBUG_TIME
			c_cache_end = cpuTime();
			tmCache += (double)(c_cache_end - c_cache_start);
#endif
		} // end CACHE_AT_AND_NODES

		
#ifdef DEBUG_TIME
		c_true_expand_start = cpuTime();
#endif
		if(!inCache)
		{
			// start Bucket Elimination if possible, if node is not in cutset
			if ( (SWITCH_TO_BE) && (node->type() == OR) && (m_wCutset[node->label()] <= cacheSize) )
			{
				double gValue;
				node->computeExactSubtree(this, m_N, gValue);
				double temp = 1.0;
				for (int i=0; i < m_descendants[node->label()].size() - 1 ; i++)
					temp *= scalingFactor;
				node->setG(gValue * temp);
			}
			else
			{
				switch (expand(node, succ, false, UNKNOWN, usePruning))
				{					
				case E_SUCCESS:
					{
						// Add successors to stack.
						while (!succ.empty())
						{
							CAONode* child = succ.top();
							succ.pop();

							OPEN.push(child);
						}
						
						// Increase number of expansions.
						++expansions;

						if (AND == node->type())
							++nodesAND;
						else
							++nodesOR;

						break;
					}
				case E_FAILURE:
					{
						result = S_FAILURE;
						goto done;
					}
				}
			}
		}
#ifdef DEBUG_TIME
		c_true_expand_end = cpuTime();
		tmTrueExpand += (double)(c_true_expand_end - c_true_expand_start);

		c_expand_end = cpuTime();
		tmExpand += (double)(c_expand_end - c_expand_start);
#endif

// Propagation.

#ifdef DEBUG_TIME
		c_propagate_start = cpuTime();
#endif

		CAONode* cand = node;
		//bool firstPassPropagation = true;
		while (0 == cand->childrenSize())
		{
			CAONode* parent = cand->parent();
			bool complete = false;

			switch (cand->type())
			{
			case AND:
				{	
					// undo list of modifications;
					cand->cleanInvalidValues(m_validValues,m_validValuesCount, m_assignment);

					// update percentage done
//					if(parent->label() == percentageDoneVariable)
//						percentageDone = updatePercentageDone();

					// Update parent's g-value.
					//double g = max(cand->getG(), parent->getG());
					double g = cand->getG() * cand->getL();
					if (cand->getG() != 0  &&  cand->getL() != 0  &&  g == 0 )
						g = exp(-500.0);
					//double g = (cand->getG() + parent->getG());
					parent->setG( g + parent->getG() );

					// for counting solutions
					double gCounting = (cand->getGCount() + parent->getGCount());					
					parent->setGCount(gCounting);

			

					if (CACHE_AT_AND_NODES)
					{
						// check if need to save in cache:
						
#ifdef DEBUG_TIME
						c_cache_start = cpuTime();
#endif
						if ( cacheSize > 0 )
						{
							int var = cand->parent()->label();
							if(m_cacheSep[var]->getArgc() <= cacheSize)
							{
								double temp = cand->getG();
								m_cacheSep[var]->setCurrentValue(&temp);
								if (COUNTING)
								{
									double tempCounting = cand->getGCount();
									m_cacheSepCounting[var]->setCurrentValue(&tempCounting);
								}
							}
						}
						
#ifdef DEBUG_TIME
						c_cache_end = cpuTime();
						tmCache += (double)(c_cache_end - c_cache_start);
#endif
					} // end CACHE_AT_AND_NODES

					break;
				}
			
			case OR:
				{
					cand->setG( cand->getG() * scalingFactor );

					// Update parent's g-value.
					if (NULL == parent)
					{
						m_probabilityOfQuery = cand->getG() ;

						removeNode(CLOSED, cand);
						delete cand;

						complete = true;	// Complete search.
						break;
					}
					else
					{

						if (CACHE_AT_OR_NODES)
						{
							// check if need to save in cache:
							
#ifdef DEBUG_TIME
							c_cache_start = cpuTime();
#endif
							if ( cacheSize > 0 )
							{
								int var = cand->label();
								if(m_cacheFull[var]->getArgc() <= cacheSize)
								{
									double temp = cand->getG() / scalingFactor;
									m_cacheFull[var]->setCurrentValue(&temp);
									
									if (COUNTING)
									{
										double tempCounting = cand->getGCount();
										m_cacheFullCounting[var]->setCurrentValue(&tempCounting);
									}
								}				
							}
							
#ifdef DEBUG_TIME
							c_cache_end = cpuTime();
							tmCache += (double)(c_cache_end - c_cache_start);
#endif
						} // end CACHE_AT_OR_NODES

						//propagate value to parent						
						double g = cand->getG() * parent->getG();
						if ( cand->getG() != 0  &&  parent->getG() != 0  &&  g == 0 )
							parent->setG ( exp(-500.0) );
						else
							parent->setG(g);

						//propagate countvalue to parent
						double gCounting = cand->getGCount() * parent->getGCount();
						parent->setGCount(gCounting);

						// prune the other children of parent, and remove from OPEN
						if (g == 0) 
						{
							for (int i=0; i < parent->childrenSize() - 1 ; i++)
							{
								CAONode* eraseNode = OPEN.top();
								OPEN.pop();
								delete eraseNode;
							}
							parent->clear();
						}
					}

					break;
				}			
			}

			//firstPassPropagation = false;
			
			if (complete)
				break;

			// Remove candidate from children set of the parent.
			if(cand)
				parent->remove(cand);

			// Remove candidate from CLOSED set.
			if(cand->type() == AND)
				m_assignment[parent->label()] = -1;
			removeNode(CLOSED, cand);
			

			m_probabilityOfQuery = root->getG();
			m_numberOfSolutions = root->getGCount();
			delete cand;

			// Get next candidate for propagation.
			cand = CLOSED.front();
		}
#ifdef DEBUG_TIME
		c_propagate_end = cpuTime();
		tmPropagate += (double)(c_propagate_end - c_propagate_start);
#endif
	}

	result = S_SUCCESS;

done:

	// Done.
	c_end = cpuTime();

//	current_clock = stopClock(0, &time_counter) ;
//	double tm1 = (double)current_clock / 1000.0 ;

	double tm2 = (double)(c_end - c_start);
	tmSearch = tm2;

	tmCpuSearch = tmSearch;

	// Output solution.
	switch (result)
	{
	case S_SUCCESS:
		{
			sprintf(tmpbuf, "\n   out of time: false");
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			globalCostFromEvidence = globalCostFromEvidence * m_probabilityOfQuery;

			sprintf(tmpbuf, "\n   component log probability of evidence: %7g", (log(m_probabilityOfQuery)));
			sprintf(tmpbuf, "\n   component probability of evidence: %7g", (m_probabilityOfQuery));
			cout << tmpbuf;			
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			//sprintf(tmpbuf, "\n   number of solutions: %g", m_numberOfSolutions);
			//cout << tmpbuf;			
			//if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			//sprintf(tmpbuf, "\n   optimal solution cost: %7g", m_bestSolutionCost);
			//cout << tmpbuf;
			//if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			sprintf(tmpbuf, "\n   CPU time (search):          %-7g", tmSearch);
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

////////////////////////////////////////////////
#ifdef DEBUG_TIME
			sprintf(tmpbuf, "\n   CACHE time :                %-7g", tmCache);
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);


			sprintf(tmpbuf, "\n   EXPAND time :               %-7g", tmExpand);
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);	

			sprintf(tmpbuf, "\n   TRUE EXPAND time :          %-7g", tmTrueExpand);
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			sprintf(tmpbuf, "\n   FALSE EXPAND time :         %-7g", tmFalseExpand);
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			sprintf(tmpbuf, "\n   EXPAND OR time :            %-7g", tmExpandOR);
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			sprintf(tmpbuf, "\n   EXPAND AND time :           %-7g", tmExpandAND);
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			sprintf(tmpbuf, "\n   PROPAGATE time :            %-7g", tmPropagate);
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
#endif
////////////////////////////////////////////////

			sprintf(tmpbuf, "\n   number of expansions:        %d", expansions);
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			sprintf(tmpbuf, "\n   number of AND expansions:    %d", nodesAND);
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			sprintf(tmpbuf, "\n   number of OR expansions:     %d", nodesOR);
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			cout << "\n END AO SEARCH    ----------------------------------------------\n";	
			if (m_outputFile) fprintf(m_outputFile, "\n END AO SEARCH    ----------------------------------------------\n");

			break;
		}
	case S_TIMEOUT:
		{
			sprintf(tmpbuf, "\n   out of time: true");
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			sprintf(tmpbuf, "\n   best solution cost: %7g", m_bestSolutionCost);
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			sprintf(tmpbuf, "\n\n   CPU time (search): %7g", tmSearch);
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			sprintf(tmpbuf, "\n   number of expansions: %d", expansions);
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			sprintf(tmpbuf, "\n   number of AND expansions: %d", nodesAND);
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			sprintf(tmpbuf, "\n   number of OR expansions: %d", nodesOR);
			cout << tmpbuf;
			if (m_outputFile) fprintf(m_outputFile, tmpbuf);

			cout << "\n END AO SEARCH    ----------------------------------------------\n";	
			if (m_outputFile) fprintf(m_outputFile, "\nEND.\n");

			break;
		}
	};

	// Clear memory.

	return result;
}

//// This function executes AND-OR tree search with MBE pruning (Alpha-Beta pruning style).
//int CProblemMixed::execAOMB(int ibound, long tmLimit, double& tmCpuSearch, 
//	long& expansions, long& nodesAND, long& nodesOR, bool silent)
//{
//	// Safety checks.
//	assert(m_N > 0);
//	assert(m_K > 0);
////robert	assert(m_C > 0);
//	assert(m_tree != NULL);
//
//	// Log file.
////	m_outputFile = fopen("log_aomb.txt", "w");
////	assert(m_outputFile != NULL);
//
//	// DFS stacks.
//	stack<CAONode*> succ;
//	list<CAONode*> OPEN;
//	list<CAONode*> CLOSED;
//	
//	// General variables.
//	double ub;
//	int result;						// Output code.
//	long pruned;
//	clock_t c_start;				// Start timer.
//	clock_t c_end;					// End timer.
//	char tmpbuf[1024];				// Buffer.
//
//	// Intermediate results.
//	long _timeNextSlice = 1024;		// Next time point for intermediate time (every 1 sec).
//	long _timeLastCheck = 0 ;		// Time when the time-check was done last time.
//	time_t time_counter ;
//	unsigned long current_clock ;
//	double tmSearch;
//
//	// Set time zone from TZ environment variable.
//	//_tzset();
//
//	// Prologue.
//	sprintf(tmpbuf, "\nAOMB started ...");
//	cout << tmpbuf;
//	if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//	sprintf(tmpbuf, "\n   time limit: %d", tmLimit);
//	cout << tmpbuf;
//	if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//	sprintf(tmpbuf, "\n   i-bound: %d", ibound);
//	cout << tmpbuf;
//	if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//	
//	sprintf(tmpbuf, "\n   silent mode: %s", (silent ? "true" : "false"));
//	cout << tmpbuf;
//	if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//	if (!silent)
//	{
//		cout << "\n   DETAILS";
//		if (m_outputFile) fprintf(m_outputFile, "\n   DETAILS");
//	}
//
//	// Start the timer (_clock)
//	c_start = clock();
//
//	// Start the local clock
//	startClock(&time_counter) ;
//
//	// Initialize search.
//	expansions = 0;
//	pruned = 0;
//	nodesAND = 0;
//	nodesOR = 0;
//
//	// Set first variable.
//	int var0; var0 = m_ordering[0];
//	CAONode* root = new CAONode(OR, var0);
//	assert(root != NULL);
//
//	// Initialize the variable stack.
//	root->buildValueCosts(this, ibound, ub);
//	OPEN.push_front(root);
//
////	sprintf(tmpbuf, "\nLOG: root node is %d", root->label());
////	if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//	result = S_FAILURE;
//
//	// Start search.
//	while (!OPEN.empty())
//	{
//		// EXPAND step : 
//		//   1a) pick top of OPEN.
//		//   1b) generate successors: 
//		//	      1b-1) AND nodes - values of the variable.
//		//        1b-2) OR nodes - child nodes of the variable in legal tree.
//
//		// PROPAGATION step :
//		//    1) propagate g-values to ancestors.
//		//    2) EXPAND.
//
////expand:
//
//		// Check for silent run.
//		if (!silent)
//		{
//			current_clock = stopClock(0, &time_counter);
//			if (current_clock >= _timeNextSlice && 
//				_timeLastCheck < _timeNextSlice) 
//			{
//				// print intemediate results
//				double current_time = ((double) current_clock) / 1000.0 ;
//				
//				sprintf(tmpbuf, "\n    -int: time %7g, nodes %d, pruned %d, opened %d, closed %d", 
//					current_time, expansions, pruned, OPEN.size(), CLOSED.size());
//
//				cout << tmpbuf;
//				if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//				
//				// Compute next time check bound as next 1024 multiple time point from the current one.
//				int _timeFraction = current_clock & 7 ;
//				_timeNextSlice = current_clock + 1024 - _timeFraction ;
//
//				// Save current time for next round.
//				_timeLastCheck = current_clock ;
//			}
//		}
//
//		// Check for time limit violation
//		if ((tmLimit > 0) && (clock() - c_start >= tmLimit))
//		{
//			result = S_TIMEOUT;		// timeout.
//
//			goto done;
//		}
//
//		// Node expansion.
//		CAONode* node = OPEN.front();
//		OPEN.pop_front();
//		CLOSED.push_front(node);
//
////		sprintf(tmpbuf, "\nLOG: expanding %s node %d ...", 
////			((node->type() == AND) ? "AND" : "OR"), node->label());
////		if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//		switch (expand(node, succ, true, ibound))
//		{
//		case E_SUCCESS:
//			{
//				// Add successors to stack.
//				while (!succ.empty())
//				{
//					CAONode* child = succ.top();
//					succ.pop();
//
//					OPEN.push_front(child);
//				}
//				
//				// Increase number of expansions.
//				++expansions;
//				
//				if (AND == node->type())
//					++nodesAND;
//				else
//					++nodesOR;
//
//				break;
//			}
//		case E_FAILURE:
//			{
//				result = S_FAILURE;
//				goto done;
//			}
//		}
//
//
//		// Propagation.
//		CAONode* cand = node;
//		while (0 == cand->childrenSize())
//		{
//			CAONode* parent = cand->parent();
//			bool complete = false;
//
////			sprintf(tmpbuf, "\nLOG: propagation at %s node %d ...", 
////				((cand->type() == AND) ? "AND" : "OR"), cand->label());
////			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			switch (cand->type())
//			{
//			case AND:
//				{
//					// Update parent's g-value.
//					double g = max(cand->getG(), parent->getG());
//					parent->setG(g);
//					parent->setUpdated(true);
//
////					sprintf(tmpbuf, "\nLOG:  -[AND] update OR parent g-value to %7g", g);
////					if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//					
//					// Prune all siblings that have an estimated cost
//					// not better than the g-value that has been just propagated.
//					vector<CAONode*> siblings;
//					parent->prune(cand, g, siblings);
//					
//					vector<CAONode*>::iterator it = siblings.begin();
//					for (; it != siblings.end(); ++it)
//					{
//						CAONode* tmp = (*it);
//						removeNode(OPEN, tmp);
//						removeNode(CLOSED, tmp);
//						
//						++pruned;
//	
////						sprintf(tmpbuf, "\nLOG:  -[AND] prune sibling %d having estimate %7g", tmp->label(), tmp->cost());
////						if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//						delete tmp;
//					}
//
//					// Free temporary buffers.
//					siblings.clear();
//
//					break;
//				}
//			case OR:
//				{
//					// Update parent's g-value.
//					if (NULL == parent)
//					{
//						m_bestSolutionCost = cand->getG();
//
//						removeNode(CLOSED, cand);
//						delete cand;
//
//						complete = true;	// Complete search.
//						break;
//					}
//					else
//					{
//						double g = cand->getG() * parent->getG();
//						parent->setG(g);
//						parent->setUpdated(true);
//
////						sprintf(tmpbuf, "\nLOG:  -[OR] update AND parent g-value to %7g", g);
////						if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//					}
//
//					break;
//				}
//			}
//			
//			if (complete)
//				break;
//
//			// Remove candidate from children set of the parent.
//			parent->remove(cand);
//
//			// Remove candidate from CLOSED set.
//			removeNode(CLOSED, cand);
//
//			delete cand;
//
//			// Get next candidate for propagation.
//			cand = CLOSED.front();
//		}
//	}
//
//	result = S_SUCCESS;
//
//done:
//
//	// Done.
//	c_end = clock();
//	current_clock = stopClock(0, &time_counter) ;
//	double tm1 = (double)current_clock / 1000.0 ;
//	double tm2 = (double)(c_end - c_start) / 1000.0;
//	tmSearch = tm2;
//
//	tmCpuSearch = tmSearch;
//
//	// Output solution.
//	switch (result)
//	{
//	case S_SUCCESS:
//		{
//			sprintf(tmpbuf, "\n   out of time: false");
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n   optimal solution cost: %7g", m_bestSolutionCost);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n\n   CPU time (search): %7g", tmSearch);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n   number of expansions: %d", expansions);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n   number of pruned: %d", pruned);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n   number of AND expansions: %d", nodesAND);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n   number of OR expansions: %d", nodesOR);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			cout << "\nEND.\n";
//			if (m_outputFile) fprintf(m_outputFile, "\nEND.\n");
//
//			break;
//		}
//	case S_TIMEOUT:
//		{
//			sprintf(tmpbuf, "\n   out of time: true");
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n   best solution cost: %7g", m_bestSolutionCost);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n\n   CPU time (search): %7g", tmSearch);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n   number of expansions: %d", expansions);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n   number of AND expansions: %d", nodesAND);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			sprintf(tmpbuf, "\n   number of OR expansions: %d", nodesOR);
//			cout << tmpbuf;
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			cout << "\nEND.\n";
//			if (m_outputFile) fprintf(m_outputFile, "\nEND.\n");
//
//			break;
//		}
//	};
//
//	// Clear memory.
//	list<CAONode*>::iterator ito = OPEN.begin();
//	for (; ito != OPEN.end(); ++ito)
//		delete (*ito);
//	OPEN.clear();
//
//	list<CAONode*>::iterator itc = CLOSED.begin();
//	for (; itc != CLOSED.end(); ++itc)
//		delete (*itc);
//	CLOSED.clear();
//
////	fclose(m_outputFile);
////	m_outputFile = NULL;
//
//	return result;
//}
//
//Checks if the scope of the function is included in the current assignment
bool CProblemMixed::isEligibleToBeChecked(CFunction* fun)
{
	int* argv = fun->getArgv();
	int argc = fun->getArgc();

	for (int i=0; i < argc; i++)
	{		
		if (m_assignment[argv[i]] < 0)
			return false;
	}

	return true;
}




// Checks the consistency of the current assignment ending in (var,val) for and/or search in mixed networks
// m_assignment is consistent for the current path, but not necessarily consistent for other variables not in path
bool CProblemMixed::consistent_old(int var, int val)
{
	int old_value = m_assignment[var];
	m_assignment[var] = val;

	function_map& funs = m_functionsConstraint;
	function_map::iterator it = funs.begin();
	for (; it != funs.end(); ++it)
	{
		CFunction* fun = (*it).second;
		if (!isEligibleToBeChecked(fun))
			continue;

		int v;
		fun->getCurrentValue(&v);

		if (v == 0)
		{
			// set old value back
			m_assignment[var] = old_value;
			return false;
		}
	}
	
	// set old value back
	m_assignment[var] = old_value;
	return true;
}



// Checks the consistency of the current assignment ending in (var,val) for and/or search in mixed networks
// m_assignment is consistent for the current path, but not necessarily consistent for other variables not in path
bool CProblemMixed::consistent(int var, int val)
{
	int old_value = m_assignment[var];
	m_assignment[var] = val;

	function_v::iterator it = m_variablesInFunctions[var].begin();
	for (; it != m_variablesInFunctions[var].end(); ++it)
	{
		CFunction* fun = (*it);
		if (!isEligibleToBeChecked(fun))
			continue;

		int v;
		fun->getCurrentValue(&v);

		if (v == 0)
		{
			// set old value back
			m_assignment[var] = old_value;
			return false;
		}
	}
	
	// set old value back
	m_assignment[var] = old_value;
	return true;
}

// Checks the consistency of the current assignment ending in var by FC for and/or search in mixed networks
bool CProblemMixed::forwardChecking_old(int var, int val)
{
	int old_value = m_assignment[var];
	m_assignment[var] = val;

	// Check all the variables which are descendants of var
	variable_v::iterator itd = m_descendants[var].begin();
	for (; itd != m_descendants[var].end(); ++itd)
	{
		int v = (*itd);

		if (v == var) continue;

		//if (isEvidence(v))
		//{
		//	if (consistent (v, m_assignment[v]))
		//		continue;
		//	else 
		//	{
		//		// set old value back
		//		m_assignment[var] = old_value;
		//		return false;
		//	}
		//}

		bool cons = false;
		for (int i=0; i < m_domains[v]; i++)
			if (consistent(v, i))
			{
				cons = true;
				break;
			}

		if (cons == false)
		{
			// set old value back
			m_assignment[var] = old_value;
			return false;
		}
	}

	// set old value back
	m_assignment[var] = old_value;
	return true;
}



// Checks the consistency of the current assignment ending in var by FC for and/or search in mixed networks
bool CProblemMixed::forwardChecking (int var, int val)
{
	int old_value = m_assignment[var];
	m_assignment[var] = val;

	// Check all the variables which are descendants of var
	variable_v::iterator itd = m_descendants[var].begin();
	for (; itd != m_descendants[var].end(); ++itd)
	{
		int v = (*itd);

		if (v == var) continue;

		if (isEvidence(v))
		{
			if (consistent (v, m_assignment[v]))
				continue;
			else 
			{
				// set old value back
				m_assignment[var] = old_value;
				return false;
			}
		}

		bool cons = false;
		for (int i=0; i < m_domains[v]; i++)
			if (consistent(v, i))
			{
				cons = true;
				break;
			}

		if (cons == false)
		{
			// set old value back
			m_assignment[var] = old_value;
			return false;
		}
	}

	// set old value back
	m_assignment[var] = old_value;
	return true;
}

bool CProblemMixed::forwardChecking(CAONode* node)
{	
	assert(node->type() == AND);

	int var = node->parent()->label();
	int val = node->label();

	queue<int> vQueue;
	vQueue.push(var);

	while(vQueue.empty())
	{	
		int vLiteral = vQueue.front();
		vQueue.pop();

		// iterate through the active functions
		function_v::iterator it = m_variablesInFunctions[vLiteral].begin();
		for (; it != m_variablesInFunctions[vLiteral].end(); ++it)
		{
			CFunction* fun = (*it);

			int argc = fun->getArgc();
			int* argv = fun->getArgv();
			
			int count = 0;
			int v;
			for (int i=0 ; i < argc ; i++)
			{
				int tempv = argv[i];
				if (m_assignment[tempv] < 0 )
				{
					v = tempv;
					count++;
				}
				if (count > 1) 
					break; // too many unassigned vars
			}

			if (count != 1)
				continue;

			// now v is the only one unassigned variable in fun
			int potentialSingleValue = -1;
			for (int i=0; i < m_domains[v]; i++)
			{
				if (!m_validValues[v][i]) 
					continue;

				setValue(v,i);
				int value;
				((CConstraintTable*)fun)->getCurrentValue(&value);
				setValue(v,-1);
				if (value == 0)
				{
					m_validValues[v][i] = false;
					--m_validValuesCount[v];
					node->addInvalidValue(v,i);
					if (m_validValuesCount[v] == 0)
						return false;
				}
				else
					potentialSingleValue = i;
			}	
			if (m_validValuesCount[v] == 1)
			{
				assert(potentialSingleValue != -1);
				setValue(v, potentialSingleValue);
				vQueue.push(v);
			}
		}
	}
	return true;

///////////////////////////////////////////////////////////////////////
	//// Check all the variables which are descendants of var
	//variable_v::iterator itd = m_descendants[var].begin();
	//for (; itd != m_descendants[var].end(); ++itd)
	//{
	//	int v = (*itd);

	//	if (v == var) continue;

	//	assert(!isEvidence(v));
	//	//if (isEvidence(v))
	//	//{
	//	//	if (consistent (v, m_assignment[v]))
	//	//		continue;
	//	//	else 
	//	//	{
	//	//		return false;
	//	//	}
	//	//}

	//	for (int i=0; i < m_domains[v]; i++)
	//	{
	//		if (!m_validValues[v][i]) 
	//			continue;

	//		if (!consistent(v, i))
	//		{
	//			m_validValues[v][i] = false;
	//			--m_validValuesCount[v];
	//			node->addInvalidValue(v,i);
	//			if (m_validValuesCount[v] == 0)
	//				return false;
	//		}
	//		else
	//			break;
	//	}		
	//}
	//
	//return true;

}

bool CProblemMixed::forwardCheckingWithMarking(CAONode* node)
{
	assert (node->type() == AND);

	int var = node->parent()->label();
	int val = node->label();

	variable_s activeVars; // collect active variables (that appear unassigned in scopes together with var)
	
	// iterate through the active functions
	function_v::iterator it = m_variablesInFunctions[var].begin();
	for (; it != m_variablesInFunctions[var].end(); ++it)
	{
		CFunction* fun = (*it);

		int argc = fun->getArgc();
		int* argv = fun->getArgv();
		
		for (int i=0 ; i < argc ; i++)
		{
			int v = argv[i];
			if (m_assignment[v] < 0 &&
				activeVars.find(v) == activeVars.end() )
				activeVars.insert(v);
		}
	}

	// now check the active variables
	variable_s::iterator itvar = activeVars.begin();
	for ( ; itvar != activeVars.end(); ++itvar )
	{
		int v = (*itvar);

		function_v::iterator itfunc = m_variablesInFunctions[v].begin();
		for (; itfunc != m_variablesInFunctions[v].end(); ++itfunc)
		{
			CFunction* f = (*itfunc);
			if ( ! supportIn(f, v, node) )
			{
				activeVars.clear();
				return false;
			}
		}
	}

	activeVars.clear();
	return true;
}

// Checks the consistency by looking at all constraints, projecting over the current assignment
bool CProblemMixed::forwardCheckingByProjection_old(int var, int val)
{
	int old_value = m_assignment[var];
	m_assignment[var] = val;

	// Check for every constraint
	function_map& funs = m_functionsConstraint;
	function_map::iterator it = funs.begin();
	for (; it != funs.end(); ++it)
	{
		CFunction* fun = (*it).second;
	
		// Check if the constraint is relevant (if scope intersects with current assignment)
		int argc = fun->getArgc();
		int* argv = fun->getArgv();
		// Safety checks.
		assert(argc > 0);
		assert(argv != NULL);

		int* isAssigned = new int[argc]; // array of 0 and 1 to indicate if the corresponding variable is in the current path
		bool isRelevantConstraint = false; // flag to indicate if the scope intersects with the current path
		bool isFullyAssigned = true; // flag to indicate if the constraint has its scope fully assigned;

		int i;
		for (i=0; i < argc; i++)
		{		
			if (m_assignment[argv[i]] < 0)
			{
				isAssigned[i] = 0;
				isFullyAssigned = false;
			}
			else
			{
				isAssigned[i] = 1;
				isRelevantConstraint = true;
			}	
		}

		if (!isRelevantConstraint || isFullyAssigned) // go to the next constraint
		{
			if (isAssigned)
				delete[] isAssigned;
			break;
		}

		// now we have a relevant constraint, and the array isAssigned which gives the assigned vars in the scope
		// start checking the projection on the assigned vars for consistency
		// enumerate all possible combinations of the non-assigned vars in the scope, for their corresponding
		// positions. stop when we find one that is consistent (because this will be in the projection)

		// find the smallest variable in the scope (rightmost one) not assigned
		int smallestVarNotYetAssigned = -1;
		for (i = argc-1; i >= 0 ; i--)
			if (isAssigned[i] == 1)
				continue;
			else
			{
				smallestVarNotYetAssigned = i;
				break;
			}

		int* parents = new int[argc]; // this is the scope of the constraint
		for (i = 0; i < argc ; i++) 
		{
			if(isAssigned[i] == 1)
				parents[i] = m_assignment[argv[i]]; // assigned values are maintained
			else 
				parents[i] = 0;	// enumerate only through those that are not assigned;
		}
		parents[smallestVarNotYetAssigned] = -1; // to initiate the enumeration properly

		// compute how many assignments we have to enumerate
		int nCount = 1;
		for (i = 0; i < argc; i++)
		{
			if (isAssigned[i] == 0)
				nCount *= getStaticDomainSize(argv[i]);
		}

		bool isConsistent = false;
		// Begin enumeration
		for (i =0; i < nCount; i++)
		{

			// Find next combination from the current.
			for (int pcai = argc - 1 ; pcai >= 0 ; --pcai) 
			{
				if (isAssigned[pcai] == 1)
					continue;

				if (++parents[pcai] < getStaticDomainSize(pcai)) break ;
				parents[pcai] = 0 ;
			}

			// now check the value of the constraint for the current assignment
			int k = argc - 1;
			int addr = parents[k];
			assert(addr >= 0);
			
			int j = getStaticDomainSize(argv[k]);
			for (--k ; k >= 0 ; k--)
			{
				// here we should add to 'adr' the product of :
				//		#-of-values-of-variable-k * current-value-of-k
				int value = parents[k];
				assert(value >= 0);

				addr += j*value ;
				j *= getStaticDomainSize(argv[k]) ;
			}

			if( ((CConstraintTable*)fun)->getValueAt(addr) == 1)
			{
				isConsistent = true;
				break;
			}
		}

		if(parents)
			delete[] parents;

		if (isAssigned)
				delete[] isAssigned;

		if (!isConsistent)
		{			
			// set old value back
			m_assignment[var] = old_value;
			return false;
		}
	}

	// set old value back
	m_assignment[var] = old_value;
	return true;
}





// Checks the consistency by looking at all constraints, projecting over the current assignment
bool CProblemMixed::forwardCheckingByProjection(int var, int val)
{
	int old_value = m_assignment[var];
	m_assignment[var] = val;

	// Check for every constraint
	function_v funs = m_variablesInFunctions[var];
	function_v::iterator it = funs.begin();
	for (; it != funs.end(); ++it)
	{
		CFunction* fun = (*it);
	
		// Check if the constraint is relevant (if scope intersects with current assignment)
		int argc = fun->getArgc();
		int* argv = fun->getArgv();
		// Safety checks.
		assert(argc > 0);
		assert(argv != NULL);

		int* isAssigned = new int[argc]; // array of 0 and 1 to indicate if the corresponding variable is in the current path
		bool isRelevantConstraint = false; // flag to indicate if the scope intersects with the current path
		bool isFullyAssigned = true; // flag to indicate if the constraint has its scope fully assigned;

		int i;
		//int varAssigned = 0;
		for (i=0; i < argc; i++)
		{		
			if (m_assignment[argv[i]] < 0)
			{
				isAssigned[i] = 0;
				isFullyAssigned = false;
			}
			else
			{
				isAssigned[i] = 1;
				//varAssigned++;
				isRelevantConstraint = true;
			}	
		}

		// if too many unassigned, break   (argc-varAssigned > 1) ||
		if ( !isRelevantConstraint || isFullyAssigned) // go to the next constraint
		{
			if (isAssigned)
				delete[] isAssigned;
			break;
		}

		// now we have a relevant constraint, and the array isAssigned which gives the assigned vars in the scope
		// start checking the projection on the assigned vars for consistency
		// enumerate all possible combinations of the non-assigned vars in the scope, for their corresponding
		// positions. stop when we find one that is consistent (because this will be in the projection)

		// find the smallest variable in the scope (rightmost one) not assigned
		int smallestVarNotYetAssigned = -1;
		for (i = argc-1; i >= 0 ; i--)
			if (isAssigned[i] == 1)
				continue;
			else
			{
				smallestVarNotYetAssigned = i;
				break;
			}

		int* parents = new int[argc]; // this is the scope of the constraint
		for (i = 0; i < argc ; i++) 
		{
			if(isAssigned[i] == 1)
				parents[i] = m_assignment[argv[i]]; // assigned values are maintained
			else 
				parents[i] = 0;	// enumerate only through those that are not assigned;
		}
		parents[smallestVarNotYetAssigned] = -1; // to initiate the enumeration properly

		// compute how many assignments we have to enumerate
		int nCount = 1;
		for (i = 0; i < argc; i++)
		{
			if (isAssigned[i] == 0)
				nCount *= getStaticDomainSize(argv[i]);
		}

		bool isConsistent = false;
		// Begin enumeration
		for (i =0; i < nCount; i++)
		{

			// Find next combination from the current.
			for (int pcai = argc - 1 ; pcai >= 0 ; --pcai) 
			{
				if (isAssigned[pcai] == 1)
					continue;

				if (++parents[pcai] < getStaticDomainSize(argv[pcai])) break ;
				parents[pcai] = 0 ;
			}

			// now check the value of the constraint for the current assignment
			int k = argc - 1;
			int addr = parents[k];
			assert(addr >= 0);
			
			int j = getStaticDomainSize(argv[k]);
			for (--k ; k >= 0 ; k--)
			{
				// here we should add to 'adr' the product of :
				//		#-of-values-of-variable-k * current-value-of-k
				int value = parents[k];
				assert(value >= 0);

				addr += j*value ;
				j *= getStaticDomainSize(argv[k]) ;
			}

			if( ((CConstraintTable*)fun)->getValueAt(addr) == 1)
			{
				isConsistent = true;
				break;
			}
		}

		if(parents)
			delete[] parents;

		if (isAssigned)
				delete[] isAssigned;

		if (!isConsistent)
		{			
			// set old value back
			m_assignment[var] = old_value;
			return false;
		}
	}

	// set old value back
	m_assignment[var] = old_value;
	return true;
}





bool CProblemMixed::load(int itype, char* filename)
{
	switch (itype)
	{
	case INP_SIMPLE:
		{
			// Declare variables.
			int argc,i;
			int* argv = NULL;
			CProbabilityTable* cpt = NULL;
			
			// Read in BN file.
			FILE     *bn_file;
			int LINE_LEN = 100000, num_vars, num_funs, *scopes, *domains;
			char     line[100000], nextc, tmp_string[1000];
			int      j, var, fun, var_num;
			map<int, VARIABLE* > variables;
			int table_mode = MSB;

			bn_file = fopen(filename, "r");
			assert(bn_file != NULL);

			// Skip the comments.
			do{
				fgets(line, LINE_LEN, bn_file);
				sscanf(line, "%c", &nextc );
			} while ( nextc == 'c' );
			//First line after comments is the network name, skip this for now as well.
			
			fgets(line, LINE_LEN, bn_file);
			sscanf(line, "%c %s", &nextc, &tmp_string);
			if (!strcmp(tmp_string, "lsb"))
				table_mode = LSB;
			else if (!strcmp(tmp_string, "msb"))
				table_mode = MSB;
			
			// Read num_vars.
			fgets(line, LINE_LEN, bn_file);
			if( !sscanf(line, "%i", &num_vars ))
			{
				printf("Problem has to start with #variables.\n");
				printf("WRONG INPUT FORMAT, ABORT !!\n");
				exit(-1);
			}
			
			// Update number of variables.
			num_funs = num_vars;

			m_N = num_vars;
			m_K = -1; // initial maximal domain size, will be kept track of.

			scopes = new int[num_vars];
			domains = new int[num_vars];

			// Read optimal MPE value if specified.
			int c = getc(bn_file);
			int tmp;
			ungetc(c, bn_file);
			if (c == 'l')
			{   // log optimal mpe value specified.
				fgets(line, LINE_LEN, bn_file);
				sscanf(line, "log10_opt_mpe = %lf", &tmp);
			}

			// Read variables.
			for (var = 0; var < num_vars; ++var)
			{
				fgets(line, LINE_LEN, bn_file);
				int num_parents;
				if (sscanf(line, "%c %d %d %d %s", &nextc, &var_num, &domains[var], &num_parents, tmp_string) < 5)
				{
					printf("Reading variable line failed:\n%s\nIt has to consist of a v, three integers (var_num, dom.size, num. parents) and a string (varname).\n", line );
					printf("WRONG INPUT FORMAT, ABORT !!\n");
					exit(-1);
				}
				
				// Safety checks.
				assert(nextc == 'v');	
				
				// Record new variable.
				VARIABLE* pvar = new VARIABLE(var, var_num, domains[var]);
				variables.insert(make_pair(var_num, pvar));

				// Keep track of maximal domain size so far.
				m_K = max(m_K, domains[var]);  
				scopes[var] = num_parents + 1;
			}

			init(); // need m_N and m_K for this. Needs to be done before constructing CPTs

			// Update domain sizes to all variables.
			map<int, VARIABLE* >::iterator itvar = variables.begin();
			for (; itvar != variables.end(); ++itvar)
			{
				VARIABLE* pvar = (*itvar).second;

				m_domains[pvar->index] = pvar->domain;
				//assert(pvar->domain == domain_sizes[pvar->index]);
			}

			// Read probabilities (potentials).
			for (fun = 0; fun < num_funs; ++fun)
			{
				//==== Read variable.
				fgets(line, LINE_LEN, bn_file);
				int chars_read, tmp;
				int child;

				// Create the scope of the function.
				switch (table_mode)
				{
				case MSB:
					{
						sscanf(line, "%c %d %n", &nextc, &var_num, &chars_read);
						assert(nextc == 'p');

						// Child is first position.
						map<int, VARIABLE* >::iterator itvar = variables.find(var_num);
						VARIABLE* pvar = (*itvar).second;
						child = pvar->index;

						//==== Read parents of variable.
						argc = scopes[child];
						argv = new int[argc];

						for (j = 1; j < argc; ++j)
						{
							int par_num;
							sscanf(line+chars_read, "%d %n", &par_num, &tmp);
							chars_read += tmp;
							
							itvar = variables.find(par_num);
							pvar = (*itvar).second;

							argv[j - 1] = pvar->index;
						}

						// Set the child of the potential.
						argv[argc - 1] = child;

						break;
					}
				case LSB:
					{
						char seps[] = " \t\n";
						char* tok = strtok(line, seps);
						assert(0 == strcmp(tok, "p"));
						
						vector<int> v_argv;
						tok = strtok(NULL, seps);
						while (tok)
						{
							int v = atoi(tok);
							v_argv.push_back(v);

							tok = strtok(NULL, seps);
						}

						// Set function scope. Child is lsb.
						argc = v_argv.size();
						argv = new int[argc];
						
						j = 0;
						vector<int>::iterator it = v_argv.begin();
						for (; it != v_argv.end(); ++it)
							argv[j++] = (*it);

						v_argv.clear();

						child = argv[argc - 1];

						break;
					}
				};

				// Create probability table (CPT).
				cpt = new CProbabilityTable(0, CPT, argc, argv);
				cpt->setOwner(this);

				// Compute number of entries (table size).
				int num_entries = 1;
				for(i = 0; i < argc; i++)
				{
					num_entries *= domains[argv[i]]; // for this, the domain sizes must be correct.
				}

				// Allocate memory for the potential.
				double *table = new double[num_entries];
				assert(table != NULL);

				// Read probability table.
				switch (table_mode)
				{
				case MSB:
					{
						chars_read = 0;
						int factor_of_child = num_entries/domains[child];
						for (j = 0; j < num_entries; ++j)
						{
							//== In .simple format, child is the most significant bit. Here, it is the 
							//== least significant one. Thus, convert the index.
							int value_of_child = j / factor_of_child;
							int index = (j-(factor_of_child*value_of_child))*domains[child]+value_of_child;

							double p;
							fscanf(bn_file, "%lf %n", &p, &tmp);
							table[index] = p;
							chars_read += tmp;
						}

						break;
					}
				case LSB:
					{
						chars_read = 0;
						for (j = 0; j < num_entries; ++j)
						{
							//== In .simple format, child is the most significant bit. Here, it is the 
							//== least significant one. Thus, convert the index.
							int index = j;
							double p;

							fscanf(bn_file, "%lf %n", &p, &tmp);
							table[index] = p;
							chars_read += tmp;
						}
					}
				};

				// Create the table (do not normalize).
				cpt->create(num_entries, table);

				// Add function to the repository.
				addFunction(cpt);
				delete[] table;
			}

			// Create the moral graph.
			m_graph->init(m_functions);
			m_connected = m_graph->isConnected();

			// Clear memory.
			map<int, VARIABLE* >::iterator itvars = variables.begin();
			for (; itvars != variables.end(); ++itvars)
			{
				VARIABLE* pvar = (*itvars).second;
				delete pvar;
			}
			variables.clear();
			delete[] domains;
			delete[] scopes;
			

			break;
		}
	default:
		break;
	}

	return true;
}


bool clst_desc(CTreeDecompositionCluster* a, CTreeDecompositionCluster* b) // a < b if highestVar(a) < highestVar(b)
{
	return (a->highestVar() < b->highestVar());
}

void CProblemMixed::makeTreeDecomposition(int type)
{
//	if (m_treeDecomposition)
//		delete (m_treeDecomposition);

//	m_treeDecomposition = new CTreeDecomposition(this);

	treeDecompositionCluster_v clusters;
	makeTreeDecompositionClusters(type, clusters);
	m_treeDecomposition = makeTreeDecompositionClustersAssemble(clusters);
	
	
	//debug
	// cout << "lllll";
	// m_treeDecomposition->print();
	// tempTree->print();

}

bool CProblemMixed::load2(int itype, char* filename)
{
	switch (itype)
	{
	case INP_SIMPLE:
		{
			// Declare variables.
			int argc,i;
			int* argv = NULL;
			CProbabilityTable* cpt = NULL;
			
			// Read in BN file.
			FILE     *bn_file;
			int LINE_LEN = 100000, num_vars, num_funs, *scopes, *domains;
			char     line[100000], nextc, tmp_string[1000];
			int      j, var, fun, var_num;
			map<int, VARIABLE* > variables;
			int table_mode = MSB;

			bn_file = fopen(filename, "r");
			assert(bn_file != NULL);

			// Skip the comments.
			do{
				fgets(line, LINE_LEN, bn_file);
				sscanf(line, "%c", &nextc );
			} while ( nextc == 'c' );
			//First line after comments is the network name, skip this for now as well.
			
			fgets(line, LINE_LEN, bn_file);
			sscanf(line, "%c %s", &nextc, &tmp_string);
			if (!strcmp(tmp_string, "lsb"))
				table_mode = LSB;
			else if (!strcmp(tmp_string, "msb"))
				table_mode = MSB;
			
			// Read num_vars.
			fgets(line, LINE_LEN, bn_file);
			if( !sscanf(line, "%i", &num_vars ))
			{
				printf("Problem has to start with #variables.\n");
				printf("WRONG INPUT FORMAT, ABORT !!\n");
				exit(-1);
			}
			
			// Update number of variables.
			num_funs = num_vars;

			m_N = num_vars;
			m_K = -1; // initial maximal domain size, will be kept track of.
			m_C = num_funs;

			scopes = new int[num_vars];
			domains = new int[num_vars];

			// Read optimal MPE value if specified.
			int c = getc(bn_file);
			int tmp;
			ungetc(c, bn_file);
			if (c == 'l')
			{   // log optimal mpe value specified.
				fgets(line, LINE_LEN, bn_file);
				sscanf(line, "log10_opt_mpe = %lf", &tmp);
			}

			// Read variables.
			for (var = 0; var < num_vars; ++var)
			{
				fgets(line, LINE_LEN, bn_file);
				int num_parents;
				if (sscanf(line, "%c %d %d %d %s", &nextc, &var_num, &domains[var], &num_parents, tmp_string) < 5)
				{
					printf("Reading variable line failed:\n%s\nIt has to consist of a v, three integers (var_num, dom.size, num. parents) and a string (varname).\n", line );
					printf("WRONG INPUT FORMAT, ABORT !!\n");
					exit(-1);
				}
				
				// Safety checks.
				assert(nextc == 'v');	
				
				// Record new variable.
				VARIABLE* pvar = new VARIABLE(var, var_num, domains[var]);
				variables.insert(make_pair(var_num, pvar));

				string key = string(tmp_string);
				m_mapInt2String[var] = key;
				m_mapString2Int[key] = var;

				// Keep track of maximal domain size so far.
				m_K = max(m_K, domains[var]);  
				scopes[var] = num_parents + 1;
			}

			init();

			// Update domain sizes to all variables.
			map<int, VARIABLE* >::iterator itvar = variables.begin();
			for (; itvar != variables.end(); ++itvar)
			{
				VARIABLE* pvar = (*itvar).second;

				m_domains[pvar->index] = pvar->domain;
			}

			// Read probabilities (potentials).
			for (fun = 0; fun < num_funs; ++fun)
			{
				//==== Read variable.
				fgets(line, LINE_LEN, bn_file);
				int chars_read, tmp;
				int child;

				// Create the scope of the function.
				switch (table_mode)
				{
				case MSB:
					{
						sscanf(line, "%c %d %n", &nextc, &var_num, &chars_read);
						assert(nextc == 'p');

						// Child is first position.
						variable_m::iterator itvar = variables.find(var_num);
						VARIABLE* pvar = (*itvar).second;
						child = pvar->index;

						//==== Read parents of variable.
						argc = scopes[child];
						argv = new int[argc];

						for (j = 1; j < argc; ++j)
						{
							int par_num;
							sscanf(line+chars_read, "%d %n", &par_num, &tmp);
							chars_read += tmp;
							
							itvar = variables.find(par_num);
							pvar = (*itvar).second;

							argv[j - 1] = pvar->index;
						}

						// Set the child of the potential.
						argv[argc - 1] = child;

						break;
					}
				case LSB:
					{
						char seps[] = " \t\n";
						char* tok = strtok(line, seps);
						assert(0 == strcmp(tok, "p"));
						
						vector<int> v_argv;
						tok = strtok(NULL, seps);
						while (tok)
						{
							int v = atoi(tok);
							v_argv.push_back(v);

							tok = strtok(NULL, seps);
						}

						// Set function scope. Child is lsb.
						argc = v_argv.size();
						argv = new int[argc];
						
						j = 0;
						vector<int>::iterator it = v_argv.begin();
						for (; it != v_argv.end(); ++it)
						{
							var_num = (*it);
							variable_m::iterator itvar = variables.find(var_num);
							VARIABLE* pvar = (*itvar).second;
							var = pvar->index;

							argv[j++] = var;
						}

						v_argv.clear();

						// Child is last in scope.
						child = argv[argc - 1];

						break;
					}
				};

				// Create probability table (CPT).
				unsigned long id = (unsigned long)child;
				cpt = new CProbabilityTable(id, CPT, argc, argv);
				cpt->setOwner(this);

				// Compute number of entries (table size).
				int num_entries = 1;
				for(i = 0; i < argc; i++)
				{
					num_entries *= domains[argv[i]]; // for this, the domain sizes must be correct.
				}

				// Allocate memory for the potential.
				double *table = new double[num_entries];
				assert(table != NULL);

				// Read probability table.
				switch (table_mode)
				{
				case MSB:
					{
						chars_read = 0;
						int factor_of_child = num_entries/domains[child];
						for (j = 0; j < num_entries; ++j)
						{
							//== In .simple format, child is the most significant bit. Here, it is the 
							//== least significant one. Thus, convert the index.
							int value_of_child = j / factor_of_child;
							int index = (j-(factor_of_child*value_of_child))*domains[child]+value_of_child;

							double p;
							fscanf(bn_file, "%lf %n", &p, &tmp);
							table[index] = p;
							chars_read += tmp;
						}

						break;
					}
				case LSB:
					{
						chars_read = 0;
						for (j = 0; j < num_entries; ++j)
						{
							//== In .simple format, child is the most significant bit. Here, it is the 
							//== least significant one. Thus, convert the index.
							int index = j;
							double p;

							fscanf(bn_file, "%lf %n", &p, &tmp);
							table[index] = p;
							chars_read += tmp;
						}
					}
				};

				// Create the table (do not normalize).
				cpt->create(num_entries, table);
				delete[] table;

				// Add function to the repository.
				addFunction(cpt, false);
			}

			// Create the moral graph.
			m_graph->init(m_functions);
			m_connected = m_graph->isConnected();

			// Clear memory.
			map<int, VARIABLE* >::iterator itvars = variables.begin();
			for (; itvars != variables.end(); ++itvars)
			{
				VARIABLE* pvar = (*itvars).second;
				delete pvar;
			}
			variables.clear();
			delete[] domains;
			delete[] scopes;

			break;
		}

	case INP_ORDER:
		{
			// Declare variables.
			int i, var;

			// Read in BN file.
			FILE     *fp;
			int LINE_LEN = 1024;
			char     line[1024], s1[128], s2[128];

			fp = fopen(filename, "r");
			assert(fp != NULL);

			// Parse the order file.
			while (!feof(fp))
			{
				s1[0] = '\0';
				s2[0] = '\0';
				fgets(line, LINE_LEN, fp);
				sscanf(line, "%d %s %s", &i, &s1, &s2);

				string key = string(s1);
				map<string, int>::iterator it;
				it = m_mapString2Int.find(key);
				if (it != m_mapString2Int.end())
				{
                    var = (*it).second;
				
					if (find(m_superlinkOrder.begin(), m_superlinkOrder.end(), var) != m_superlinkOrder.end())
						continue;

					m_superlinkOrder.push_back(var);
					if (0 == strcmp(s2, "breakNode"))
						m_superlinkCutset.push_back(var);
				}
			}
			
			int ss1 = (int)m_superlinkOrder.size();
			int ss2 = (int)m_superlinkCutset.size();

			// Note: The superlink order does not contain the evidence
			// variables (namely singleton variables).

			break;
		}
	default:
		break;
	}

	return true;
}


bool CProblemMixed::parser_simple2(char *filename)
{
	// Declare variables.
	char buf[100000], lsb[4], ev[4], c, name[128], *tok;
	int sz = 100000, ne, n, i, id, d, e, np, k = 0, r = 0, var, ch;
	bool flag = true;
	vector<int> domains, parents, evidence, scope;
	vector<double> table;
	map<int,int> id2var;

	// Read in BN file.
	ifstream in(filename);
	//m_id2name.clear();

	// Prologue.
	in.getline(buf, sz);	// network
	in.getline(buf, sz);	// format
	sscanf(buf, "%c %s", &c, &lsb);
	in.getline(buf, sz);	// evidence
	sscanf(buf, "%c %s %d", &c, &ev, &ne);
	if (!ne) flag = false;

	// Variables
	in.getline(buf, sz);
	sscanf(buf, "%d", &n);
	domains.resize(n, 0);
	parents.resize(n, 0);
	evidence.resize(n, -1);
	for (i = 0; i < n; ++i)
	{
		in.getline(buf, sz);
		if (flag) sscanf(buf, "%c %d %d %d %d %s", &c, &id, &d, &e, &np, &name);
		else sscanf(buf, "%c %d %d %d %s", &c, &id, &d, &np, &name);

		id2var[id] = i;
		domains[i] = d;
		parents[i] = np;
		evidence[i] = (flag) ? e : -1;
		//m_id2name[i] = string(name);

		k = max(k, d);
		r = max(r, np+1);
	}

	m_N = n;
	m_K = k;
	m_r = r;
	init();

	// Update domain sizes.
	for (size_t ii = 0; ii < domains.size(); ++ii)
		m_domains[ii] = domains[ii];

	// Functions
	for (i = 0; i < n; ++i)
	{
		in.getline(buf, sz);
		tok = strtok(buf, " ");
		assert(0 == strcmp(tok, "p"));
	
		// Read in the scope.
		scope.clear();
		tok = strtok(NULL, " ");
		while (tok)
		{
			id = atoi(tok);
			map<int,int>::iterator mi = id2var.find(id);
			var = (*mi).second;

			scope.push_back(var);
			tok = strtok(NULL, " ");
		}
		
		ch = scope.back();
		assert(scope.size() == parents[ch]+1);

		// Create a function and add it to the network.
		int argc = (int)scope.size();
		int *argv = new int[argc];
		int tab_size = 1;
		for (size_t ii = 0, j = 0; ii < scope.size(); ++ii)
		{
			tab_size *= domains[scope[ii]];
			argv[j++] = scope[ii];
		}
		// Read in the table.
		table.resize(tab_size, 0);
		for (int ti = 0; ti < tab_size; ++ti) in >> table[ti];
		in.getline(buf, sz);

		assert(tab_size == (int)table.size());
		double *tab = new double[tab_size];
		for (size_t ii = 0, j = 0; ii < table.size(); ++ii)
			tab[j++] = table[ii];

		
		CProbabilityTable *fn = new CProbabilityTable(ch, CPT, argc, argv, tab_size, tab);
		fn->setOwner(this);
		addFunction(fn);
	}

	// Update the number of functions.
	m_C = (int)m_functions.size();

	// Set the evidence, if any.
	for (i = 0; i < n; ++i)
		if (evidence[i] >= 0)
			setEvidence(i, evidence[i]);

	// Free memory.
	table.clear();
	domains.clear();
	scope.clear();
	parents.clear();
	evidence.clear();
	in.close();

	return true;
}

// mark the connected components of the problems and pack each into a subproblem
void CProblemMixed::makeConnectedComponents()
{
	// Assert as evidence all singleton variables. Reindex variables, scopes, moral graph.
	// removeEvidence();

	int connectedComponent;

	// cout << "\n number of variables: " << m_N;

	// Create new moral graph.
	m_graph->init(m_functions);
	m_connected = m_graph->isConnected();

	//// debug
	//function_map::iterator it = m_functions.begin();
	//for( ; it != m_functions.end(); ++it)
	//{
	//	CFunction* func = (*it).second;
	//	int argc = func->getArgc();
	//	int* argv = func->getArgv();
	//	for (int i=0; i<argc; i++)
	//		cout << argv[i] << " ";
	//	cout << endl;
	//}

	//if(m_graph->isConnected())
	//{
	//	independentComponentProblems.push_back(this);
	//	cout << "\n graph has only one component";
	//}
	//else
	//{
		//make a dfs traversal to mark the connected components
		int* visited = new int[m_N];
		for(int i=0; i<m_N; i++)
			visited[i] = -1;			

		stack<int> dfsStack;
		connectedComponent = -1;

		for (int v=0; v < m_N; v++)
		{
			if (visited[v] != -1)
				continue;
			else
			{
				connectedComponent++;
				// cout << "\n component number: " << connectedComponent;
				visited[v] = connectedComponent;
				dfsStack.push(v);

				while (!dfsStack.empty())
				{
					int var = dfsStack.top();
					dfsStack.pop();

					// Get the not yet visited neighbors.
					variable_s neighbors;
					m_graph->neighborhood(var, neighbors);

					variable_s::iterator it = neighbors.begin();
					for (; it != neighbors.end(); ++it)
					{
						int n = (*it);
						if (visited[n] == -1)
						{
							visited[n] = connectedComponent;
							dfsStack.push(n);
						}
					}

					neighbors.clear();
				}
			}
		}

		// each component is a vector with its variables
		variable_v* components = new variable_v[connectedComponent+1];
		for (int v=0 ; v < m_N; v++)
		{
			int c = visited[v];
			components[c].push_back(v);
		}

		// debug
		//for (int c = 0; c < connectedComponent+1; c++)
		//{
		//	cout << "\n component " << c << "; contains " << components[c].size() << " variables " << endl ;

		//	variable_v::iterator it = components[c].begin();
		//	for ( ; it != components[c].end(); ++it )
		//		cout << (*it) << " ";			
		//}

		// mark the functions for the creation of components
		function_map::iterator itf = m_functions.begin();
		for ( ; itf != m_functions.end(); ++itf )
		{
			CFunction* f = (*itf).second;
			f->setOriginal(false);
		}

		// now create a new problem for each component
		for (int c = connectedComponent; c >= 0; c--)
		{
			CProblemMixed* newProb = new CProblemMixed();
			
			// debug
			// cout << "\n now processing component number: " << c;

			newProb->createProblemFromComponent(components[c], m_functions, m_domains, visited, c);
			independentComponentProblems.push_back(newProb);
		}
	
		m_functions.clear();
	//}
}

// given a connected component of a bigger problem, create a new smaller problem
void CProblemMixed::createProblemFromComponent(variable_v& vars, function_map& functions, int* domains, int* visited, int componentNumber)
{
	map<int,int> mapOld2New;
	function_v newFuns;

	// set number of nodes
	m_N = vars.size();

	// Initialize the problem.
	init();

	int newK = 0;
	for (int i=0; i<m_N; i++)
	{
		mapOld2New[vars[i]] = i;
		m_domains[i] = domains[vars[i]];
		newK = max(newK, m_domains[i]);
	}
	// set K
	m_K = newK;	

	int fCount = 0;
	function_map::iterator it = functions.begin();
	for ( ; it != functions.end(); ++it )
	{
		CFunction* f = (*it).second;
		if( f->isOriginal() )
			continue;

		int* argv = f->getArgv();
		int firstVar = argv[0];
		if (visited[firstVar] != componentNumber)
			continue; // the function does not belong to current component\

		int argc = f->getArgc();
		
		// debug
		// cout << "\n function number: " << fCount << "     last var: " << argv[argc-1];
				
		f->setOriginal(true);
		f->setID(fCount);

		// Reindex the scope.
		for (int i = 0; i < argc; ++i)
		{
			int ovar = argv[i];
			map<int,int>::iterator itv = mapOld2New.find(ovar);
			assert(itv != mapOld2New.end());
			int nvar = (*itv).second;

			argv[i] = nvar;
		}
		
		f->setOwner(this);
		addFunction(f);
		fCount++;
	}
	
	m_C = (int)m_functions.size();
	
	mapOld2New.clear();	
}



void CProblemMixed::makeTreeDecompositionClusters(int type, treeDecompositionCluster_v &clusters)
{
	m_graph->order(type, m_ordering, m_position);		// order the graph

	// create the induced graph (induced edges = 3, initial =1)
	int width = m_graph->width(m_ordering, m_position);							

	// now create all the clusters (one for each variable: var and parents)

 	for (int i = m_N - 1; i >= 1; --i)
	{
		int var = m_ordering[i];
		variable_s clst_vars;

		for (int j = i - 1; j >= 0; --j)
		{
			int parent = m_ordering[j];
			if (m_graph->isConnected(var, parent))
			{
				clst_vars.insert(parent);
			}
		}
		clst_vars.insert(var);
		CTreeDecompositionCluster* clst = new CTreeDecompositionCluster(var, clst_vars);
		clst->setOwnerProb(this);
		

		clusters.push_back(clst);
		clst_vars.clear();
	}


	// //temp debug
	// cout << "\nClusters before reversal:\n";
	// treeDecompositionCluster_v::reverse_iterator ttt = clusters.rbegin();
	// for ( ; ttt != clusters.rend(); ttt++)
	//	(*ttt)->print();
	// //end temp debug

	// // reverse clusters, so that it is ordered according to ordering (induced from last to first)
	reverse(clusters.begin(), clusters.end());
	
/*
	// now pick only maximal clusters from first in ordering to last	
	treeDecompositionCluster_v temp_clusters;
	treeDecompositionCluster_v::iterator it1 = clusters.begin();
	treeDecompositionCluster_v::iterator it2;
	for ( ; it1 != clusters.end() ; it1++)
	{
		bool isMaximal = true;
		for ( it2 = it1 + 1; it2 != clusters.end(); it2++)
		{			
			variable_s vars1 = (*it1)->getVariables();
			variable_s vars2 = (*it2)->getVariables();
		
			if( includes(vars2.begin(), vars2.end(), vars1.begin(), vars1.end()) )
			{
				isMaximal = false;
				break;
			}
		}

		if(isMaximal)
			temp_clusters.push_back(*it1);
		else
			delete(*it1);
	}

	clusters.clear();
	copy(temp_clusters.begin(), temp_clusters.end(), back_inserter(clusters));
	temp_clusters.clear();
*/

	// temp debug
	// cout << "\nClusters before sorting:\n";
	// treeDecompositionCluster_v::reverse_iterator tt1 = clusters.rbegin();
	// for ( ; tt1 != clusters.rend(); tt1++)
	//	(*tt1)->print();
	// end temp debug
	
	// sort clusters
	sort(clusters.begin(), clusters.end(), clst_desc);

	 //temp debug
	 //cout << "\nClusters after sorting:\n";
	 //treeDecompositionCluster_v::reverse_iterator tt = clusters.rbegin();
	 //for ( ; tt != clusters.rend(); tt++)
	//	(*tt)->print();
	 //end temp debug
}



// assembles the clusters in a tree decomposition, computes separators
CTreeDecomposition* CProblemMixed::makeTreeDecompositionClustersAssemble(treeDecompositionCluster_v &clusters)		
{
	CTreeDecomposition* treeDecomposition = new CTreeDecomposition(this);

//	sort(clusters.begin(), clusters.end(), clst_desc);
	
	treeDecompositionCluster_v::iterator it = clusters.begin();
	for ( ; it != clusters.end() ; it ++ )
	{
		(*it)->setOwner(treeDecomposition);
		(*it)->setOwnerProb(this);
		(*it)->children().clear();
		(*it)->setParent(NULL);
	}

	// assemble join tree and form separators
	// connect each cluster to a predecessor in the ordering, with which it shares the max number of variables	
	treeDecomposition->setRootCluster(*(clusters.begin()));
	if (clusters.size() > 1)
	{		
		// start from end of ordering to beginning
		treeDecompositionCluster_v::reverse_iterator cl_child = clusters.rbegin();
		treeDecompositionCluster_v::reverse_iterator cl_parent;
		CTreeDecompositionCluster* best_parent;		
		variable_s sharedVars, bestSharedVars;
		for ( ; cl_child != clusters.rend() ; cl_child++ )
		{			
			// find earlier cluster with which cl_child shares max number of variables			
			best_parent = NULL;
			bestSharedVars.clear();
			sharedVars.clear();

			if ( (*cl_child) == (*clusters.begin()) )
				break;

			for ( cl_parent = cl_child + 1 ; cl_parent != clusters.rend() ; cl_parent++ )
			{
				//copy((*cl_parent)->getVariables().begin(), (*cl_parent)->getVariables().end(), ostream_iterator<int>(cout, " "));
				//cout << "\n";
				//copy((*cl_child)->getVariables().begin(), (*cl_child)->getVariables().end(), ostream_iterator<int>(cout, " "));
				//cout << "\n";
				set_intersection(	(*cl_child)->getVariables().begin(), (*cl_child)->getVariables().end(),
									(*cl_parent)->getVariables().begin(), (*cl_parent)->getVariables().end(),
									inserter(sharedVars, sharedVars.begin()) );
				if ( sharedVars.size() > bestSharedVars.size() ) 
				{
					best_parent = (*cl_parent);
					bestSharedVars.clear();
					copy(sharedVars.begin(), sharedVars.end(), inserter(bestSharedVars, bestSharedVars.begin()));
				}
				sharedVars.clear();
			}


			// if we can't connect to lower cluster, then try to connect up to one which is not a child
			if(best_parent == NULL)
			{
				bestSharedVars.clear();
				sharedVars.clear();

				for (cl_parent = clusters.rbegin() ; cl_parent != cl_child ; cl_parent++)
				{
					// verify all ancestors of cl_parent, make sure they are different than cl_child
					CTreeDecompositionCluster* pClst = (*cl_parent);
					bool isAncestor = false;
					while(pClst != NULL)
					{
						if( pClst->parent() == (*cl_child) )
						{
							isAncestor = true;
							break;
						}

						pClst = pClst->parent();
					}

					if ( isAncestor )
						continue;

					set_intersection(	(*cl_child)->getVariables().begin(), (*cl_child)->getVariables().end(),
										(*cl_parent)->getVariables().begin(), (*cl_parent)->getVariables().end(),
										inserter(sharedVars, sharedVars.begin()) );
					if ( sharedVars.size() > bestSharedVars.size() ) 
					{
						best_parent = (*cl_parent);
						bestSharedVars.clear();
						copy(sharedVars.begin(), sharedVars.end(), inserter(bestSharedVars, bestSharedVars.begin()));
					}
					sharedVars.clear();
				}				
			}

		
			assert(best_parent != NULL);

			// make parent child links
			best_parent->addChild( (*cl_child) );
			(*cl_child)->setParent(best_parent);

			// make separator
			CTreeDecompositionSeparator *sep = new CTreeDecompositionSeparator(bestSharedVars);
			sep->setClusterChild(*cl_child);
			sep->setClusterParent(best_parent);
			sep->setClusterIDChild( (*cl_child)->getID() );
			sep->setClusterIDParent( best_parent->getID() );
			sep->setOwner(treeDecomposition);
			treeDecomposition->getSeparators().push_back(sep);

			bestSharedVars.clear();

			// cout << "[ c=" << (*cl_child)->getID() << " --- p=" << best_parent->getID() << " ]\n";
		}	
	}

	// make the tree decomposition
	treeDecomposition->setClusters(clusters);
	treeDecomposition->setSize(clusters.size());
	clusters.clear();
	
	// debug
	// cout << "\n\n test cyclicity -----------------------------------------------********************--\n";
	// treeDecomposition->print();


	// if ( treeDecomposition->cyclic() )
	// {
	// 	cout << "\n\n -------------------------------------------------------------cyclic-----------";
	// 	exit(0);
	// }
	// else
	// 	cout << "\n\n -------------------------------------------------------NOT---cyclic-----------";


	return treeDecomposition;
}


void CProblemMixed::makeAndOrCutsetTree(int scopeSize)
{	
	assert(m_tree_cutset == NULL);

	// make a copy of m_treeDecomposition, because it will be altered in the next step
	treeDecompositionCluster_v clusters;
	treeDecompositionCluster_v origClusters = m_treeDecomposition->getClusters();
	treeDecompositionCluster_v::iterator it = origClusters.begin();
	for ( ; it != origClusters.end() ; it++ )
	{
		CTreeDecompositionCluster* tempClust = new CTreeDecompositionCluster( (*it)->getID(), (*it)->getVariables() );
		clusters.push_back(tempClust);	
	}
	
	CTreeDecomposition* startTreeDecomp = makeTreeDecompositionClustersAssemble(clusters);

		
	// debug
	// cout << "\n\n Start tree decomp ---------------------------\n";
	// startTreeDecomp->print();

	// this creates all the nodes of the pseudo tree, with the old ordering and position!
	m_tree_cutset = new CLegalTree(m_N, m_ordering, m_position);
	
	// create dummy root node
	CLegalTreeNode* dummyRoot = new CLegalTreeNode(-1);
	dummyRoot->setHeight(-1);
	startTreeDecomp->setCutsetTreeNode(dummyRoot);

	//debug
	m_tree_cutset->setRoot(dummyRoot);

	treeDecomposition_v treeDecStack;
	treeDecStack.push_back(startTreeDecomp);

	while (!treeDecStack.empty())
	{
		CTreeDecomposition* curTreeDec = treeDecStack.back();
		treeDecStack.pop_back();


		// debug
		// cout << "\n\n ---------------------------------------------------------\n";
		// curTreeDec->print();
		// m_tree_cutset->print();

		// if the curTreeDec is made of just one cluster		
		if (curTreeDec->getClusters().size() == 1)
		{
			int size = ( * curTreeDec->getClusters().begin() )->getVariables().size();
			int* remWidth = new int[size];
			for (int i = 0 ; i < size ; i++)
				remWidth[i] = size - i - 1 ;

			variable_v vars;

			variable_s origVars = ( *curTreeDec->getClusters().begin() )->getVariables();

			copy(origVars.begin(), origVars.end(), inserter(vars, vars.begin() ) );

			m_tree_cutset->addVariables(curTreeDec->getCutsetTreeNode(), vars, remWidth );

			//if (curTreeDec != m_treeDecomposition)
				delete curTreeDec;
					
			//debug
			//m_treeDecomposition->print();

			vars.clear();
			delete[] remWidth;
			continue;
		}

		curTreeDec->computeSeparatorWeights(scopeSize);

		CTreeDecompositionSeparator* breakSeparator = curTreeDec->highestWeightSeparator();
		
		variable_v breakVars;		
		copy(breakSeparator->getVariables().begin(), breakSeparator->getVariables().end(), inserter(breakVars, breakVars.begin()));

		// order the separator vars; the one that leaves min remaining width towards .begin()
		variable_v breakVarsOrdered;
		int* remainingWidth = new int[breakVars.size()];

		int var = -1;

		int count = 0;

		// let's try eliminating one variable at a time (remove while loop)
		while(breakVars.size() > 0)
		{
			variable_v::iterator it = breakVars.begin();
			int minWidth = INT_MAX;
			for ( ; it != breakVars.end() ; it++ )
			{
				int tempWidth = curTreeDec->computeRemainingWidth( (*it) );
				if (tempWidth < minWidth )
				{
					var = (*it);
					minWidth = tempWidth;
				}
			}


			breakVarsOrdered.push_back(var);
			remainingWidth[count] = minWidth;

			curTreeDec->eliminateVar(var);

			it = find(breakVars.begin(), breakVars.end(), var);
			assert ( it != breakVars.end() );

			
			breakVars.erase(it);

			count ++;
		}

		
		breakVars.clear();

		// add breakSeparator->getVariables() to the pseudo-tree, increase nodeCounts each time
		// add the var that leaves smallest width first.
		CLegalTreeNode* latestNode = m_tree_cutset->addVariables(curTreeDec->getCutsetTreeNode(), 
																	breakVarsOrdered, remainingWidth);

		// debug
		// cout << "\n\n -------------------------after adding to pseudo tree--------------------------------\n";
		// curTreeDec->print();
		// m_tree_cutset->print();

		// now the tree decomposition already eliminated the vars 
		// mark connected componnents on m_visited field
		int compNum = curTreeDec->connectedComponents();
		
		for (int i = 0; i < compNum; i++)
		{
			treeDecompositionCluster_v tempClust;
			curTreeDec->getComponent( i, tempClust );

			CTreeDecomposition* newTreeDec = makeTreeDecompositionClustersAssemble(tempClust);
		
			newTreeDec->setCutsetTreeNode(latestNode); //set pointer to LegalTreeNode
			
			treeDecStack.push_back(newTreeDec);
		}

		if (curTreeDec != m_treeDecomposition)
			delete curTreeDec;

		//debug
		//m_treeDecomposition->print();

		delete[] remainingWidth;


	}

	m_tree_cutset->setRoot( *dummyRoot->children().begin() );
	delete dummyRoot;



//		m_descendants = new variable_v[m_N];
//		assert(m_descendants != NULL);
//		for (int var = 0; var < m_N; ++var)	
//			m_tree_cutset->descendants(var, m_descendants[var]);

}

void CProblemMixed::makeGWCutsetTree(int scopeSize)
{
	assert(m_tree_cutset == NULL);

	// make a copy of m_tree_cutset, because it will be altered in the next step
	treeDecompositionCluster_v clusters;
	treeDecompositionCluster_v origClusters = m_treeDecomposition->getClusters();
	treeDecompositionCluster_v::iterator it = origClusters.begin();
	for ( ; it != origClusters.end() ; it++ )
	{
		CTreeDecompositionCluster* tempClust = new CTreeDecompositionCluster( (*it)->getID(), (*it)->getVariables() );
		clusters.push_back(tempClust);	
	}
	
	CTreeDecomposition* startTreeDecomp = makeTreeDecompositionClustersAssemble(clusters);
		
	// debug
	// cout << "\n\n Start tree decomp ---------------------------\n";
	// startTreeDecomp->print();

	// this creates all the nodes of the pseudo tree, with the old ordering and position!
	m_tree_cutset = new CLegalTree(m_N, m_ordering, m_position);
	
	// create dummy root node
	CLegalTreeNode* dummyRoot = new CLegalTreeNode(-1);
	dummyRoot->setHeight(-1);
	startTreeDecomp->setCutsetTreeNode(dummyRoot);

	//debug
	m_tree_cutset->setRoot(dummyRoot);

	treeDecomposition_v treeDecStack;
	treeDecStack.push_back(startTreeDecomp);

	while (!treeDecStack.empty())
	{
		CTreeDecomposition* curTreeDec = treeDecStack.back();
		treeDecStack.pop_back();


		// debug
		// cout << "\n\n ---------------------------------------------------------\n";
		// curTreeDec->print();
		// m_tree_cutset->print();

		// if the curTreeDec is made of just one cluster
		
		
		if (curTreeDec->getClusters().size() == 1)
		{
			int size = ( * curTreeDec->getClusters().begin() )->getVariables().size();
			int* remWidth = new int[size];
			for (int i = 0 ; i < size ; i++)
				remWidth[i] = size - i - 1 ;

			variable_v vars;

			variable_s origVars = ( *curTreeDec->getClusters().begin() )->getVariables();

			copy(origVars.begin(), origVars.end(), inserter(vars, vars.begin() ) );

			m_tree_cutset->addVariables(curTreeDec->getCutsetTreeNode(), vars, remWidth );

			//if (curTreeDec != m_treeDecomposition)
				delete curTreeDec;
					
			//debug
			//m_treeDecomposition->print();

			vars.clear();
			delete[] remWidth;
			continue;
		}

		// here compute the best variable to be eliminated, according to GWC
		int var = curTreeDec->computeBestVariableGWC(scopeSize);
		variable_v breakVarsOrdered;
		breakVarsOrdered.push_back(var);

		int* remainingWidth = new int[1];
		remainingWidth[0] = curTreeDec->computeRemainingWidth(var);

		curTreeDec->eliminateVar(var);

		// add breakSeparator->getVariables() to the pseudo-tree, increase nodeCounts each time
		// add the var that leaves smallest width first.
		CLegalTreeNode* latestNode = m_tree_cutset->addVariables(curTreeDec->getCutsetTreeNode(), 
																	breakVarsOrdered, remainingWidth);

		// debug
		// cout << "\n\n -------------------------after adding to pseudo tree--------------------------------\n";
		// curTreeDec->print();
		// m_tree_cutset->print();

		// now the tree decomposition already eliminated the vars 
		// mark connected componnents on m_visited field
		int compNum = curTreeDec->connectedComponents();
		
		for (int i = 0; i < compNum; i++)
		{
			treeDecompositionCluster_v tempClust;
			curTreeDec->getComponent( i, tempClust );

			CTreeDecomposition* newTreeDec = makeTreeDecompositionClustersAssemble(tempClust);
		
			newTreeDec->setCutsetTreeNode(latestNode); //set pointer to LegalTreeNode
			
			treeDecStack.push_back(newTreeDec);
		}

		if (curTreeDec != m_treeDecomposition)
			delete curTreeDec;

		//debug
		//m_treeDecomposition->print();

		delete[] remainingWidth;


	}

	m_tree_cutset->setRoot( *dummyRoot->children().begin() );
	delete dummyRoot;
}

void CProblemMixed::makeCutsetOr(int scopeSize)
{

	// m_tree has to be w annotated
//	wCutsetInit();

	CLegalTreeNode* currNode;
	CLegalTreeNode* parentNode;

	currNode = m_tree_cutset->getRoot();
	parentNode = NULL;
	currNode->setParent(NULL);


	stack<CLegalTreeNode*> dfsStack;

	if (currNode->getWCutsetHeight() >= scopeSize)
		dfsStack.push(currNode);

	int count = 0 ;
	while(!dfsStack.empty())
	{		
		currNode = dfsStack.top();
		dfsStack.pop();
	
		// change parent if necessary	
		CLegalTreeNode* oldParent = currNode->parent();

		if (oldParent != parentNode)
		{
			// remove currNode from the children list of oldParent
			legaltreenode_v::iterator it = oldParent->children().begin();
			for ( ; it != oldParent->children().end() ; it++ )
				if ( (*it) == currNode ) break;
			assert ( it != oldParent->children().end() );
			oldParent->children().erase(it);

			// make new links
			currNode->setParent(parentNode);			
			parentNode->children().push_back(currNode);

			// make hard connection in the graph
			int var1 = currNode->variable();
			int var2 = parentNode->variable();
			m_graph->connect(var1, var2);
		}

		// now stack children which are part of the cutset
		legaltreenode_v::iterator it = currNode->children().begin();
		for ( ; it != currNode->children().end(); it++)
		{
			if ((*it)->getWCutsetHeight() >= scopeSize)
				dfsStack.push( (*it) );
		}

		parentNode = currNode;

		cout << "\n nodes in cutset " << count++;
	}

	// update nodes height
	currNode = m_tree_cutset->getRoot();
	currNode->setHeight(0);

	dfsStack.push(currNode);
	while (!dfsStack.empty())
	{
		currNode = dfsStack.top();
		dfsStack.pop();

		if(currNode != m_tree_cutset->getRoot())
		{
			currNode->setHeight( currNode->parent()->getHeight() + 1 );
		}

		// now stack children which are part of the cutset
		legaltreenode_v::iterator it = currNode->children().begin();
		for ( ; it != currNode->children().end(); it++)
			dfsStack.push( (*it) );
	}
}


// extracts flat CPTs and adds them to m_functionsConstraint
void CProblemMixed::extractDeterminism()
{	
	function_map::iterator it = m_functions.begin();
	for (; it != m_functions.end(); ++it) 
	{
		bool isStrictlyPositive = true; // to determine if the flat table has any zeros

		CFunction* fun = (*it).second;
		int size = ((CProbabilityTable*) fun)->getTableSize();
		double* tableProb = ((CProbabilityTable*) fun)->getTable();

		int* tableConstraint = new int[size];

		int i;
		for (i=0 ; i<size; i++)
		{
			if (tableProb[i] == 0)
			{
				tableConstraint[i] = 0;
				isStrictlyPositive = false;
			}
			else
				tableConstraint[i] = 1;
		}

		if (isStrictlyPositive)
		{
			delete[] tableConstraint;
			continue;
		}

		int argc = fun->getArgc();
		int* argv = new int[argc];
		assert(argv != NULL);
		for ( i=0; i<argc; i++)
			argv[i] = fun->getArgv()[i];

		int id = m_Constraints;

		CConstraintTable* constraint = new CConstraintTable(id, CONSTRAINT, argc, argv);
		constraint->setOwner(this);
		constraint->create(size, tableConstraint);
		addConstraintFunction(constraint);
	}
	 
	// Reset original functions.
	function_map& funs = m_functionsConstraint;
	function_map::iterator it_f = funs.begin();
	for (; it_f != funs.end(); ++it_f)
	{
		CFunction* fun = (*it_f).second;
		fun->setUsed(false);
	}
}

// needs parentSepSet initialized
// for each variable computes the max w, for which the var belongs to the min height w-cutset
//void CProblemMixed::wCutsetInit()
//{
//	assert(m_parentSepSet != NULL);
//
//	if (m_wCutset)
//		delete[] m_wCutset;
//
//	m_wCutset = new int[m_N];
//	memset(m_wCutset, -1, m_N * sizeof(int));
//
//	// mark the reduced cutsets
//	// here we take into account separator size (context) - not cluster size !!!
//	int var;
//	for (var = 0; var < m_N; var ++)
//	{
//		variable_v::reverse_iterator it = m_parentSepSet[var].rbegin();
//		int width = 0;
//		for ( ; it != m_parentSepSet[var].rend() ; it++)
//		{
//			m_wCutset[(*it)] = max ( width, m_wCutset[(*it)] );
//			width ++ ;
//		}
//	}	
//
//	// for every var, propagate values up to the ancestors, by max
//	for(var = 0 ; var < m_N; var++)
//	{
//		for(int descendant = 0; descendant < m_N;  descendant++)
//		{	
//			if(isDescendant(var, descendant))
//				m_wCutset[var] = max ( m_wCutset[descendant] , m_wCutset[var] );
//		}
//	}
//
//	setWCutsetHeight();
//}
//

// this is the new functions that assumes only caching at OR levels
void CProblemMixed::wCutsetInit()
{
	assert(m_parentSet != NULL);

	if (m_wCutset)
		delete[] m_wCutset;

	m_wCutset = new int[m_N];
	for(int i=0; i<m_N; i++)
		m_wCutset[i] = -1;

	// mark the reduced cutsets
	// here we take into account separator size (context) - not cluster size !!!
	int var;
	for (var = 0; var < m_N; var ++)
	{
		variable_v::reverse_iterator it = m_parentSet[var].rbegin();
		int width = 0;
		for ( ; it != m_parentSet[var].rend() ; it++)
		{
			m_wCutset[(*it)] = max ( width, m_wCutset[(*it)] );
			width ++ ;
		}
	}	

	// for every var, propagate values up to the ancestors, by max
	for(var = 0 ; var < m_N; var++)
	{
		for(int descendant = 0; descendant < m_N;  descendant++)
		{	
			if(isDescendant(var, descendant))
				m_wCutset[var] = max ( m_wCutset[descendant] , m_wCutset[var] );
		}
	}

	setWCutsetHeight();
}

void CProblemMixed::setWCutsetHeight()
{
	if (m_wCutsetHeight)
		delete[] m_wCutsetHeight;

	m_wCutsetHeight = new int[m_N];
	for(int i=0; i<m_N; i++)
		m_wCutsetHeight[i] = -1;

	CLegalTreeNode** nodes = m_tree->getNodes();

	for (int i = 0 ; i < m_N ; i++)
	{
		int height = nodes[i]->getHeight();
		int var  = nodes[i]->variable();
		int w = m_wCutset[var];

		m_wCutsetHeight[w] = max ( m_wCutsetHeight[w] , height );
	}
}	

// This function removes evidence in the network. It must be called right after
// loading the pedigree network and before preprocessing the network.
void CProblemMixed::removeEvidence(bool* barrenVariables)
{
	// Safety checks.
	assert(m_N > 0);

	// Declare variables.
	int i, d, k, var, idx, newN, newK, newID, ovar, nvar;
	map<int,int> evidence;
	map<int,int> mapOld2New;
	map<int,int> domains;
	vector<int> order;
	vector<int> cutset;
	function_v newFuns;

	// Identify evidence variables.
	idx = newN = newK = 0;
	for (var = 0; var < m_N; ++var)
	{
		if(barrenVariables!= NULL && barrenVariables[var])
				continue;

		// Check for evidence variable.
		if ((1 == m_domains[var]) || isEvidence(var))
		{
			evidence[var] = getValue(var);
			setEvidence(var, getValue(var));
		}
		else
		{			
			mapOld2New[var] = idx;
			k = getStaticDomainSize(var);
			newK = max(newK, k);
			domains[idx] = k;

			idx++;
			newN++;
		}
	}

	// Recompute the scopes of the functions.
	globalCostFromEvidence = 1.0;

	function_map::iterator itf = m_functions.begin();
	for (; itf != m_functions.end(); ++itf)
	{		
		CProbabilityTable* cpt = (CProbabilityTable*)((*itf).second);
		// check for barren variable 
		int argc = cpt->getArgc();
		int* argv = cpt->getArgv();
		if (barrenVariables != NULL && barrenVariables[argv[argc-1]])
			continue;

		bool isConstant = false;
		bool isCpt = false;
		if (m_assignment[argv[argc-1]] == -1)
			isCpt = true;

		CFunction* ncpt = NULL;
		map<int,int>::iterator it = evidence.begin();
		for (; it != evidence.end(); ++it)
		{
			int evar = (*it).first;
			int eval = (*it).second;
			
			if (cpt->isMemberOf(evar))
			{
				// Substitute evidence variable.
				cpt->substitute(evar, eval, ncpt);
				delete cpt;
				cpt = (CProbabilityTable*)ncpt;
				if(isCpt)
					cpt->setType(CPT);
				else
					cpt->setType(unknown);

				// Check for constant.
				if (cpt->isConstant())
				{
					double cval;
					cpt->getCurrentValue(&cval);
					globalCostFromEvidence *= cval;
					
					isConstant = true;
					delete cpt;

					break;
				}
			}
		}

		if (!isConstant)
			newFuns.push_back(cpt);
	}

	// Reset the problem.
	reset();

	// Build the updated repository.
	newID = 0;
	m_functions.clear();
	function_v::iterator itnf = newFuns.begin();
	for (; itnf != newFuns.end(); ++itnf)
	{
		CFunction* f = (*itnf);
		f->setOriginal(true);
		f->setID(newID);
		int* argv = f->getArgv();
		int  argc = f->getArgc();

		// Reindex the scope.
		for (i = 0; i < argc; ++i)
		{
			ovar = argv[i];
			map<int,int>::iterator itv = mapOld2New.find(ovar);
			assert(itv != mapOld2New.end());
			nvar = (*itv).second;

			argv[i] = nvar;
		}

		addFunction(f);
		newID++;
	}

	m_C = (int)m_functions.size();
	m_N = newN;
	m_K = newK;

	// Reinitialize the problem.
	init();

	// Update new domains.
	for (i = 0; i < m_N; ++i)
	{
		map<int,int>::iterator it = domains.find(i);
		if (it != domains.end())
		{
			d = (*it).second;
			m_domains[i] = d;
		}
	}

	// Check the scope of the functions.
	itf = m_functions.begin();
	for (; itf != m_functions.end(); ++itf)
	{
		CFunction* f = (*itf).second;
		int* argv = f->getArgv();
		int  argc = f->getArgc();

		// Check the scope.
		for (i = 0; i < argc; ++i)
		{
			assert(argv[i] < m_N);
		}

	}

	//// Reindex superlink order.
	//for (i = 0; i < (int)m_superlinkOrder.size(); ++i)
	//{
	//	ovar = m_superlinkOrder[i];
	//	map<int,int>::iterator itv = mapOld2New.find(ovar);
	//	assert(itv != mapOld2New.end());
	//	nvar = (*itv).second;

	//	order.push_back(nvar);
	//}

	//m_superlinkOrder.clear();
	//copy(order.begin(), order.end(), back_inserter(m_superlinkOrder));

	//// Reindex superlink cutset.
	//for (i = 0; i < (int)m_superlinkCutset.size(); ++i)
	//{
	//	ovar = m_superlinkCutset[i];
	//	map<int,int>::iterator itv = mapOld2New.find(ovar);
	//	assert(itv != mapOld2New.end());
	//	nvar = (*itv).second;

	//	cutset.push_back(nvar);
	//}

	//m_superlinkCutset.clear();
	//copy(cutset.begin(), cutset.end(), back_inserter(m_superlinkCutset));

	// Create new moral graph.
//	m_graph->init(m_functions);
//	m_connected = m_graph->isConnected();

	// Free memory.
	mapOld2New.clear();
	newFuns.clear();
	domains.clear();
	order.clear();
	cutset.clear();
}


void CProblemMixed::reset()
{
	//destroy();

	if (m_graph) delete m_graph;
	if (m_tree)	delete m_tree;
	if (m_descendants)
	{
		for (int i = 0; i < m_N; ++i)
			m_descendants[i].clear();
		delete[] m_descendants;
	}

	// Deallocate memory.
	if (m_domains) delete[] m_domains;
	if (m_ordering)	delete[] m_ordering;
	if (m_position)	delete[] m_position;
	if (m_assignment) delete[] m_assignment;
	if (m_evidence)	delete[] m_evidence;
	if (m_backupAssignment)	delete[] m_backupAssignment;
	if (m_backupEvidence) delete[] m_backupEvidence;

	// Clear the adjacency matrix
	if (m_adjacencies)
	{
		for (int i = 0; i < m_N; ++i)
			delete m_adjacencies[i];
		delete[] m_adjacencies;
	}

	m_domains = NULL;
	m_ordering = NULL;
	m_position = NULL;
	m_backupAssignment = NULL;
	m_backupEvidence = NULL;
	m_evidence = NULL;
	m_assignment = NULL;
	m_adjacencies = NULL;
	m_graph = NULL;



	if (m_graph) delete m_graph;
	if (m_bayesGraph) delete m_bayesGraph;
	if (m_constraintGraph) delete m_constraintGraph;

	if (m_tree_cutset) delete m_tree_cutset;

	if (m_buckets) delete m_buckets;

	if (m_bestSolution)	delete[] m_bestSolution;

	if (m_ordering_backup) delete[] m_ordering_backup;

	if (m_position_backup) delete[] m_position_backup;

	if (m_orderingDFS) delete[] m_orderingDFS;

	if (m_positionDFS) delete[] m_positionDFS;

	if (m_ordering_cutsetDFS) delete[] m_ordering_cutsetDFS;

	if (m_position_cutsetDFS) delete[] m_position_cutsetDFS;


//#ifdef CACHE_AT_OR_NODES
	if (m_parentSet)
	{
		for (int i = 0; i < m_N; ++i)
			m_parentSet[i].clear();
		delete[] m_parentSet;
	}

	if (m_cacheFull) delete[] m_cacheFull;

	if (m_cacheFullCounting) delete[] m_cacheFullCounting;
//#endif CACHE_AT_OR_NODES


//#ifdef CACHE_AT_AND_NODES
	if (m_parentSepSet)
	{
		for (int i = 0; i < m_N; ++i)
			m_parentSepSet[i].clear();
		delete[] m_parentSepSet;
	}

	if (m_cacheSep)	delete[] m_cacheSep;

	if (m_cacheSepCounting) delete[] m_cacheSepCounting;
//#endif CACHE_AT_AND_NODES

	if (m_advancedCacheFlagAND)	delete[] m_advancedCacheFlagAND;

	if (m_advancedCacheFlagOR) delete[] m_advancedCacheFlagOR;

	function_map::iterator it = m_functionsConstraint.begin();
	for ( ; it != m_functionsConstraint.end(); ++it)
		delete (*it).second;
	m_functionsConstraint.clear();

	if (m_treeDecomposition) delete m_treeDecomposition;	

	if (m_wCutset) delete[] m_wCutset;

	if (m_wCutsetHeight) delete[] m_wCutsetHeight;
}

  
void CProblemMixed::createShannonTrees()
{
	function_map::iterator it = m_functionsConstraint.begin();
	for ( ; it != m_functionsConstraint.end(); ++it) 
	{
		CFunction* fun = (*it).second;

		// build Shannon Tree (for constraint propagation)
		((CConstraintTable*)fun)->createShannonTree();
	}

}

// mark the barren Variables
void CProblemMixed::markBarrenVariables(bool*& barrenVariables)
{
	CGraphHash graphHash;
	// initially, mark all as barren. then propagate evidence up
	for(int i=0; i<m_N; i++)
		barrenVariables[i] = true; 


	// Create the graph hash. links will point bottom up, from leaves to root
	function_map::iterator itf = m_functions.begin();
	for (; itf != m_functions.end(); ++itf)
	{
		CFunction* f = (*itf).second;
		int argc = f->getArgc();
		int* argv = f->getArgv();

		graphHash.addNode(argv[argc-1]);
		for(int i=0; i<argc-1; i++)
		{
			// add nodes
			graphHash.addNode(argv[i]);

			// add directed edge
			graphHash.addDirectedEdge(argv[argc-1], argv[i]);
			// cout << "\n edge: [" << argv[argc-1] << " , " << argv[i] << " ] ";
		}
	}

	// now make a dfs starting from evidence nodes, and mark the unbarren variables
	stack<int> dfsStack;
//	cout << "\n initializing barren variables";
	for (int var=0; var<m_N; var++)
	{
		if ((1 == m_domains[var]) || isEvidence(var))
		{
			dfsStack.push(var);
//			cout << "\n is not barren:" << var;
		}
	}

	while (!dfsStack.empty())
	{
		int var = dfsStack.top();
		dfsStack.pop();

		barrenVariables[var] = false;

		variable_s s(*(graphHash.neighbors(var)));
		variable_s::iterator it = s.begin();
		for (; it != s.end(); ++it)
		{
			int v = (*it);
			if (barrenVariables[v])
				dfsStack.push(v);
		}
	}
}

int CProblemMixed::findVariableToComputePercentageDone(variable_v& vars)
{
	int var0 = m_ordering[0];
	int maxDepth = 100;
	vars.clear();

	int currVar = var0;
	for(int i=0; i<maxDepth; i++)
	{
		vars.push_back(currVar);
		variable_v children;
		m_tree->getChildren(currVar, children);

		if(children.size() != 1)
			break;
		
		currVar = children[0];		
	}
	
	return vars.back();
}

double CProblemMixed::updatePercentageDone()
{
	int size = topChainVariables.size();
	double total = 1.0;
	for (int i=0; i<size; i++)
		total *= (double) getStaticDomainSize(topChainVariables[i]);
	
	int k = size - 1;
	double addr = (double) (getValue(topChainVariables[k]));

	int j = getStaticDomainSize(topChainVariables[k]);
	for (--k ; k >= 0 ; k--)
	{
		// here we should add to 'addr' the product of :
		//		#-of-values-of-variable-k * (current-value-of-k + 1)
		int value = getValue(topChainVariables[k]);
		assert(value >= 0);

		addr += (double) j*value ;
		j *= getStaticDomainSize(topChainVariables[k]) ;
	}

	return ( addr/(total-1) ) * 100.0;
} 

// mark invalid values for unassigned variables and check support for var
bool CProblemMixed::supportIn(CFunction* fun, int var, CAONode* node)
{
	assert(m_assignment[var] < 0);
	assert(node->type() == AND);
	assert(m_validValuesCount[var] > 0);

	int argc = fun->getArgc();
	int* argv = fun->getArgv();

	// check each valid value of var
	for ( int val = 0; val < m_domains[var] ; val++ )
	{
		if ( !m_validValues[var][val] ) continue;

		setValue(var, val);

		int* isAssigned = new int[argc]; // array of 0 and 1 to indicate if the corresponding variable is in the current path
		bool isFullyAssigned = true; // flag to indicate if the constraint has its scope fully assigned;

		for (int j=0; j < argc; j++)
		{		
			if (m_assignment[argv[j]] < 0)
			{
				isAssigned[j] = 0;
				isFullyAssigned = false;
			}
			else
				isAssigned[j] = 1;
		}

		if (isFullyAssigned)
		{
			int value;
			((CConstraintTable*)fun)->getCurrentValue(&value);
			if( value == 0 )
			{	
				m_validValues[var][val] = false;
				--m_validValuesCount[var];
				node->addInvalidValue(var,val);
				if (m_validValuesCount[var] == 0)
				{
					setValue(var,-1);
					if (isAssigned) delete[] isAssigned;					
					return false;
				}			
			}

			// go to next value
			continue;
		}

		// find the smallest variable in the scope (rightmost one) not assigned
		int smallestVarNotYetAssigned = -1;
		for (int j = argc-1; j >= 0 ; j--)
			if (isAssigned[j] == 1)
				continue;
			else
			{
				smallestVarNotYetAssigned = j;
				break;
			}
		assert(smallestVarNotYetAssigned != -1);

		int* parents = new int[argc]; // this is the scope of the constraint
		for (int j = 0; j < argc ; j++) 
		{
			if(isAssigned[j] == 1)	
				parents[j] = m_assignment[argv[j]]; // assigned values are maintained
			else 
				parents[j] = 0;	// enumerate only through those that are not assigned;
		}
		parents[smallestVarNotYetAssigned] = -1; // to initiate the enumeration properly

		// compute how many assignments we have to enumerate; just valid values
		int nCount = 1;
		for (int j = 0; j < argc; j++)
		{
			if (isAssigned[j] == 0)
				nCount *= m_validValuesCount[argv[j]];
		}
		
		bool isConsistent = false;
		// Begin enumeration
		for (int p = 0; p < nCount; p++)
		{
			// Find next combination from the current.
			for (int pcai = argc - 1 ; pcai >= 0 ; --pcai) 
			{
				if (isAssigned[pcai] == 1)
					continue;

				++parents[pcai];
				while( !m_validValues[argv[pcai]][parents[pcai]] && parents[pcai] < getStaticDomainSize(argv[pcai]) )
				{
					++parents[pcai];
				}

				if (parents[pcai] < getStaticDomainSize(argv[pcai])) break ;

				parents[pcai] = 0 ;
				while( !m_validValues[argv[pcai]][parents[pcai]] )
					++parents[pcai];
			}

			// now check the value of the constraint for the current assignment
			int k = argc - 1;
			int addr = parents[k];
			assert(addr >= 0);
			
			int j = getStaticDomainSize(argv[k]);
			for (--k ; k >= 0 ; k--)
			{
				// here we should add to 'adr' the product of :
				//		#-of-values-of-variable-k * current-value-of-k
				int value = parents[k];
				assert(value >= 0);

				addr += j*value ;
				j *= getStaticDomainSize(argv[k]) ;
			}

			if( ((CConstraintTable*)fun)->getValueAt(addr) == 1)
			{
				isConsistent = true;
				break;
			}
		}

		if (!isConsistent)
		{
			m_validValues[var][val] = false;
			--m_validValuesCount[var];
			node->addInvalidValue(var,val);
			if (m_validValuesCount[var] == 0)
			{
				setValue(var,-1);
				if (isAssigned) delete[] isAssigned;
				//if (parents) delete[] parents;
				return false;
			}
		}

		if (isAssigned) delete[] isAssigned;
		if (parents) delete[] parents;
	}

	setValue(var,-1);
	return true;
}

bool CProblemMixed::findDeterministicVariables()
{
	m_isDeterministic = new bool[m_N];

	for (int i=0 ; i<m_N ; i++)
		m_isDeterministic[i] = false;

	bool existsDeterminism = false;
	function_map::iterator itf = m_functions.begin();
	for ( ; itf != m_functions.end(); ++itf)
	{
		CFunction* f = (*itf).second;
		if (f->type() != CPT)
			continue;

		int argc = f->getArgc();
		int* argv = f->getArgv();
		int var = argv[argc-1];

		int parents = 1;
		for (int i=0; i<argc-1; i++)
			parents *= m_domains[argv[i]];

		int addr = 0;
		for (int i=0; i<parents ; i++)
		{
			int count = 0;
			for (int j=0; j<m_domains[var]; j++)
			{
				if ( ((CProbabilityTable*) f )->getValueAt(addr) > 0 )
					count++;

				addr++;
			}
			if (count > 1) 
				break;
			else
				if ( i == parents -1 )
				{
					m_isDeterministic[var] = true;
					existsDeterminism = true;
				}
		}
	}

	return existsDeterminism;
}
void CProblemMixed::fixOrderingDeterministic(variable_s* parents)
{
	// create watch lists for non deterministic variables
	variable_s* watchList = new variable_s[m_N];
	variable_s deterministicRoot;

	for (int i=0; i<m_N; i++)
	{
		if (m_isDeterministic[i])
		{
			// find last parent in ordering
			int lastParent = -1;
			int lastPosition = -1;
			variable_s::iterator it = parents[i].begin();
			for ( ; it != parents[i].end(); ++it)
			{
				int p = (*it);
				if (m_position[p] > lastPosition)
				{
					lastParent = p;
					lastPosition = m_position[p];
				}
			}
			if(lastParent == -1)
				deterministicRoot.insert(i);  // this is a deterministic root created from eliminating evidence
			else
				watchList[lastParent].insert(i);
		}
	}

	//for (int i=0; i<m_N; i++)
	//{
	//	cout << "  _w[" << i << "]: ";
	//	variable_s::iterator it = watchList[i].begin();
	//	for ( ; it != watchList[i].end(); it++)
	//	{
	//		cout << (*it) << " ";
	//	}
	//}

	stack<int> Stack;
	for(int i=m_N-1; i>=0; i--)
	//for(int i=0; i<m_N; i++)
	{
		int v = m_ordering[i];
		if (!m_isDeterministic[v])
			Stack.push(v);
	}

	bool* used = new bool[m_N];
	for (int i=0; i<m_N; i++)
		used[i] = false;

	// first place the deterministic roots
	variable_s::iterator itroot = deterministicRoot.begin();
	for ( ; itroot != deterministicRoot.end(); itroot++ )
	{
		int v = (*itroot);
		Stack.push(v);
	}
	deterministicRoot.clear();

	// here we maintain the relative ordering of the non-deterministic variable, instantiate them in turn, and
	// follow each one by as many deterministic as possible. 
	// so this version never instantiates a deterministic variable (when its parents are not assigned)
	int count = 0;
	while(!Stack.empty())
	{
		int v = Stack.top();
		Stack.pop();

		m_ordering[count] = v;
		m_position[v] = count;
		used[v] = true;
		count++;

		variable_s::iterator it = watchList[v].begin();
		for ( ; it != watchList[v].end(); it++)
		{
			int detV = (*it);
			assert(!used[detV]);
			if (used[detV])
				cout << "uu ";
			
			// check if parents are all used (i.e. already in ordering)
			bool parentsInOrder = true;
			variable_s::iterator itp = parents[detV].begin();
			for ( ; itp != parents[detV].end(); itp++ )
			{
				int p = (*itp);
				if (!used[p])
				{
					watchList[p].insert(detV);
					parentsInOrder = false;
					break;
				}
			}

			if (parentsInOrder)
				Stack.push(detV);
		}
	}

	delete[] used;
	for (int i=0; i<m_N; i++)
		watchList[i].clear();
	delete[] watchList;
}


// in this function we maintain the initial ordering (e.g. min-fill), and after each instantiation we bring up as many
// deterministic variables as possible. This is different from previous function, which only chose to instantiate 
// non-deterministic variables when nothing else was avaialable. Here, we allow to instantiate deterministic variables
// even when their parents are not instantiated yet (to maintain the work of min-fill).
void CProblemMixed::fixOrderingDeterministic2(variable_s* parents, int* substitution)
{
	// create watch lists for non deterministic variables
	variable_s* watchList = new variable_s[m_N];
	variable_s deterministicRoot;

	int* ordering = new int[m_N];	

	for (int i=0; i<m_N; i++)
	{
		if (m_isDeterministic[i])
		{
			// find last parent in ordering
			int lastParent = -1;
			int lastPosition = -1;
			variable_s::iterator it = parents[i].begin();
			for ( ; it != parents[i].end(); ++it)
			{
				int p = (*it);
				if (m_position[p] > lastPosition)
				{
					lastParent = p;
					lastPosition = m_position[p];
				}
			}
			if(lastParent == -1)
				deterministicRoot.insert(i);  // this is a deterministic root created from eliminating evidence
			else
				watchList[lastParent].insert(i);
		}
	}

	//for (int i=0; i<m_N; i++)
	//{
	//	cout << "  _w[" << i << "]: ";
	//	variable_s::iterator it = watchList[i].begin();
	//	for ( ; it != watchList[i].end(); it++)
	//	{
	//		cout << (*it) << " ";
	//	}
	//}

	stack<int> Stack;

	//for(int i=m_N-1; i>=0; i--)
	////for(int i=0; i<m_N; i++)
	//{
	//	int v = m_ordering[i];
	//	if (!m_isDeterministic[v])
	//		Stack.push(v);
	//}

	bool* used = new bool[m_N];
	for (int i=0; i<m_N; i++)
		used[i] = false;

	// first place the deterministic roots
	variable_s::iterator itroot = deterministicRoot.begin();
	for ( ; itroot != deterministicRoot.end(); itroot++ )
	{
		int v = (*itroot);
		Stack.push(v);
	}
	deterministicRoot.clear();

	// here we maintain the relative ordering of the non-deterministic variable, instantiate them in turn, and
	// follow each one by as many deterministic as possible. 
	// so this version never instantiates a deterministic variable (when its parents are not assigned)
	int count = 0;
	int subst = -1;
	int positionInOldOrdering = 0;
	while(true)
	{
		while(!Stack.empty())
		{
			int v = Stack.top();
			Stack.pop();

			ordering[count] = v;
			used[v] = true;
			count++;
			substitution[v] = subst;

			variable_s::iterator it = watchList[v].begin();
			for ( ; it != watchList[v].end(); it++)
			{
				int detV = (*it);	
				if (used[detV])
					continue;
				
				// check if parents are all used (i.e. already in ordering)
				bool parentsInOrder = true;
				variable_s::iterator itp = parents[detV].begin();
				for ( ; itp != parents[detV].end(); itp++ )
				{
					int p = (*itp);
					if (!used[p])
					{
						watchList[p].insert(detV);
						parentsInOrder = false;
						break;
					}
				}

				if (parentsInOrder)
					Stack.push(detV);
			}
		}

		while(used[ m_ordering[positionInOldOrdering] ])
		{
			positionInOldOrdering++;
			if (positionInOldOrdering == m_N)
				break;
		}

		if( positionInOldOrdering < m_N )
		{
			Stack.push( m_ordering[positionInOldOrdering] );
			subst = m_ordering[positionInOldOrdering];
		}
		else
			break;
	}

	for (int i=0; i<m_N; i++)
	{
		int v = ordering[i];
		m_ordering[i] = v;
		m_position[v] = i;
	}

	delete[] ordering;
	delete[] used;
	for (int i=0; i<m_N; i++)
		watchList[i].clear();
	delete[] watchList;
}

void CProblemMixed::findTopologicalOrdering(variable_v& order)
{
	CGraphHash graphHash;
	stack<int> Stack;
	variable_v orderTemp;

	int* indegree = new int[m_N];
	for (int i=0; i<m_N; i++)
		indegree[i] = 0;

	// Create the graph hash. links will point top down, from root to leaves
	function_map::iterator itf = m_functions.begin();
	for (; itf != m_functions.end(); ++itf)
	{
		CFunction* f = (*itf).second;
		int argc = f->getArgc();
		int* argv = f->getArgv();		

		graphHash.addNode(argv[argc-1]);
		for(int i=0; i<argc-1; i++)
		{
			// add nodes
			graphHash.addNode(argv[i]);

			// add directed edge
			if(f->type()==CPT && !graphHash.containsEdge(argv[i], argv[argc-1]) && !graphHash.containsEdge(argv[argc-1],argv[i]))
			{
				graphHash.addDirectedEdge(argv[i], argv[argc-1]);
				indegree[argv[argc-1]]++;
			}
		}		
	}

	bool* preferred = new bool[m_N];
	for (int i=0; i<m_N; i++)
		preferred[i] = false;

	for (int i=0; i<m_N; i++)
		if (indegree[i] == 0)
			Stack.push(i);

	while (!Stack.empty())
	{
		int var = Stack.top();
		Stack.pop();
		orderTemp.push_back(var);

		variable_s s(*(graphHash.neighbors(var)));
		variable_s::iterator it = s.begin();
		// first put non preferred on the stack
		for (; it != s.end(); ++it)
		{
			int v = (*it);
			if (indegree[v] == 1  &&  !preferred[v])
				Stack.push(v);
		}
		// then put preferred on the stack
		it = s.begin();
		for (; it != s.end(); ++it)
		{
			int v = (*it);
			if (indegree[v] == 1  &&  preferred[v])
				Stack.push(v);
		}
		// then mark all neighbors as preferred and reduce indegree
		it = s.begin();
		for (; it != s.end(); ++it)
		{
			int v = (*it);
			indegree[v]--;
			preferred[v] = true;
		}

		// eliminate var from the graph
		graphHash.removeNode(var);
	}

	// need to feed back the reversed order;
	for (int i=m_N-1; i>=0; i--)
		order.push_back(orderTemp[i]);

	orderTemp.clear();
	delete[] indegree;
	delete[] preferred;
}

void CProblemMixed::computeFractionDone(int currentVar)
{
	double total;
	//for (int i=0; i<m_N; i++)
	//	logTotal += log ((double) m_domains[i]);

	//int k = size - 1;
	//double addr = (double) (getValue(topChainVariables[k]));


	//int j = getStaticDomainSize(currentVar);
	//for (--k ; k >= 0 ; k--)
	//{
	//	// here we should add to 'addr' the product of :
	//	//		#-of-values-of-variable-k * (current-value-of-k + 1)
	//	int value = getValue(currentVar);
	//	assert(value >= 0);

	//	addr += (double) j*value ;
	//	j *= getStaticDomainSize(topChainVariables[k]) ;
	//}
}
void CProblemMixed::createGridProblem(int nSide)
{
	init();

	int fCount = 0;

	for (int i=0 ; i < nSide ; i++)
		for (int j=0; j < nSide ; j++)
		{
			// (i,j) are the coordinates in grid, i is row, j is column
			// compute current variable 
			
			int var = (i * nSide) + j;

			// see if we can have a constraint to the right
			if (j < nSide - 1)
			{
				// Create a function and add it to the network.
				int argc = 2;
				int *argv = new int[argc];
				argv[0] = var;
				argv[1] = var+1;

				int tab_size = 4;
				double *tab = new double[tab_size];
				tab[0] = 1.0;
				tab[1] = 1.0;
				tab[2] = 1.0;
				tab[3] = 0.0;
				
				CProbabilityTable *fn = new CProbabilityTable(fCount, CPT, argc, argv, tab_size, tab);
				fn->setOwner(this);
				addFunction(fn);
				fCount++;
			}

			// see if we can have a constraint below
			if (i < nSide - 1)
			{
				// Create a function and add it to the network.
				int argc = 2;
				int *argv = new int[argc];
				argv[0] = var;
				argv[1] = var+nSide;

				int tab_size = 4;
				double *tab = new double[tab_size];
				tab[0] = 1.0;
				tab[1] = 1.0;
				tab[2] = 1.0;
				tab[3] = 0.0;
				
				CProbabilityTable *fn = new CProbabilityTable(fCount, CPT, argc, argv, tab_size, tab);
				fn->setOwner(this);
				addFunction(fn);
				fCount++;			
			}
		}
		
	m_r = fCount;
}