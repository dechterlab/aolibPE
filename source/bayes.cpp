//// bayes.cpp -- Bayesian Networks.
//
///*
// * Copyright (C) 2004 Radu Marinescu
// *
// * This program is free software; you can redistribute it and/or modify
// * it under the terms of the GNU General Public License as published by
// * the Free Software Foundation; either version 2, or (at your option)
// * any later version.
// *
// * This program is distributed in the hope that it will be useful,
// * but WITHOUT ANY WARRANTY; without even the implied warranty of
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// * GNU General Public License for more details.
// *
// * You should have received a copy of the GNU General Public License
// * along with this program; if not, write to the Free Software
// * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// */
//
///*
// * NOTE: This is an interal header file.
// * You should not attempt to use it directly.
// */
//
//
//#include "bayes.h"
//#include "ConstraintTable.h"
//#include "cpt.h"
//#include <memory.h>
//
///////////////////////////////////////////////////////////////////////////
//// CProblemBayes class implementation.
//CProblemBayes::CProblemBayes() : CProblem()
//{
//	m_type = UNKNOWN;
//	m_connected = false;
//
//	m_graph = NULL;
//	m_tree = NULL;
//	m_buckets = NULL;
//
//	m_paramParents		= UNKNOWN;
//	m_paramConnectivity = UNKNOWN;
//	m_paramHeight		= UNKNOWN;
//	m_paramWidth		= UNKNOWN;
//	m_paramSigma		= UNKNOWN;
//
//	m_bestSolutionCost = 0;
//	m_bestSolution = NULL;
//
//	m_outputFile = NULL;
//
//	m_levelInfo = NULL;
//
//	m_orderingDFS = NULL;
//	m_positionDFS = NULL;
//
//	m_descendants = NULL;
//}
//
//CProblemBayes::CProblemBayes(int t, int n, int k, int p) : CProblem(n, k, p), m_type(t)
//{
//	m_connected = false;
//
//	m_graph = NULL;
//	m_tree = NULL;
//	m_buckets = NULL;
//
//	m_paramParents		= UNKNOWN;
//	m_paramConnectivity = UNKNOWN;
//	m_paramHeight		= UNKNOWN;
//	m_paramWidth		= UNKNOWN;
//	m_paramSigma		= UNKNOWN;
//
//	m_bestSolutionCost = 0;
//	m_bestSolution = NULL;
//
//	m_outputFile = NULL;
//
//	m_levelInfo = NULL;
//
//	m_orderingDFS = NULL;
//	m_positionDFS = NULL;
//
//	m_descendants = NULL;
//}
//
//// Destroy the problem instance.
//CProblemBayes::~CProblemBayes()
//{
//	if (m_graph)
//		delete m_graph;
//
//	if (m_tree)
//		delete m_tree;
//
//	if (m_buckets)
//		delete m_buckets;
//
//	if (m_bestSolution)
//		delete[] m_bestSolution;
//
//	if (m_orderingDFS)
//		delete[] m_orderingDFS;
//
//	if (m_positionDFS)
//		delete[] m_positionDFS;
//
//	if (m_descendants)
//	{
//		for (int i = 0; i < m_N; ++i)
//			m_descendants[i].clear();
//		delete[] m_descendants;
//	}
//}
//
//int CProblemBayes::getTreeHeight()
//{
//	if (m_tree)
//		return m_tree->height();
//
//	return 0;
//}
//
//variable_v& CProblemBayes::getTreeDescendants(int var)
//{
//	return m_descendants[var];
//}
//
//int CProblemBayes::getInducedWidth()
//{
//	if (m_graph)
//		return m_graph->getWidth();
//
//	return 0;
//}
//
//int* CProblemBayes::getOrdering(bool dfsOrd)
//{
//	return ((!dfsOrd) ? m_ordering : m_orderingDFS);
//}
//
//int* CProblemBayes::getPosition(bool dfsOrd)
//{
//	return ((!dfsOrd) ? m_position : m_positionDFS);
//}
//
//function_v& CProblemBayes::getBucket(int var, bool dfsOrd)
//{
//	// Safety checks.
//	assert(m_buckets != NULL);
//
//	int pos = (dfsOrd) ? m_positionDFS[var] : m_position[var];
//	CBucket* bkt = m_buckets->getBucketAt(pos);
//	assert(bkt != NULL);
//
//	return bkt->functions();
//}
//
//bool CProblemBayes::init()
//{
//	// Safety checks.
//	assert(m_evidence == NULL);
//	assert(m_assignment == NULL);
//	assert(m_adjacencies == NULL);
//	assert(m_graph == NULL);
//
//	m_assignment = new int[m_N];
//	for(int i=0; i<m_N; i++)
//		m_assignment[i] = -1;
//
//	m_evidence = new bool[m_N];
//	for(int i=0; i<m_N; i++)
//		m_evidence[i] = false;
//
//	m_backupAssignment = new int[m_N];
//	for(int i=0; i<m_N; i++)
//		m_backupAssignment[i] = -1;
//
//	m_backupEvidence = new bool[m_N];
//	for(int i=0; i<m_N; i++)
//		m_backupEvidence[i] = false;
//
//	m_domains = new int[m_N];
//	for (int k = 0; k < m_N; ++k) m_domains[k] = m_K;
//
//	// Initialize the moral graph.
//	m_graph = new CGraph(m_N);
//
//	return true;
//}
//
//bool CProblemBayes::initialize(function_map& functions_bayes, function_map& functions_constraint, CProblemMixed* prob, int algType)				// Initialize an auxiliary network from a mixed one
//{
//	init();
//
//	// copy the bayes functions
//	function_map::iterator it = functions_bayes.begin();
//	for (; it != functions_bayes.end(); ++it)
//	{		
//		CFunction* fun = (*it).second;
//		int argc = fun->getArgc();
//		int* argv = new int[argc];
//		memcpy(argv, fun->getArgv(), argc * sizeof(int));
//
//		CProbabilityTable* cpt = new CProbabilityTable(fun->getID(), CPT, argc, argv);
//		cpt->setOwner(this);
//		cpt->init(); //now m_tableSize and m_table are initialized
//		memcpy(cpt->getTable(), ((CProbabilityTable*)fun)->getTable(), ((CProbabilityTable*)fun)->getTableSize() * sizeof(double)); //copy the cpt
//
//		cpt->verify();
//		addFunction(cpt);
//	}
//
//	// create the new auxiliary functions from the constraints
//	function_map::iterator iter = functions_constraint.begin();
//	int index = 0; // to define the new aux vars 
//	for (; iter != functions_constraint.end(); ++iter)
//	{
//		CFunction* constraint = (*iter).second;
//		int argc = 1 + constraint->getArgc(); // add the new aux var in the count
//		int* argv = new int[argc];
//		memcpy(argv, constraint->getArgv(), (argc - 1)  * sizeof(int));
//		// robert
//		int new_var = prob->getN() + index;
//		argv[argc-1] = new_var;
//		m_domains[new_var] = 2;
//		setEvidence(new_var, 1);
//		index++;
//
//		CProbabilityTable* cpt = new CProbabilityTable(new_var, CPT, argc, argv);
//		cpt->setOwner(this);
//		cpt->init();
//
//		for(int i=0; i < ((CConstraintTable*)constraint)->getTableSize(); i++)
//		{
//			if( ((CConstraintTable*)constraint)->getValueAt(i) == 1)
//			{
//				cpt->setValueAt(2*i, (double) 0.0);
//				cpt->setValueAt(2*i + 1, (double) 1.0);
//			}
//			else
//			{
//				cpt->setValueAt(2*i, (double) 1.0);
//				cpt->setValueAt(2*i + 1, (double) 0.0);
//			}
//		}
//		cpt->verify();
//		addFunction(cpt);
//	}
//
//
//
//	// Create the moral graph.
//	m_graph->init(m_functions);
//	m_connected = m_graph->isConnected();
////	return true;
//
//
//	//////////////////////
//	// Safety checks.
//	assert(m_graph != NULL);
//
//	if (!m_connected)
//		return false;	// not connected.
//
//	// Copy variable ordering from prob, put aux vars first
//	
//	// Allocate memory.
//	int* ordering = new int[m_N];
//	int* position = new int[m_N];
//	assert(ordering != NULL);
//	assert(position != NULL);
//
//	// if we do AO search, put the aux vars first in the ordering (to be instantiated first by and/or search)
//	if(algType == AO_AUX)
//	{
//		int Nprob = prob->getN();
//		int Ncons = m_N - Nprob;
//
//		int i;
//		for (i=0; i < Ncons; i++)
//			ordering[i] =  Nprob + i;
//		for (i = 0; i < Nprob ; i++)
//			ordering[Ncons + i] = prob->getOrdering(false)[i];
//
//		for (i = 0; i < Nprob; i++)
//			position[i] = prob->getPosition(false)[i] + Ncons;
//		for (i = 0; i < Ncons; i++)
//			position[Nprob + i] = i;
//
//		m_ordering = ordering;
//		m_position = position;
//
//		// add a clique of connections between auxiliary vars (this helps to create good legal tree
//		// with auxiliary vars at the beginning
//		for (i = 0; i < Ncons - 2 ; i++)
//			for (int j = i+1; j < Ncons ; j++)
//				m_graph->connect(Nprob + i, Nprob + j);
//	}
//
//	// if we do BE, then put the auxiliary vars last in the ordering, to be in the first buckets to be processed
//	if(algType == BE_AUX)
//	{
//		int Nprob = prob->getN();
//		int Ncons = m_N - Nprob;
//
//		int i;
//		for (i=0; i<Nprob; i++)
//		{
//			ordering[i] = prob->getOrdering(false)[i];
//			position[i] = prob->getPosition(false)[i];
//		}
//
//		for (i=0; i<Ncons; i++)
//		{
//			ordering[Nprob + i] = Nprob + i;
//			position[Nprob + i] = Nprob + i;
//		}
//
//		m_ordering = ordering;
//		m_position = position;
//	}
//
//	
////	printf("\n");
////	for (i=0; i<m_N; i++)
////		printf("%d ", m_ordering[i]);
////	printf("\n");
////	for (i=0; i<m_N; i++)
////		printf("%d ", m_position[i]);
//
//
//	m_graph->setWidth( m_graph->width(ordering, position) );
//
//	{
//		// Create the rooted tree arrangement from the induced graph.
//		assert(m_tree == NULL);
//		m_tree = new CLegalTree(m_N, m_ordering, m_position);
//
//		//Robert
//		m_tree->create(m_graph);
//		//m_tree->createChain(m_graph);
//
//		m_tree->order(m_orderingDFS, m_positionDFS);
//
//		m_descendants = new variable_v[m_N];
//		assert(m_descendants != NULL);
//		for (int var = 0; var < m_N; ++var)	
//			m_tree->descendants(var, m_descendants[var]);
//
//		memcpy(m_ordering, m_orderingDFS, m_N*sizeof(int));
//		memcpy(m_position, m_positionDFS, m_N*sizeof(int));
//	}
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
//	m_buckets->init(true);
//
//	return true;
//
//}
//
//
//// This function asserts likely evidence for the Bayesian network.
//bool CProblemBayes::assertEvidence(int ne)
//{
//	// Safety checks.
//	assert(m_N > 0);
//
//	if (0 == ne)
//		return false;
//
//	// Reset previous evidence.
//	for(int i=0; i<m_N; i++)
//	{
//		m_evidence[i] = false;
//		m_assignment[i] = -1;
//	}
//
//	// Declare variables.
//	int numVarWithValues = 0, i;
//
//	// Generate a value for every variable.
//	while (numVarWithValues < m_N) 
//	{
//		function_map::iterator it = m_functions.begin();
//		for (; it != m_functions.end(); ++it) 
//		{
//			CFunction* fun = (*it).second;
//			CProbabilityTable* pt = dynamic_cast<CProbabilityTable*>(fun);
//			assert(fun != NULL);
//
//			int child = pt->getChild();
//			if (m_assignment[child] >= 0) continue ;
//			
//			// Check if all parents have a value.
//			int n = pt->getNumberOfParents() ;
//			for (i = 0 ; i < n ; i++) 
//			{
//				int p = pt->getParent(i);
//				if (!isEvidence(p)) break;
//			}
//
//			// If all parents have a value, set a value for the child.
//			if (i < n) continue;
//
//			// Compute the total sum of the PT for this variable given fixed parent values.
//			double x = 0.0;
//			for (i = 0; i < m_domains[child]; ++i) 
//			{
//				setEvidence(child, i);
//
//				double f;
//				pt->getCurrentValue(&f);
//				
//				x += f;
//			}
//
//			double y = x * RandUniformDouble();
//			x = 0.0;
//			for (i = 0; i < m_domains[child]; ++i) 
//			{
//				setEvidence(child, i);
//
//				double f;
//				pt->getCurrentValue(&f);
//				
//				x += f ;
//				
//				if (x >= y) 
//				{
//					setEvidence(child, i);
//
//					break ;
//				}
//			}
//
//			++numVarWithValues;
//		}
//	}
//
//	// Pick randomly evidence variables.
//	variable_s evidence;
//	for (int ei = 0 ; ei < ne ; ei++) 
//	{
//		while (1) 
//		{
//			int evar = RandUniform(m_N);
//
//			// Check if already assigned.
//			if (evidence.find(evar) == evidence.end())
//			{
//				evidence.insert(evar);
//				break;
//			}
//		}
//	}
//
//	// Set evidence variables.
//	for (int evar = 0 ; evar < m_N ; ++evar) 
//	{
//		if (evidence.find(evar) == evidence.end())
//		{
//			resetValue(evar);
//			resetEvidence(evar);
//		}
//	}
//
//	// Free temporary buffers.
//	evidence.clear();
//
//	return true;
//}
//
//bool CProblemBayes::create()
//{
//	// Safety checks.
//	assert(m_ordering == NULL);
//	assert(m_position == NULL);
//
//	// Init the instance.
//	init();
//
//	switch (m_type)
//	{
//	case GRAPH_RANDOM:
//		{
//			// Create a random bayesian structure.
//			// For each child set the number of parents.
//			int i, j;
//			int* numParents = new int[m_N];
//			for (i = 0; i < m_N; ++i) 
//			{
//				double r = RandUniformDouble();
//				if (r <= .2)
//					numParents[i] = 0;
//				else if (r <= .3)
//					numParents[i] = 1;
//				else if (r <= .5)
//					numParents[i] = 2;
//				else if (r <= .75)
//					numParents[i] = 3;
//				else if (r <= .95)
//					numParents[i] = 4;
//				else
//					numParents[i] = 5;
//				
//				numParents[i] = min(numParents[i], i);
//				numParents[i] = min(numParents[i], m_paramConnectivity);
//			}
//
//			// Create graph structure.
//			for (int child = 0; child < m_N; ++child) 
//			{
//				variable_s parents;
//				for (j = 0; j < numParents[child];) 
//				{
//					int parent;
//					double r = RandUniformDouble();
//					if (child > m_paramConnectivity)
//						parent = (int)(i - 1 - floor(m_paramConnectivity * r));
//					else
//						parent = (int) floor(child * r);
//
//					if (parents.find(parent) == parents.end())
//					{
//						parents.insert(parent);
//						++j;
//					}
//				}
//
//				// Create a function with child and parents.
//				int argc = 1 + numParents[child];
//				int* argv = new int[argc];
//				assert(argv != NULL);
//		
//				// Set the function's parents.
//				int pos = 0;
//				variable_s::iterator it = parents.begin();
//				for (; it != parents.end(); ++it) argv[pos++] = (*it);
//
//				// Set the function's child (last position).
//				argv[pos] = child;
//
//				// Add the function to the problem instance.
//				int id = child;
//				CProbabilityTable* cpt = new CProbabilityTable(id, CPT, argc, argv);
//				cpt->setOwner(this);
//				cpt->create(PT_UNIFORM);
//				addFunction(cpt);
//
//				parents.clear();
//			}
//
//			// Create the moral graph.
//			m_graph->init(m_functions);
//			m_connected = m_graph->isConnected();
//
//			break;
//		}
//	case GRAPH_RANDOM_FIXED:
//		{
//			// Create a random bayesian structure 
//			// with fixed number of parents per child.
//
//			int N = m_N;
//			int P = m_paramParents;
//			int C = m_C;
//
//			int* ordering = new int[N];
//			int* position = new int[N];
//
//			// Create a random ordering of the variables.
//			int i;
//			for (i = 0; i < N; ++i)
//			{
//				ordering[i] = i;
//				position[i] = i;
//			}
//
//			// Randomly, switch pairs of variables.
//			for (i = 0 ; i < N ; ++i) 
//			{
//				int j = RandUniform(N);
//			
//				// Switch variable i and j.
//				int k = ordering[j];
//				ordering[j] = ordering[i];
//				ordering[i] = k;
//
//				position[ordering[j]] = j;
//				position[ordering[i]] = i;
//			}
//
//			bool* cpts = new bool[N];
//			for(int i=0; i<N; i++)
//				cpts[i] = false;
//
//			int count = 0;
//			while (count < C)
//			{
//				// Pick the child.
//				int child = ordering[RandUniform(N-P)];
//
//				// Check if variable was visited
//				if (cpts[child]) continue;
//				cpts[child] = 1;
//
//				// Pick P parents for the child.
//				int numHigherVars = N - position[child] - 1;
//				// Notice : number of parents cannot be larger than 'numHigherVars'.
//				int numParents = P;
//				if (numParents > numHigherVars) 
//					numParents = numHigherVars;
//				
//				variable_s parents;
//				for (i = 0; i < numParents;) 
//				{
//					int parent = ordering[position[child] + 1 + RandUniform(numHigherVars)];
//					// Check that this parent is not the same as the child.
//					if (child == parent) continue ;
//					// Check other parents.
//					if (parents.find(parent) == parents.end())
//					{
//						parents.insert(parent);
//						++i;
//					}
//				}
//
//				// Create a function with child and parents.
//				int argc = 1 + numParents;
//				int* argv = new int[argc];
//				assert(argv != NULL);
//
//				// Set the function's parents.
//				int pos = 0;
//				variable_s::iterator it = parents.begin();
//				for (; it != parents.end(); ++it) argv[pos++] = (*it);
//
//				// Set the function's child (last position).
//				argv[pos] = child;
//
//				// Add the function to the problem instance.
//				int id = child;
//				CProbabilityTable* cpt = new CProbabilityTable(id, CPT, argc, argv);
//				cpt->setOwner(this);
//				cpt->create(PT_UNIFORM);
//				addFunction(cpt);
//
//				parents.clear();
//
//				++count;
//			}
//
//			// Create priors.
//			for (i = 0; i < N; ++i)
//			{
//				if (cpts[i]) continue;
//
//				int argc = 1;
//				int* argv = new int[argc];
//				assert(argv != NULL);
//
//				// Set the function's child (last position).
//				argv[0] = i;
//					
//				// Add the function to the problem instance.
//				int id = i;
//				CProbabilityTable* cpt = new CProbabilityTable(id, CPT, argc, argv);
//				cpt->setOwner(this);
//				cpt->create(PT_UNIFORM);
//				addFunction(cpt);
//
//			}
//
//			// Create the moral graph.
//			m_graph->init(m_functions);
//			m_connected = m_graph->isConnected();
//
//			delete[] ordering;
//			delete[] position;
//			delete[] cpts;
//
//			break;
//		}
//	case GRAPH_GRID:
//		{
//			break;
//		}
//	case GRAPH_CODING:
//		{
//			break;
//		}
//	};
//
//	return true;
//}
//
//
//void CProblemBayes::createCustom()
//{
//	// Create a problem.
//	assert(m_N > 0);
//	assert(m_K > 0);
//
//	m_N = 7;
//	m_K = 2;
//
//	init();
//
//	// Create CPTs.
//	int argc;
//	int* argv = NULL;
//	CProbabilityTable* cpt = NULL;
//
//	// CPT 0.
//	argc = 1;
//	argv = new int[argc];
//	argv[0] = 0;
//
//	cpt = new CProbabilityTable(0, CPT, argc, argv);
//	cpt->setOwner(this);
//	cpt->create(PT_UNIFORM);
//	addFunction(cpt);
//
//	// CPT 1.
//	argc = 2;
//	argv = new int[argc];
//	argv[0] = 0;
//	argv[1] = 1;
//
//	cpt = new CProbabilityTable(1, CPT, argc, argv);
//	cpt->setOwner(this);
//	cpt->create(PT_UNIFORM);
//	addFunction(cpt);
//
//	// CPT 2.
//	argc = 2;
//	argv = new int[argc];
//	argv[0] = 0;
//	argv[1] = 2;
//
//	cpt = new CProbabilityTable(2, CPT, argc, argv);
//	cpt->setOwner(this);
//	cpt->create(PT_UNIFORM);
//	addFunction(cpt);
//
//	// CPT 3.
//	argc = 2;
//	argv = new int[argc];
//	argv[0] = 1;
//	argv[1] = 3;
//
//	cpt = new CProbabilityTable(3, CPT, argc, argv);
//	cpt->setOwner(this);
//	cpt->create(PT_UNIFORM);
//	addFunction(cpt);
//
//	// CPT 4.
//	argc = 3;
//	argv = new int[argc];
//	argv[0] = 0;
//	argv[1] = 1;
//	argv[2] = 4;
//
//	cpt = new CProbabilityTable(4, CPT, argc, argv);
//	cpt->setOwner(this);
//	cpt->create(PT_UNIFORM);
//	addFunction(cpt);
//
//	// CPT 5.
//	argc = 2;
//	argv = new int[argc];
//	argv[0] = 2;
//	argv[1] = 5;
//
//	cpt = new CProbabilityTable(5, CPT, argc, argv);
//	cpt->setOwner(this);
//	cpt->create(PT_UNIFORM);
//	addFunction(cpt);
//
//	// CPT 6.
//	argc = 2;
//	argv = new int[argc];
//	argv[0] = 2;
//	argv[1] = 6;
//
//	cpt = new CProbabilityTable(6, CPT, argc, argv);
//	cpt->setOwner(this);
//	cpt->create(PT_UNIFORM);
//	addFunction(cpt);
//
//	// Create the moral graph.
//	m_graph->init(m_functions);
//	m_connected = m_graph->isConnected();
//
//}
//
//// This function preprocess the network's graph. It creates the
//// variable ordering, computes the induced width and from the induced
//// graph it creates the rooted tree arrangement [Bayardo96].
//bool CProblemBayes::preprocess(int voType, bool withTree)
//{
//	// Safety checks.
//	assert(m_graph != NULL);
//
//	if (!m_connected)
//		return false;	// not connected.
//
//	// Create variable ordering, induced graph.
//	m_graph->order(voType, m_ordering, m_position);
//
//	// Check if we have to create the rooted tree.
//	if (withTree)
//	{
//		// Create the rooted tree arrangement from the induced graph.
//		assert(m_tree == NULL);
//		m_tree = new CLegalTree(m_N, m_ordering, m_position);
//
//		//Robert
//		m_tree->create(m_graph);
//		//m_tree->createChain();
//
//		m_tree->order(m_orderingDFS, m_positionDFS);
//
//		m_descendants = new variable_v[m_N];
//		assert(m_descendants != NULL);
//		for (int var = 0; var < m_N; ++var)	
//			m_tree->descendants(var, m_descendants[var]);
//
//		memcpy(m_ordering, m_orderingDFS, m_N*sizeof(int));
//		memcpy(m_position, m_positionDFS, m_N*sizeof(int));
//	}
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
//	m_buckets->init(withTree);
//
//	return true;
//}
//
//void CProblemBayes::print()
//{
//	int i;
//
//	for (i = 0; i < m_N; ++i)
//	{
//		if (isEvidence(i))
//		{
//			cout << "\n - evidence: " << i
//				<< " set to: " << getValue(i);
//		}
//	}
//
//	// Ordering:
//	cout << "\n Ordering: ";
//	for (i = 0; i < m_N; ++i)
//	{
//		cout << m_ordering[i] << " ";
//	}
//	cout << endl;
//
//	m_buckets->print();
//
//	// Legal Tree.
//	m_tree->print();
//
//	// Descendants in legal tree.
//	cout << "\n Legal Tree Descendants: ";
//	for (i = 0; i < m_N; ++i)
//	{
//		cout << endl;
//		copy(m_descendants[i].begin(), m_descendants[i].end(),
//			ostream_iterator<int>(cout, " "));
//	}
//
//	function_map::iterator it = m_functions.begin();
//	for (; it != m_functions.end(); ++it)
//	{		
//		CFunction* fun = (*it).second;
//		((CProbabilityTable*) fun)->print(true, true);
//	}
//	cout << "this was problem bayes";
//
//	cout << endl;
//}
//
//bool CProblemBayes::precompute(int ibound, double& cpuPrecompute)
//{
//	// Safety checks.
//	assert(m_buckets != NULL);
//
//	// Declare variables.
//	int k;
//	double* bounds = NULL;
//
//	// Run MBE with augmentation.
//	clock_t c_start = clock();
//	bool bOk = m_buckets->process(MBE_AUGMENT, ibound, k, bounds);
//	clock_t c_end = clock();
//
//	// Record pre-computing CPU time.
//	cpuPrecompute = (double)(c_end - c_start) / 1000.0;
//
//	delete[] bounds;
//
//	return bOk;
//}
//
//void CProblemBayes::test()
//{
//	int k;
//	double* bounds = NULL;
//
//	m_buckets->process(MBE_SIMPLE, 2, k, bounds);
//
//	cout << "\n MBE results: [";
//	for (int i = 0; i < k; ++i)
//	{
//		cout << bounds[i] << " ";
//	}
//	cout << "]" << endl;
//
//	// Free temporary buffers.
//	delete[] bounds;
//}
//
//// This function initialize the domains.
//bool CProblemBayes::initLevelInfo()
//{
//	// Destroy previous level manager.
//	destroyLevelInfo();
//
//	// Create searh level manager.
//	m_levelInfo = new SLINFO;
//	m_levelInfo->prev = NULL;
//	m_levelInfo->level = -1;
//	m_levelInfo->status = UNKNOWN;
//
//	int N = m_N;
//	int K = m_K;
//
//	// Create domain list/cost.
//	m_levelInfo->domainSize = new int[m_N];
//	m_levelInfo->domainCost = new double[m_N * m_K];
//	m_levelInfo->domainState = new bool[m_N * m_K];
//	
//	int addr;
//	for (int var = 0; var < N; ++var) 
//	{
//		m_levelInfo->domainSize[var] = getStaticDomainSize(var);
//		addr = var*K;
//		for (int j = 0; j < m_levelInfo->domainSize[var]; ++j, ++addr) 
//		{
//			m_levelInfo->domainCost[addr] = 0.0;		// Value cost.
//			if (isEvidence(var))
//			{
//				if (j == getValue(var))
//					m_levelInfo->domainState[addr] = true;
//				else
//					m_levelInfo->domainState[addr] = false;
//			}
//			else
//			{
//				m_levelInfo->domainState[addr] = true;		// Value state.
//			}
//		}
//	}
//
//	return true;
//}
//
//// This function destroys the search level information.
//void CProblemBayes::destroyLevelInfo()
//{
//	// Safety checks.
//	if (m_levelInfo)
//	{
//		while (NULL != m_levelInfo)
//		{
//			SLINFO* temp = m_levelInfo;
//			m_levelInfo = temp->prev;
//			
//			temp->destroy();
//			delete temp;
//		}
//
//		m_levelInfo = NULL;
//	}
//}
//
//// This function saves the current search level information
//// onto domain stack.
//void CProblemBayes::pushLevelInfo()
//{
//	// Safety checks.
//	assert(m_levelInfo != NULL);
//
//	// Create next level info structure.
//	SLINFO* temp = new SLINFO;
//	temp->prev = m_levelInfo;
//	temp->level = m_levelInfo->level + 1;
//	temp->status = UNKNOWN;
//	temp->domainSize = new int[m_N];
//	temp->domainCost = new double[m_N * m_K];
//    temp->domainState = new bool[m_N * m_K];
//
//	memcpy(temp->domainSize, m_levelInfo->domainSize, m_N * sizeof(int));
//	memcpy(temp->domainCost, m_levelInfo->domainCost, m_N * m_K * sizeof(double));
//	memcpy(temp->domainState, m_levelInfo->domainState, m_N * m_K * sizeof(bool));
//
//	// Update pointer link.
//	m_levelInfo = temp;
//}
//
//// This function restores the previous search level info.
//void CProblemBayes::popLevelInfo()
//{
//	// Safety checks.
//	assert(m_levelInfo != NULL);
//
//	SLINFO* temp = m_levelInfo->prev;
//	
//	// Destroy current top of the stack.
//	m_levelInfo->destroy();
//	delete m_levelInfo;
//
//	// Update pointer link.
//	m_levelInfo = temp;
//}
//
//// This function saves the current assignment as the best one.
//void CProblemBayes::saveCurrentSolution()
//{
//	// Safety checks.
//	assert(m_bestSolution != NULL);
//	assert(m_assignment != NULL);
//	assert(m_N > 0);
//
//	memcpy(m_bestSolution, m_assignment, m_N * sizeof(int));
//}
//
//// This function computes the cost of the current (partial) assignment.
//// If the variable in the argument is -1 then we have a full assignment,
//// otherwise we compute the utility of the partial assignment up to the
//// variable.
//double CProblemBayes::utility(int var)
//{
//	// Safety checks.
//	assert(m_assignment != NULL);
//	assert(m_N > 0);
//
//	double util = 1.0;
//	
//	if (-1 == var)
//	{
//		// Evaluate a full assignment.
//		function_map::iterator it = m_functions.begin();
//		for (; it != m_functions.end(); ++it)
//		{
//			double p;
//			CFunction* fun = (*it).second;
//
//			fun->getCurrentValue(&p);
//
//			util *= p;
//		}
//	}
//	else
//	{
//		// Get current search level.
//		int level = m_position[var];
//
//		// Evaluate a partial assignment.
//		for (int pos = 0; pos <= level; ++pos)
//		{
//			// Process original functions only.
//			CBucket* bucket = m_buckets->getBucketAt(pos);
//			function_v& funs = bucket->functions();
//			function_v::iterator it = funs.begin();
//			for (; it != funs.end(); ++it)
//			{
//				CFunction* fun = (*it);
//				if (!fun->isOriginal())
//					continue;	// skip non-original functions.
//
//				double p;
//				assert(-1 != m_assignment[var]);
//				fun->getCurrentValue(&p);
//
//				util *= p;
//			}
//		}
//	}
//
//	return util;
//}
//
//// This function evalues the heuristic value costs at current search level.
//void CProblemBayes::evalStaticHeuristic(int var)
//{
//	// Safety checks.
//	assert(m_ordering != NULL);
//	assert(m_position != NULL);
//	assert((var >= 0) && (var < m_N));
//
//	// Get current search tree level.
//	int level = m_levelInfo->level;
//	int k = getStaticDomainSize(var);
//
//	double* domainCost = m_levelInfo->domainCost;
//	assert(domainCost != 0);
//
//	// Current variable's address index.
//	int baseAddr = var * m_K;
//	
//	// Heuristic's f,g,h components.
//	double* g = new double[k];
//	double* h = new double[k];
//	double* f = new double[k];
//	for(int i=0; i<k; i++)
//	{
//		g[i] = 0;
//		f[i] = 0;
//		h[i] = 0;
//	}
//
//	// Get the current level's bucket.
//	CBucket* bucket = m_buckets->getBucketAt(level);
//
//	// Compute "h" component.
//	bucket->evalStaticHeuristic(k, h);
//
//	// Check for evidence and keep the value.
//	if (isEvidence(var))
//	{
//		// Get evidence value.
//		int i = getValue(var);
//
//		// Update the address.
//		int addr = baseAddr + i;
//	
//		// Compute "g" component.
//		g[i] = utility(var);
//
//		// Heuristic evaluation function.
//		f[i] = g[i] * h[i];
//
//		// Update current value's cost.
//		domainCost[addr] = f[i];
//	}
//	else
//	{
//		// Backup current assignment.
//		backupAssignment();
//
//		// Compute "g" component.
//		for (int i = 0; i < k; ++i)
//		{
//			// Update the address.
//			int addr = baseAddr + i;
//
//			setValue(var, i);
//			g[i] = utility(var);
//
//			// Heuristic evaluation function.
//			f[i] = g[i] * h[i];
//
//			// Update current value's cost.
//			domainCost[addr] = f[i];
//		}
//
//		// Restore current assignment.
//		restoreAssignment();
//	}
//
//	delete[] g;
//	delete[] h;
//	delete[] f;
//}
//
//void CProblemBayes::evalDynamicHeuristic(int var, int ibound)
//{
//	// Safety checks.
//	assert(m_buckets != NULL);
//
//	// Declare variables.
//	int k;
//	double* h = NULL;
//
//	// Get current variable's level (position).
//	int level = m_position[var];
//
//	// Run MBE with current variable as first in the ordering.
//	m_buckets->process(MBE_PARTIAL, ibound, k, h, level);
//
//	// Update current value costs.
//	double* domainCost = m_levelInfo->domainCost;
//	assert(domainCost != 0);
//
//	// Current variable's address index.
//	int baseAddr = var * m_K;
//	double* g = new double[k];
//	double* f = new double[k];
//	for(int i=0; i<k; i++)
//	{
//		f[i] = 0;
//		g[i] = 0; 
//	}
//	
//	// Check for evidence and keep the value.
//	if (isEvidence(var))
//	{
//		// Get evidence value.
//		int i = getValue(var);
//
//		// Update the address.
//		int addr = baseAddr + i;
//	
//		// Compute "g" component. Need to refer to the previous variable.
//		int prev = (level > 0) ? m_ordering[level - 1] : var;
//		g[i] = (level > 0) ? utility(prev) : 1.0;
//
//		// Heuristic evaluation function.
//		f[i] = g[i] * h[i];
//
//		// Update current value's cost.
//		domainCost[addr] = f[i];
//	}
//	else
//	{
//		// Backup current assignment.
//		backupAssignment();
//
//		for (int i = 0; i < k; ++i)
//		{
//			// Update the address.
//			int addr = baseAddr + i;
//
//			// Compute "g" component. Need to refer to previous variable.
//			int prev = (level > 0) ? m_ordering[level - 1] : var;
//			g[i] = (level > 0) ? utility(prev) : 1.0;
//
//			// Heuristic evaluation function.
//			f[i] = g[i] * h[i];
//
//			// Update current value's cost.
//			domainCost[addr] = f[i];
//		}
//
//		// Restore current assignment.
//		restoreAssignment();
//	}
//
//	delete[] h;
//	delete[] g;
//	delete[] f;
//}
//
//// This function selects the next variable for a static ordering.
//int CProblemBayes::selectVariable(int svType)
//{
//	// Safety checks.
//	assert(m_assignment != NULL);
//	assert(m_ordering != NULL);
//	assert(m_position != NULL);
//
//	int var;
//
//	switch (svType)
//	{
//	case STATIC:
//		{
//			// Get current search tree level.
//			int level = m_levelInfo->level;
//
//			var = m_ordering[level + 1];
//
//			break;
//		}
//	case DYNAMIC:
//		{
//			break;
//		}
//	}
//
//	return var;
//}
//
//// This function selects the next value for the current variable.
//int CProblemBayes::selectValue(int var, int& val, bool prune)
//{
//	// Safety checks.
//	assert(m_assignment != NULL);
//
//	int result = SV_FAILURE;		// Result of the function.
//
//	int* domainSize = m_levelInfo->domainSize;
//	double* domainCost = m_levelInfo->domainCost;
//	bool* domainState = m_levelInfo->domainState;
//	
//	assert(domainSize != NULL);
//	assert(domainCost != NULL);
//	assert(domainState != NULL);
//
//	// Select next best value.
//	int bestValue = -1;
//	double bestCost = -1.0;
//	int baseAddr = var * m_K;
//
//	int K = getStaticDomainSize(var);
//	for (int k = 0; k < K; ++k)
//	{
//		int addr = baseAddr + k;
//		if (false == domainState[addr])
//			continue;	// skip already visited values.
//
//		int value = k;
//		double cost = domainCost[addr];
//
//		if (prune)
//		{
//			if (cost <= m_bestSolutionCost)
//			{
//				// Current assignment cannot be extended 
//				// to a better solution. Disable the value.
//				domainState[addr] = false;
//				continue;	
//			}
//			else
//			{
//				// Current assignment is a valid one.
//				if (cost > bestCost)
//				{
//					bestValue = value;
//					bestCost = cost;
//				}
//			}
//		}
//		else
//		{
//			bestValue = value;
//			break;
//		}
//	}
//
//	// Next value selection.
//	val = bestValue;
//
//	if (-1 == val)
//	{
//		// No value found.
//		result = SV_NOVALUE;
//	}
//	else
//	{
//		// Found a value. Mark it as used.
//		domainState[baseAddr + val] = false;
//
//		result = SV_SUCCESS;
//	}
//
//	return result;
//}
//
//// This function runs a Branch and Bound search algorithm. It is guided
//// by a set of pre-computed heuristics for each variable-value pair. Assumes
//// static variable ordering and dynamic value ordering.
//int CProblemBayes::execBBMBs(int ibound, long tmLimit, double& tmCpuSearch,
//	long& backtracks, long& expansions, bool silent)
//{
//	// Safety checks.
//	assert(m_N > 0);
//	assert(m_K > 0);
//	assert(m_C > 0);
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
//	for(int i=0; i<m_N; i++)
//		m_bestSolution[i] = -1;
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
//// This function runs a Branch and Bound search algorithm. It is guided
//// by a set of heuristics that are dynamically computed using MBE. Assumes
//// static variable ordering and dynamic value ordering.
//int CProblemBayes::execBBMBd(int ibound, long tmLimit, double& tmCpuSearch,
//	long& backtracks, long& expansions, bool silent)
//{
//	// Safety checks.
//	assert(m_N > 0);
//	assert(m_K > 0);
//	assert(m_C > 0);
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
//	for(int i=0; i<m_N; i++)
//		m_bestSolution[i] = -1;
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
//// This function runs a Backtracking search algorithm. It has no heuristic
//// guidance and assumes static variable ordering. It is used for comparison
//// with the more advanced AND-OR search (without pruning) algorihm.
//int CProblemBayes::execBT(long tmLimit, double& tmCpuSearch,
//	long& backtracks, long& expansions, bool silent)
//{
//	// Safety checks.
//	assert(m_N > 0);
//	assert(m_K > 0);
//	assert(m_C > 0);
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
//	sprintf(tmpbuf, "\nBT started ...");
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
//	for(int i=0; i<m_N; i++)
//		m_bestSolution[i] = -1;
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
//		switch (selectValue(var, val, false))
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
//
//	return result;
//}
//
//
//void CProblemBayes::removeNode(list<CAONode*>& l, CAONode* node)
//{
//	list<CAONode*>::iterator it = l.begin();
//	while (it != l.end())
//	{
//		CAONode* tmp = (*it);
//		if (tmp == node)
//		{
//			l.erase(it);
//			break;
//		}
//
//		++it;
//	}
//}
//
//bool CProblemBayes::forwardPrune(CAONode* node, double estimate)
//{
//	// Safety checks.
//	assert(node->type() == AND);
//
//	CAONode* parentOR = node->parent();			// OR parent (consistent with the legal tree).
//	assert(parentOR != NULL);
//
//	if (parentOR->isUpdated())
//	{
//		double pseudoMax = estimate * node->getG();
//		double gValueOR = parentOR->getG();
//
//		if (pseudoMax <= gValueOR)
//			return true;
//		else
//			return false;
//	}
//
//	return false;
//}
//
//// This function expands and AND-OR search tree node.
//int CProblemBayes::expand(CAONode* node, stack<CAONode*>& succ, bool heuristic, int ibound)
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
//			int pos = m_position[var];
//			CBucket* bucket = m_buckets->getBucketAt(pos);
//			assert(bucket != NULL);
//
//			// Multiply all original functions in the bucket.
//			double g = bucket->multiply();
//			node->setG(g);
//
////			sprintf(tmpbuf, "\nLOG:  -[AND] initial g-value is set to %7g", g);
////			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//
//			// Check for zero l-value. Prune under it.
//			if (0.0 == g)
//			{
//				result = E_SUCCESS;
//				zeroDeadEnds +=1;
//				break;
//			}
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
//			break;
//		}
//	case OR:	// Expand an OR node (variable).
//		{
//			// Check for heuristic computation (full AO tree expansion).
//			if (!heuristic)
//			{
//				node->setG(0.0);
//
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
//				vector<int>::iterator it = values.begin();
//				for (; it != values.end(); ++it)
//				{
//					int val = (*it);
//
//					// Check consistency of the current value.
//					// If no values have been found consistent,
//					// then set the g-value of the terminal OR
//					// node to 1, otherwise leave it to -1.
//
//					CAONode* child = new CAONode(AND, val);
//					assert(child != NULL);
//
//					// Set up links.
//					child->setParent(node);
//					node->addChild(child);
//
//					succ.push(child);
//				}
//
//				result = E_SUCCESS;
//
//				// Free temporary buffers.
//				values.clear();
//			}
//			else
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
//
//			break;
//		}
//	}
//
//	return result;
//}
//
//// This function executes AND-OR tree search (full tree expansion).
//int CProblemBayes::execAO(long tmLimit, double& tmCpuSearch, long& expansions, long& nodesAND, long& nodesOR, bool silent)
//{
//	// Safety checks.
//	assert(m_N > 0);
//	assert(m_K > 0);
//	assert(m_C > 0);
//	assert(m_tree != NULL);
//
//	// DFS stacks.
//	stack<CAONode*> succ;
//	stack<CAONode*> OPEN;
//	list<CAONode*> CLOSED;
//	
//	// General variables.
//	int result;						// Output code.
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
//	sprintf(tmpbuf, "\nAO started ...");
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
//	expansions = 0;
//	nodesAND = 0;
//	nodesOR = 0;
//	zeroDeadEnds = 0;
//	m_probabilityOfQuery = -1;
//
//	// Set first variable.
//	int var0; var0 = m_ordering[0];
//	CAONode* root = new CAONode(OR, var0);
//	assert(root != NULL);
//
//	// Initialize the variable stack.
//	OPEN.push(root);
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
//				sprintf(tmpbuf, "\n    -int: time %7g, nodes %d", current_time, expansions);
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
//		CAONode* node = OPEN.top();
//		OPEN.pop();
//		CLOSED.push_front(node);
//
//		switch (expand(node, succ, false, UNKNOWN))
//		{
//		case E_SUCCESS:
//			{
//				// Add successors to stack.
//				while (!succ.empty())
//				{
//					CAONode* child = succ.top();
//					succ.pop();
//
//					OPEN.push(child);
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
//			switch (cand->type())
//			{
//			case AND:
//				{
//					// Update parent's g-value.
//					//double g = max(cand->getG(), parent->getG());
//					double g = (cand->getG() + parent->getG());
//					parent->setG(g);
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
//
//						// prune the other children of parent, and remove from OPEN
//						if (g == 0) 
//						{
//							for (int i=0; i < parent->childrenSize() - 1 ; i++)
//							{
//								CAONode* eraseNode = OPEN.top();
//								OPEN.pop();
//								delete eraseNode;
//							}
//							parent->clear();
//						}
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
//			m_probabilityOfQuery = root->getG();
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
//			sprintf(tmpbuf, "\n   probability of the query: %7g", m_probabilityOfQuery);
//			cout << tmpbuf;			
//			if (m_outputFile) fprintf(m_outputFile, tmpbuf);
//		
//		//	sprintf(tmpbuf, "\n   optimal solution cost: %7g", m_bestSolutionCost);
//		//	cout << tmpbuf;
//		//	if (m_outputFile) fprintf(m_outputFile, tmpbuf);
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
//
//	return result;
//}
//
//// This function executes AND-OR tree search with MBE pruning (Alpha-Beta pruning style).
//int CProblemBayes::execAOMB(int ibound, long tmLimit, double& tmCpuSearch, 
//	long& expansions, long& nodesAND, long& nodesOR, bool silent)
//{
//	// Safety checks.
//	assert(m_N > 0);
//	assert(m_K > 0);
//	assert(m_C > 0);
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
//int CProblemBayes::execBE_AUX(double& timeCpuBe, double& result)
//{
//	double* temp = NULL;
//	int k;
//	clock_t c_start;				// Start timer.
//	clock_t c_end;					// End timer.
//
//	c_start = clock();
//
//	m_buckets->process(MBE_SIMPLE, m_N, k, temp);
//
//	result = 0;
//	for(int i = 0; i < k; i++)
//		result += temp[i];
//
//	c_end = clock();
//	timeCpuBe = (double)(c_end - c_start) / 1000.0;
//
//	delete[] temp;
//
//	return S_SUCCESS;
//}