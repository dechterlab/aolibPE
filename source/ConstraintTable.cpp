// ConstraintTable.cpp: Constraint Table

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


#include "ConstraintTable.h"
#include "problem.h"
#include "constraintnode.h"

CConstraintTable::CConstraintTable() : CFunction()
{
}

CConstraintTable::CConstraintTable(int id, CTYPE type, int argc, int* argv)
	: CFunction(id, type, argc, argv)
{
	m_tableSize = 0;
	m_table = NULL;
	m_highestVarInScope = UNKNOWN;
	m_rootShannonTree = NULL;
	m_argvShannonTree = NULL;
}


CConstraintTable::~CConstraintTable()
{
	if (m_table)
		delete[] m_table;
	if (m_argvShannonTree)
		delete[] m_argvShannonTree;
}

///////////////////////////////////////////////////////////////////
/////////////////////                         /////////////////////
/////////////////////   copied from cpt.cpp   /////////////////////
/////////////////////                         /////////////////////
///////////////////////////////////////////////////////////////////

int CConstraintTable::getValueAt(int addr)
{
	// Safety checks.
	assert(m_table != NULL);
	assert(addr < m_tableSize);

	return m_table[addr];
}

void CConstraintTable::setValueAt(int addr, int val)
{
	// Safety checks.
	assert(m_table != NULL);
	assert(addr < m_tableSize);

	m_table[addr] = val;
}

// This function returns the table entry, based on
// the current scope assignment. Requires that 
// all variables in the constraint's scope have
// already been instantiated.
int CConstraintTable::getCurrentValue(void* val)
{
	assert(m_table != NULL);

	int addr = (m_constant) ? 0 : getCurrentAddress();
	int v = m_table[addr];

	*((int*) val) = v;

	return addr;
}

// This function sets the table entry, based on
// the current scope assignment. Requires that 
// all variables in the constraint's scope have
// already been instantiated.
int CConstraintTable::setCurrentValue(void* val)
{
	assert(m_table != NULL);
	
	int addr = (m_constant) ? 0 : getCurrentAddress();
	int v = *((int*) val);
		
	m_table[addr] = v;

	return addr;
}

// This function allocates/initializes an empty soft
// constraint. Requires the scope to have been already 
// initialized.
void CConstraintTable::init()
{
	if (!m_constant)
	{
		assert(m_argv != NULL);
		assert(m_argc > 0);

		m_tableSize = computeTableSize();
		m_table = new int[m_tableSize];
		for(int i=0; i<m_tableSize; i++)
			m_table[i] = 0;	
	}
	else
	{
		m_tableSize = 1;
		m_table = new int[m_tableSize];
		for(int i=0; i<m_tableSize; i++)
			m_table[i] = 0;	
	}
}

void CConstraintTable::destroy()
{
	if (m_table)
		delete[] m_table;

	m_table = NULL;
	m_tableSize = 0;
	m_constant = false;
	destroyShannonTree();

	if (m_argvShannonTree)
		delete[] m_argvShannonTree;
}

/*
int CConstraintTable::getChild()
{
	// Safety checks.
	
	assert(m_argv != NULL);
	assert(m_argc > 0);

	return m_argv[m_argc - 1];
}
*/

int CConstraintTable::getParent(int idx)
{
	// Safety checks.
	assert(m_argv != NULL);
	assert(m_argc > 0);
	assert((idx >= 0) && (idx < (m_argc - 1)));

	return m_argv[idx];
}

int CConstraintTable::getNumberOfParents()
{
	// Safety checks.
	assert(m_argv != NULL);
	assert(m_argc > 0);

	return (m_argc - 1);
}

// This function substitutes a variable in the scope with
// specified value. The result is a new (non-original) function.
void CConstraintTable::substitute(int var, int val, CFunction*& fun)
{
	// Safety checks.
	assert(m_argv != NULL);
	assert(m_argc > 0);
	assert(isMemberOf(var));

	// Get problem instance.
	CProblem* prob = m_owner;
	assert(prob != NULL);

	// Compute new scope.
	variable_v scope;
	for (int i = 0; i < m_argc; ++i)
	{
		int arg = m_argv[i];
		if (arg != var)
		{
			scope.push_back(arg);
		}
	}

	if (scope.empty())
	{
		// A constant is recorded.
		CConstraintTable* pt = new CConstraintTable(NOID, CPT, 0, NULL);
		pt->setOwner(prob);
		pt->setOriginal(false);
		pt->setConstant(true);
		pt->init();		// special init if constant.

		// Backup current assignment.
		prob->backupAssignment();

		if (!prob->isEvidence(var))
			prob->setValue(var, val);

		double crtVal;
		getCurrentValue(&crtVal);

		// Save current entry.
		pt->setCurrentValue(&crtVal);

		// Restore current assignment.
		prob->restoreAssignment();

		// Record the new function.
		fun = pt;
	}
	else
	{
		// Some variables.
		int i, k, j, s_pos = 0;

		// "Separator" temporary buffers.
		int  s_argc = scope.size();
		int* s_argv = new int[s_argc];
		int* s_vals = new int[s_argc];

		// Create new scope.
		variable_v::iterator its;
		for (its = scope.begin(); its != scope.end(); ++its) s_argv[s_pos++] = (*its);
		scope.clear();

		// Create new table.
		CConstraintTable* pt = new CConstraintTable(NOID, CPT, s_argc, s_argv);
		pt->setOwner(prob);
		pt->setOriginal(false);
		pt->init();

		// Backup current assignment.
		prob->backupAssignment();

		// Reset "separator" variables.
		for (k = 0; k < s_argc; ++k) 
		{
			s_vals[k] = 0;
		}
		s_vals[s_argc - 1] = -1;

		// Start elimating variables.
		while (1)
		{
			// Enumerate "separator" variables.
			for (i = s_argc - 1; i >= 0; --i)
			{
				int var = s_argv[i];
				int lastVal = prob->getStaticDomainSize(var) - 1;
				if (s_vals[i] < lastVal) break;
				s_vals[i] = 0;
			}

			if (i < 0) break;	// done;
			++s_vals[i];
			
			// NOW: all "separator" arguments have a specific value combination.
			for (j = 0; j < s_argc; ++j) prob->setValue(s_argv[j], s_vals[j]);

			// Substitute variable (eliminator).
			if (!prob->isEvidence(var))
				prob->setValue(var, val);

			double crtVal;
			getCurrentValue(&crtVal);

			// Save current entry.
			pt->setCurrentValue(&crtVal);
		}

		// Clear temporary buffers.
		delete[] s_vals;

		// Restore current assignment.
		prob->restoreAssignment();

		// Record the new function.
		fun = pt;
	}
}


void CConstraintTable::create(int tableSize, int* table)
{
	init();
	m_tableSize = tableSize;
	for (int i = 0; i < m_tableSize; ++i) m_table[i] = table[i];
}

// This function creates the constraint table based on tightness (% of allowed tuples)
void CConstraintTable::create(double tightness)
{
	// Safety checks.
	assert(m_argc > 0);
	assert(m_argv != NULL);

	// Allocate the probability table.
	init();

	// Get the problem instance.
	CProblem* prob = m_owner;
	assert(prob != NULL);

	//Find total number of possible tuples
	int i;
	int ntuples = 1;
	for (i=0; i<m_argc; i++)
		ntuples *= prob->getStaticDomainSize(m_argv[i]);
	
	//Decide tuples which are allowed, according to tightness
	int cflag = 0; //at least one tuple should be allowed
	for (i=0; i<ntuples; i++)
	{
		double p;
		p = RandUniformDouble();
		if (p <= tightness)
		{
			m_table[i] = 1;
			cflag = 1;
		}
		else
			m_table[i] = 0;
	}
//	if (cflag == 0)
//	{
//		long t = RandUniform(ntuples);
//		m_table[t] = 1;
//	}
	// make sure the problem has at least one solution ([0,0,...,0])
	// m_table[0] = 1;
}


bool CConstraintTable::verify()
{
	assert(m_table != NULL);
	assert(m_tableSize > 0);

	for (int i = 0; i < m_tableSize; ++i)
	{
		if (m_table[i] < 0)
			return false;
	}

	return true;
}

void CConstraintTable::print(bool scope, bool table)
{
	if (scope)
	{
		cout << endl;
		cout << "(";
		int i;
		for (i = 0; i < m_argc - 1; ++i)
			cout << m_argv[i] << ",";
		cout << m_argv[i] << ")";
	}

	if (table)
	{
		cout << endl;
		for (int i=0; i < m_tableSize; i++)
			cout << "["<<i<<"] = "<< m_table[i] << endl;
	}
}

// this creates a Shannon tree to represent the constraint table; it is used in constraint propagation 
// (e.g., unit resolution, forward checking)
// the tree maintains at each node a 0 if no solutions below, or a 1 if at least one solution.
// the order of variables in the tree is according to the ordering of the problem!
void CConstraintTable::createShannonTree()
{
	// safety checks
	assert(m_argc > 0);
	assert(m_argv != NULL);

	// get problem instance (we need the ordering)
	CProblem* prob = m_owner;
	assert(prob != NULL);

	// create the ordered scope for the Shannon Tree (m_argvShannonTree)
	m_argvShannonTree = new int[m_argc];
	for(int i=0; i < m_argc; i++)
		m_argvShannonTree[i] = m_argv[i];
	

	// order argv according to problem ordering
	int* ordering = prob->getOrdering();

	// bubble sort the argv
	for(int i=0; i < m_argc; i++)
	{
		for(int j=i+1; j < m_argc; j++)
		{
			if( ordering[m_argvShannonTree[j]] < ordering[m_argvShannonTree[i]])
			{
				int temp = m_argvShannonTree[i];
				m_argvShannonTree[i] = m_argvShannonTree[j];
				m_argvShannonTree[j] = temp;
			}
		}
	}

	// build the Shannon Tree by dfs traversal
	stack<CConstraintNode*> dfsStack;	
	
	// create dummy root for the tree
	CConstraintNode* root = new CConstraintNode(-1, -1); 
	m_rootShannonTree = root;
	
	int var = m_argvShannonTree[0];
	int domainSize = prob->getStaticDomainSize(var);
	
	m_rootShannonTree->createChildren(var,domainSize);
	

	for(int i = domainSize-1; i >= 0; i--)
	{
		dfsStack.push( m_rootShannonTree->children()[i] );
	}

	int varIndex = -1;
	while(!dfsStack.empty())
	{
		// search is at the parent of the top of dfsStack, so we advance to the variable of dfsStack.top()
		varIndex++; 

		CConstraintNode* node = dfsStack.top();
		dfsStack.pop();

		// Set current value.
		var = node->variable();
		int val = node->value();
		prob->setValue(var, val);

		
		if(varIndex < m_argc-1) 
		{
			// Expand the current node (scope not fully assigned yet)
			var = m_argvShannonTree[varIndex+1];
			domainSize = prob->getStaticDomainSize(var);
			node->createChildren(var,domainSize);
			for(int i = domainSize-1; i >= 0; i--)
			{
				dfsStack.push( node->children()[i] );
			}
		}
		else
		{
			int currentValue;
			getCurrentValue(&currentValue);
			node->setG(currentValue);

			// Propagate solution count
			bool continuePropagation = true;
			while (continuePropagation)
			{
				if(node == m_rootShannonTree)
					break;
				node->parent()->addToG( node->getG() );
				if( node->value() < prob->getStaticDomainSize(m_argvShannonTree[varIndex]) - 1 )
					continuePropagation = false;
				else
					node = node->parent();
				varIndex--;
			}
		}
	}
}

void CConstraintTable::destroyShannonTree()
{
	CProblem* prob = m_owner;
	assert(prob != NULL);

	stack<CConstraintNode*> dfsStack;	
	dfsStack.push(m_rootShannonTree);

	int varIndex = -1;
	while(!dfsStack.empty())
	{
		varIndex++; 
		CConstraintNode* node = dfsStack.top();
		dfsStack.pop();

		if (node->children() != NULL)
		{
			int domainSizeNextVar = m_argvShannonTree[varIndex+1];
			for(int i = domainSizeNextVar -1; i>=0 ; i--)
			{
				dfsStack.push(node->children()[i]);
			}
		}

		delete node;
	}			
}

void CConstraintTable::printShannonTree()
{
	CProblem* prob = m_owner;
	assert(prob != NULL);

	assert(m_rootShannonTree != NULL);

	queue<CConstraintNode*> breathQueue;	
	breathQueue.push(m_rootShannonTree);

	int var = -1;
	while(!breathQueue.empty())
	{
		CConstraintNode* node = breathQueue.front();
		breathQueue.pop();

		if(node->variable() != var)
		{
			var = node->variable();
			cout << "\n";
		}

		int val = node->value();
		double gValue = node->getG();
		cout << "[" << var << "," << val << "|" << gValue << "] ";
		
		if (node->children() != NULL)
		{
			int nextVar = node->children()[0]->variable();
			int domainSizeNextVar = prob->getStaticDomainSize( nextVar );
			for(int i = 0; i < domainSizeNextVar ; i++)
				breathQueue.push(node->children()[i]);
		}	
	}
}

