// TreeDecompositionCluster.cpp: implementation of the CTreeDecompositionNode class.
//
//////////////////////////////////////////////////////////////////////

#include "TreeDecomposition.h"
#include "TreeDecompositionCluster.h"
#include "defs.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CTreeDecompositionCluster::CTreeDecompositionCluster()
{
	m_id = NOID;
	m_parent = NULL;
}

CTreeDecompositionCluster::CTreeDecompositionCluster(int id, variable_s& variables)
{
	m_id = id;
	copy(variables.begin(), variables.end(), inserter(m_variables, m_variables.begin()));
	m_parent = NULL;
}

CTreeDecompositionCluster::~CTreeDecompositionCluster()
{
	destroy();
}


void CTreeDecompositionCluster::destroy()
{
	m_children.clear();
	m_variables.clear();
}

/**
 * Test if a given variable is part of the cluster's scope
 * @param [var] variable to look for
 * @return 'true' if so, 'false' otherwise
 */
bool CTreeDecompositionCluster::isMemberOf(int var)
{	
	assert(m_variables.size() != 0);
	
	if (m_variables.find(var) == m_variables.end())		
		return false;
	else
		return true;
}

void CTreeDecompositionCluster::print()
{
	cout << "[ id " << m_id << " ] --- "; // m_id is not the var of the cluster!!!
	cout << "(";

	variable_s::iterator it = m_variables.begin();
	if(it != m_variables.end() )
	{
		cout << (*it);
		it++;
	}
	for (; it != m_variables.end(); ++it)
	{
		int var = (*it);
		cout << "," << var ;
	}

	cout << ")\n";
}

// adds one more variable to the scope
void CTreeDecompositionCluster::addVar(int var)
{
	m_variables.insert(var);
}


// deletes variable var from the scope
void CTreeDecompositionCluster::deleteVar(int var)
{
	m_variables.erase(var);
}


void CTreeDecompositionCluster::addChild(CTreeDecompositionCluster* cluster)
{
	// Safety checks.
	assert(cluster != NULL);
	
	m_children.push_back(cluster);
}

int CTreeDecompositionCluster::highestVar()
{
	int* positions = m_owner_prob->getPosition();
	int maxPos = -1;

	variable_s::iterator it = m_variables.begin();
	for ( ; it != m_variables.end(); it++ )
		maxPos = max( maxPos, positions[(*it)] );

	return maxPos;
}

void CTreeDecompositionCluster::eliminateVar(int var)
{
	variable_s::iterator it = m_variables.begin();
	for ( ; it != m_variables.end() ; it++ )
		if ( (*it) == var )
		{
			m_variables.erase(it);
			break;
		}
}

