// TreeDecompositionSeparator.cpp: implementation of the CTreeDecompositionSeparator class.
//
//////////////////////////////////////////////////////////////////////

#include "TreeDecompositionSeparator.h"
#include "TreeDecomposition.h"
#include "defs.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CTreeDecompositionSeparator::CTreeDecompositionSeparator()
{
//	m_variables = NULL;
	m_clusterParent = NULL;
	m_clusterChild = NULL;
	m_clusterIDParent = -1;
	m_clusterIDChild = -1;
	m_weight = -1;
//	m_weightParent = -1;
//	m_weightChild = -1;
}

CTreeDecompositionSeparator::CTreeDecompositionSeparator(variable_s& variables)
{
	copy(variables.begin(), variables.end(), inserter(m_variables, m_variables.begin()));
	m_clusterParent = NULL;
	m_clusterChild = NULL;
	m_clusterIDParent = -1;
	m_clusterIDChild = -1;
	m_weight = -1;
}

CTreeDecompositionSeparator::~CTreeDecompositionSeparator()
{
	destroy();
}

void CTreeDecompositionSeparator::destroy()
{	
	m_variables.clear();
	m_clusterParent = NULL;
	m_clusterChild = NULL;
	m_clusterIDParent = -1;
	m_clusterIDChild = -1;
	m_weight = -1;
}

bool CTreeDecompositionSeparator::isMemberOf(int var)
{	
	assert(m_variables.size()!= 0);
		
	if (m_variables.find(var) == m_variables.end())		
		return false;
	else
		return true;
}

void CTreeDecompositionSeparator::print()
{
	cout << "\n[ p= " << m_clusterIDParent << " --- c= " << m_clusterIDChild << " ]   ";
	
	cout << "(";

	variable_s::iterator it = m_variables.begin();
	if (it != m_variables.end() )
	{
		cout << (*it);
		it++;
	}
	for (; it != m_variables.end(); ++it)
	{
		int var = (*it);
		cout << "," << var ;
	}

	cout << ")";

	cout << "   weight = " << m_weight;
}

// adds one more variable to the scope
void CTreeDecompositionSeparator::addVar(int var)
{
	m_variables.insert(var);
}

// deletes variable var from the scope
void CTreeDecompositionSeparator::deleteVar(int var)
{
	m_variables.erase(var);
}

// compute the weight of a separator
// function to compute the separtor weight, to decide where to break the tree decomposition
// can only be called when the tree decomposition is well defined
// for each separator, it computes a < 1 weight, w_left/w_right or w_right/w_left, and picks the largest
// w_right = Sum_{i=1}^{log(number_of_clusters)}  [higest_w(cluster)],
// where cluster is in the right portion after the separator is eliminated

double CTreeDecompositionSeparator::computeWeight(int scopeSize)
{
	assert(m_clusterChild != NULL);
	assert(m_clusterParent != NULL);

	double weight;

	double weight_parent = m_owner->computeSubtreeWeight(m_clusterChild, NULL, m_variables, scopeSize);
	double weight_child  = m_owner->computeSubtreeWeight(m_owner->getRootCluster(), m_clusterChild, m_variables, scopeSize);

	assert( weight_child != 0 );
	assert( weight_parent != 0 );

//	if ( weight_parent < weight_child )
//		weight = weight_parent / weight_child;
//	else
//		weight = weight_child / weight_parent;
//	m_weight = weight;

	weight = 1.0 / max (weight_child, weight_parent);
	m_weight = weight;

	return weight;
}

void CTreeDecompositionSeparator::eliminateVar(int var)
{
	variable_s::iterator it = m_variables.begin();
	for ( ; it != m_variables.end() ; it++ )
		if ( (*it) == var )
		{
			m_variables.erase(it);
			break;
		}
}