// aonode.cpp -- AND-OR search tree node.

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

#include "aonode.h"
#include "bucketstruct.h"

////////////////////////////////////////////////////////////////////////
// CAONode class implementation.

CAONode::CAONode(NTYPE type, int label)
	: m_type(type), m_label(label)
{
	m_parent = NULL;
	m_gvalue = 0;
	m_gcountvalue = 0;
	m_lvalue = 0;
	m_solved = false;
	m_cost = UNKNOWN;
	m_updated = false;
}

CAONode::~CAONode()
{
	destroy();
}

void CAONode::destroy()
{
	m_parent = NULL;
	m_children.clear();

	vector<VALUECOST*>::iterator it = m_estimates.begin();
	for (; it != m_estimates.end(); ++it)
	{
		VALUECOST* vc = (*it);
		delete vc;
	}

	m_estimates.clear();
}

void CAONode::clear()
{
	m_children.clear();
}

void CAONode::addChild(CAONode* n)
{
	// Safety checks.
	assert(n != NULL);
	
	m_children.push_back(n);
}

void CAONode::setParent(CAONode* n)
{
	m_parent = n;
}

// This function removes a child of the current node.
void CAONode::remove(CAONode* node)
{
	aonode_v::iterator it = m_children.begin();
	while (it != m_children.end())
	{
		CAONode* child = (*it);
		if (child == node)
		{
			m_children.erase(it);
			break;
		}

		++it;
	}
}

// This function removes all children of current node that have
// a heuristic estimate not better than current g-value propagated
// from the culprit node (which is a child on the current OR node).
void CAONode::prune(CAONode* culprit, double g, vector<CAONode*>& siblings)
{
	if (OR == m_type)
	{
		aonode_v::iterator it = m_children.begin();
		while (it != m_children.end())
		{
			CAONode* child = (*it);
			if (child != culprit)
			{
				double ub = child->cost();
				if (ub <= g)
				{
					it = m_children.erase(it);
					siblings.push_back(child);
					continue;
				}
			}

			++it;
		}
	}
}

void CAONode::buildValueCosts(CProblem* prob, int ibound, double& maxCost)
{
	// Safety checks.
	assert(OR == m_type);
	assert(prob != NULL);

	// Estimate value-costs for the subproblem rooted at current OR node.
	maxCost = -1.0;
	int k;
	int var = m_label;
	double* costs = NULL;
	vector<int> values;

	// Check for evidence variable.
	if (prob->isEvidence(var))
	{
		int val = prob->getValue(var);
		values.push_back(val);
	}
	else
	{
		for (int val = 0; val < prob->getStaticDomainSize(var); ++val)
		{
			values.push_back(val);
		}
	}
	
	// Create list of descendants of current OR node in legal tree.
	variable_v& descendants = prob->getTreeDescendants(var); //m_descendants[ch];
	int n = descendants.size();

	// Reset values for all descendants.
	variable_v::iterator itd = descendants.begin();
	for (; itd != descendants.end(); ++itd)
	{
		int v = (*itd);
		if (prob->isEvidence(v))
			continue;	// skip evidence.

		prob->resetValue(v);
	}

	// Create bucket structure associated with descendats list.
	CBucketStruct* buckets = new CBucketStruct(prob);
	assert(buckets != NULL);

	buckets->init(n, descendants);
	buckets->process(MBE_SIMPLE, ibound, k, costs);

	// Set up heuristic estimates.
	int i = 0;
	vector<int>::iterator it = values.begin();
	for (; it != values.end(); ++it)
	{
		int val = (*it);
		double cost = costs[val];
		
		if (cost > maxCost)
			maxCost = cost;

		m_estimates.push_back(new VALUECOST(val, cost));
	}

	// Sort in descending order the estimates.
	sort(m_estimates.begin(), m_estimates.end(), vc_desc);

	// Free temporary buffers.
	delete[] costs;
	delete buckets;
	values.clear();

}


void CAONode::computeExactSubtree(CProblem* prob, int ibound, double& sumCost)
{
	// Safety checks.
	assert(OR == m_type);
	assert(prob != NULL);

	// Estimate value-costs for the subproblem rooted at current OR node.
	sumCost = 0.0;
	int k;
	int var = m_label;
	double* costs = NULL;
	vector<int> values;

	// Check for evidence variable.
	if (prob->isEvidence(var))
	{
		int val = prob->getValue(var);
		values.push_back(val);
	}
	else
	{
		for (int val = 0; val < prob->getStaticDomainSize(var); ++val)
		{
			values.push_back(val);
		}
	}
	
	// Create list of descendants of current OR node in legal tree.
	variable_v& descendants = prob->getTreeDescendants(var); //m_descendants[ch];
	int n = descendants.size();

	// Reset values for all descendants.
	variable_v::iterator itd = descendants.begin();
	for (; itd != descendants.end(); ++itd)
	{
		int v = (*itd);
		if (prob->isEvidence(v))
			continue;	// skip evidence.

		prob->resetValue(v);
	}

	// Create bucket structure associated with descendats list.
	CBucketStruct* buckets = new CBucketStruct(prob);
	assert(buckets != NULL);

	buckets->init(n, descendants);
	buckets->process(MBE_SIMPLE, ibound, k, costs);

	// Set up heuristic estimates.
	int i = 0;
	vector<int>::iterator it = values.begin();
	for (; it != values.end(); ++it)
	{
		int val = (*it);
		double cost = costs[val];
		
		sumCost += cost;

//		m_estimates.push_back(new VALUECOST(val, cost));
	}

	// Free temporary buffers.
	delete[] costs;
	delete buckets;
	values.clear();

}

void CAONode::addInvalidValue(int var, int val)
{
	m_invalidValues.insert(make_pair(var,val));
}

void CAONode::cleanInvalidValues(bool** table, int* count, int* assignment)
{
	variable_multimap::iterator it = m_invalidValues.begin();
	for ( ; it != m_invalidValues.end() ; ++it )
	{
		int var = (*it).first;
		int val = (*it).second;

		table[var][val] = true;
		count[var]++;
		if (count[var] > 1)
			assignment[var] = -1;			
	}

	m_invalidValues.clear();
}


/*
// This function sets value estimates for all successors of an OR node.
void CAONode::setEstimates(vector<VALUECOST*>& estimates)
{
	if (OR == m_type)
	{
		copy(estimates.begin(), estimates.end(), back_inserter(m_estimates));
	}
}
*/

// Local Variables:
// mode: C++
// End:
