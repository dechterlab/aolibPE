// bucket.cpp -- Bucket structure.

/*
 * Copyright (C) 2003 Radu Marinescu
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


#include "bucket.h"
#include "minibucket.h"
#include "function.h"
#include "problem.h"
#include "bucketstruct.h"

CBucket::CBucket(int var)
{
	m_variable = var;
	m_owner = NULL;
	m_output = NULL;
	m_index = NOID;
}

CBucket::~CBucket()
{
	destroy();
}

// This function destroys the bucket and all
// data structures associated with it.
void CBucket::destroy()
{
	destroyPartition();
	m_augment.clear();
	m_output = NULL;
	
	function_v::iterator it = m_functions.begin();
	for (; it != m_functions.end(); ++it)
	{
		CFunction* fun = (*it);
		if (!fun->isOriginal())
		{
			long ref = fun->_release();
			if (ref < 0)
			{
				delete fun;
				fun = NULL;
			}
		}
		else
			fun->setUsed(false);
	}

	m_functions.clear();	
}

void CBucket::clean()
{
	destroyPartition();
	m_augment.clear();
	m_output = NULL;

	function_v::iterator it = m_functions.begin();
	while (it != m_functions.end())
	{
		CFunction* fun = (*it);
		if (!fun->isOriginal())
		{
			it = m_functions.erase(it);
			delete fun;
		}
		else
		{
			++it;
		}
	}
}

int CBucket::size()
{
	return m_functions.size();
}

void CBucket::print()
{
	cout << "\n - [" << m_variable << "]: ";
	function_v::iterator it = m_functions.begin();
	for (; it != m_functions.end(); ++it)
	{
		CFunction* c = (*it);
		c->print(true, false);
		cout << " ";
	}
}

// This function creates a "greedy" partition of mini-buckets.
// Each mini-bucket contains at most "i-bound" variables.
void CBucket::createPartition(int ibound)
{
	// Clear previous partition, if any.
	if (m_minibuckets.size())
		destroyPartition();

	// Sort the bucket's functions on scope size (descending).
	sort(m_functions.begin(), m_functions.end(), varity_desc);

	// Create the partition.
	function_v::iterator itf = m_functions.begin();
	for (; itf != m_functions.end(); ++itf)
	{
		// Look for a mini-bucket that could fit the function.
		CFunction* c = (*itf);
		int* argv = c->getArgv();
		int argc = c->getArgc();

		bool bFound = false;
		minibucket_v::iterator itb = m_minibuckets.begin();
		for (; itb != m_minibuckets.end(); ++itb)
		{
			CMiniBucket* mb = (*itb);
			if (mb->isScopeInBound(ibound, argc, argv))
			{
				mb->add(c);
				bFound = true;

				break;	// done.
			}
		}

		// If no mini-bucket was found.
		if (!bFound)
		{
			// Create a new mini-bucket
			int mbcnt = (int) m_minibuckets.size();
			CMiniBucket* mb = new CMiniBucket(mbcnt, this);
			mb->add(c);

			m_minibuckets.push_back(mb);
		}
	}
}

// Clear a partition of mini-buckets.
void CBucket::destroyPartition()
{
	minibucket_v::iterator it = m_minibuckets.begin();
	for (; it != m_minibuckets.end(); ++it)
	{
		CMiniBucket* mb = (*it);
		if (NULL != mb)
			delete mb;
	}
	m_minibuckets.clear();
}

void CBucket::add(CFunction* c)
{
	// Safety checks.
	assert(c != NULL);

	m_functions.push_back(c);
}

void CBucket::remove(CFunction* c)
{
	// Safety checks.
	assert(c != NULL);

	function_v::iterator it;
	it = find(m_functions.begin(), m_functions.end(), c);
	if (it != m_functions.end())
		m_functions.erase(it);
}

void CBucket::augment(CFunction* c)
{
	// Safety checks.
	assert(c != NULL);

	m_augment.push_back(c);
}

bool CBucket::evalStaticHeuristic(int k, double*& costs)
{
	// Get problem instance.
	CBucketStruct* buckets = m_owner;
	assert(buckets != NULL);

	CProblem* prob = buckets->getOwner();
	assert(prob != NULL);

	int var = m_variable;

	// Check for evidence.
	if (prob->isEvidence(var))
	{
		double cost = 1.0;
		int i = prob->getValue(var);

		function_v::iterator it = m_augment.begin();
		for (; it != m_augment.end(); ++it)
		{
			CFunction* fun = (*it);

			double p;
			fun->getCurrentValue(&p);

			cost *= p;
		}

		costs[i] = cost;
	}
	else
	{
		// Backup current assignment.
		prob->backupAssignment();

		for (int i = 0; i < k; ++i)
		{
			// Set current value.
			prob->setValue(var, i);

			double cost = 1.0;
			function_v::iterator it = m_augment.begin();
			for (; it != m_augment.end(); ++it)
			{
				CFunction* fun = (*it);

				double p;
				fun->getCurrentValue(&p);

				cost *= p;
			}

			costs[i] = cost;
		}

		// Restore current assignment.
		prob->restoreAssignment();
	}

	return true;
}

bool CBucket::evalDynamicHeuristic(int ibound, int k, double*& costs)
{
	return true;
}

// This function process the current bucket. It either passes along
// the functions if the bucket's variable is evidence or approximates
// the functions in the bucket via mini-buckets.
bool CBucket::process(int type, int ibound, int bottom)
{
	// Get problem instance.
	CBucketStruct* buckets = m_owner;
	assert(buckets != NULL);

	CProblem* prob = buckets->getOwner();
	assert(prob != NULL);

	int highestVar, highestPos;

	// Check if the bucket's variable is evidence.
	if (prob->isEvidence(m_variable))
	{
		function_v::iterator it = m_functions.begin();
		for (; it != m_functions.end(); ++it)
		{
			CFunction* fun = (*it);
			
			if (0 != m_index)
			{
				CFunction* pt = NULL;

				// Check for constant.
				if (fun->isConstant())
				{
					highestPos = bottom;
				}
				else
				{
					// Substitute evidence and record a new function.
					int value = prob->getValue(m_variable);
					
					fun->substitute(m_variable, value, pt);
					
//					if (!pt->verify())
//						int bbb = 1;

					if (pt->isConstant())
					{
						highestPos = bottom;
					}
					else
					{
						buckets->getHighestNode(pt, highestVar, highestPos);
						//pt->getHighestNode(highestVar, highestPos);
					}
				}

				// Destination bucket.
				CBucket* bkt = buckets->getBucketAt(highestPos);
				bkt->add(pt);
			}
		}

		return true;
	}

	// Otherwise, create mini-buckets.
	createPartition(ibound);

	// Process each mini-bucket and record the output function(s).
	minibucket_v::iterator it = m_minibuckets.begin();
	for (; it != m_minibuckets.end(); ++it)
	{
		// Current mini-bucket.
		CMiniBucket* mb = (*it);
		
		// New function to be generated.
		CFunction* fun = NULL;

		if (!mb->process(fun))
		{
			if (fun) delete fun;
			return false;
		}

		if (fun)
		{
//			if (!fun->verify())
//				int bbb = 1;

			if (0 != m_index)
			{
				if (fun->isConstant())
				{
					highestPos = bottom;		// add constants to the first bucket.
				}
				else
				{
					buckets->getHighestNode(fun, highestVar, highestPos);
					//fun->getHighestNode(highestVar, highestPos);
				}

				// Destination bucket.
				CBucket* bkt = buckets->getBucketAt(highestPos);
				bkt->add(fun);

			}
			else
			{
				m_output = fun;
			}
		}
	}

	return true;
}

// This function multiplys all functions in the bucket.
// Variables in the bucket must be all instantiated except
// the bucket's variable (if it is not evidence).
bool CBucket::solve(int k, double*& bounds)
{
	// Get problem instance.
	CBucketStruct* buckets = m_owner;
	assert(buckets != NULL);

	CProblem* prob = buckets->getOwner();
	assert(prob != NULL);

	// Backup current assignment.
	prob->backupAssignment();

	// Bucket's variable.
	int var = m_variable;

	for (int val = 0; val < k; ++val)
	{
		// Set current value.
		if (!prob->isEvidence(var))
			prob->setValue(var, val);

		double cost = 1;
		function_v::iterator it = m_functions.begin();
		for (; it != m_functions.end(); ++it)
		{
			CFunction* fun = (*it);

			double crt;
			fun->getCurrentValue(&crt);

			cost *= crt;
		}

		bounds[val] = cost;
	}

	// Restore current assignment.
	prob->restoreAssignment();

	return true;
}

// This function multiplies all original functions in the bucket.
// It requires that all arguments have been previously instantiated.
double CBucket::multiply()
{
	double cost = 1.0;

	function_v::iterator it = m_functions.begin();
	for (; it != m_functions.end(); ++it)
	{
		CFunction* fun = (*it);
		if (fun->isOriginal())
		{
			double crt;
			fun->getCurrentValue(&crt);
			
			double temp = cost*crt;
			if (cost != 0 && crt != 0  &&  temp == 0)
				cost = exp( -500.0 );
			else
				cost = temp;
		}
	}

	return cost;
}
