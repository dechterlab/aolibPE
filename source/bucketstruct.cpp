// bucketstruct.cpp -- Bucket structure associated with a (sub)problem.

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

#include "bucketstruct.h"

CBucketStruct::CBucketStruct(CProblem* owner)
	: m_owner(owner)
{
	m_buckets = NULL;
	m_size = 0;
	m_position = NULL;
	m_ordering = NULL;
	m_locals = false;
}

CBucketStruct::~CBucketStruct()
{
	destroy();
}

// This function initializes the bucket structure associated
// with the current problem instance.
void CBucketStruct::init(bool dfsOrdering)
{
	// Safety checks.
	assert(m_owner != NULL);

	// Local variables.
	int i;

	// Get the problem instance.
	CProblem* prob = m_owner;
	int N = prob->getN();
	m_size = N;

	// Get variable ordering/position.
	int* ordering = prob->getOrdering(dfsOrdering);
	int* position = prob->getPosition(dfsOrdering);

	// Setup local variable ordring/position.
	m_ordering = ordering;
	m_position = position;
	m_locals = false;

	// Destroy previous structure, if any.
	destroy();

	// Allocate memory.
	m_buckets = new CBucket*[N];
	assert(m_buckets != NULL);
	for (i = 0; i < N; ++i)
		m_buckets[i] = NULL;

	// Distribute the original functions into buckets.
	for (i = N - 1; i >= 0; --i)
	{
		int var = ordering[i];
		CBucket* bkt = new CBucket(var);
		m_buckets[i] = bkt;
		bkt->setOwner(this);
		bkt->setIndex(i);

		// Look for all functions that have "highest"
		// in their scope and place them in the bucket.
		function_map& funs = prob->functions();
		function_map::iterator it = funs.begin();
		for (; it != funs.end(); ++it)
		{
			CFunction* fun = (*it).second;
			if (fun->isUsed()) continue;

			if (fun->isMemberOf(var))
			{
				bkt->add(fun);
				fun->setUsed(true);
			}
		}
	}

	// Reset original functions.
	function_map& funs = prob->functions();
	function_map::iterator it = funs.begin();
	for (; it != funs.end(); ++it)
	{
		CFunction* fun = (*it).second;
		fun->setUsed(false);
	}
}

void CBucketStruct::init(int n, variable_v& ordering)
{
	// Safety checks.
	assert(m_buckets == NULL);
	assert(m_size == 0);

	// Get problem instance.
	CProblem* prob = m_owner;
	assert(prob != NULL);

	m_size = n;
	m_buckets = new CBucket*[n];
	assert(m_buckets != NULL);

	// Allocate local variable ordering/position.
	m_ordering = new int[n];
	m_position = new int[prob->getN()];
	m_locals = true;

	int pos = n - 1;
	variable_v::reverse_iterator it = ordering.rbegin();
	for (; it != ordering.rend(); ++it)
	{
		int var = (*it);
		CBucket* bkt = new CBucket(var);
		m_buckets[pos] = bkt;
		bkt->setOwner(this);
		bkt->setIndex(pos);
		m_ordering[pos] = var;
		m_position[var] = pos;

		int cnt1 = 0, cnt2 = 0;

/*
		// Look for all functions that have "highest"
		// in their scope and place them in the bucket.
		function_map& funs1 = prob->functions();
		function_map::iterator it1 = funs1.begin();
		for (; it1 != funs1.end(); ++it1)
		{
			CFunction* fun = (*it1).second;
			if (fun->isUsed()) continue;

			if (fun->isMemberOf(var))
			{
				bkt->add(fun);
				++cnt1;
				fun->setUsed(true);
			}
		}
*/

		// Populate bucket with original functions.
		function_v& funs = prob->getBucket(var);	
		function_v::iterator it = funs.begin();
		for (; it != funs.end(); ++it)
		{
			CFunction* fun = (*it);
			if (fun->isOriginal())
			{
				++cnt2;
				bkt->add(fun);
			}
		}

		// Go to next bucket.
		--pos;
	}
}

void CBucketStruct::print()
{
	cout << "\n Bucket structure:";
	for (int i = m_size - 1; i >= 0; --i)
	{
		CBucket* bkt = m_buckets[i];
		bkt->print();
	}
}

// This function destroys the bucket structure.
void CBucketStruct::destroy()
{
	if (m_buckets)
	{
		for (int i = 0; i < m_size; ++i)
			delete m_buckets[i];
		delete[] m_buckets;
		m_buckets = NULL;
		m_size = 0;
	}

	if (m_locals)
	{
		delete[] m_ordering;
		delete[] m_position;
	}
}

// This function cleans the structure. It preserves only the
// original functions in allocated buckets.
void CBucketStruct::clean()
{
	if (m_buckets)
	{
		for (int i = 0; i < m_size; ++i)
			m_buckets[i]->clean();
	}
}

CBucket* CBucketStruct::getBucketAt(int pos)
{
	// Safety checks.
	assert(m_buckets != NULL);
	assert((pos >= 0) && (pos < m_size));

	return m_buckets[pos];
}

void CBucketStruct::getHighestNode(CFunction* fun, int& highestVar, int& highestPos)
{
	// Safety checks.
	assert(m_ordering != NULL);
	assert(m_position != NULL);

	if (fun->isConstant())
	{
		highestVar = -1;
		highestPos = 0;
	}
	else
	{
		int* argv = fun->getArgv();
		int argc = fun->getArgc();

		int hPos = -1, hVar = -1;
		for (int i = 0; i < argc; ++i)
		{
			int arg = argv[i];

			int pos = m_position[arg];
			if (pos > hPos)
			{
				hPos = pos;
				hVar = arg;
			}
		}

		highestPos = hPos;
		highestVar = hVar;
	}
}

// This function processes the current bucket structure. At the lowest
// bucket we get the desired upper-bounds for each value of the variable.
bool CBucketStruct::process(int type, int ibound, int& k, double*& result, int level)
{
	// Safety checks.
	assert(m_buckets != NULL);
	assert(m_owner != NULL);

	// Clean previous round.
	if (MBE_PARTIAL == type)
	{
		clean();
	}

	// Get problem instance.
	CProblem* prob = m_owner;
	CBucket* bkt0 = m_buckets[level];
	int var0 = bkt0->getVariable();

	// Allocate space for the result.
	if (prob->isEvidence(var0))
		k = 1;
	else
		k = prob->getStaticDomainSize(var0);

	result = new double[k];
	for(int i=0; i<k; i++)
		result[i] = 0;

	// Execute Mini-Buckets approximation.
	int N = m_size;
	for (int pos = N - 1; pos > level; --pos)
	{
		CBucket* bkt = m_buckets[pos];
		if (!bkt->process(type, ibound, level))
			return false;
	}

	// At the lowest bucket we have the upper bounds.
	bkt0->solve(k, result);

	return true;
}



// Local Variables:
// mode: C++
// End:
