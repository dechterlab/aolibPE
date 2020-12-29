// minibucket.cpp -- Mini-Bucket structure.

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


#include "minibucket.h"
#include "problem.h"
#include "bucketstruct.h"
#include "cpt.h"

//////////////////////////////////////////////////////////////////////////////
// CMiniBucket class implementation.

CMiniBucket::CMiniBucket(int id /*NOID*/, CBucket* owner /*NULL*/)
{
	m_id = id;
	m_owner = owner;
}

CMiniBucket::~CMiniBucket()
{
	m_functions.clear();
	m_variables.clear();
}

int CMiniBucket::size()
{
	return m_functions.size();
}

void CMiniBucket::add(CFunction* c)
{
	// Safety checks.
	assert(c != NULL);

	m_functions.push_back(c);

	// Get problem instance.
	CBucket* bucket = m_owner;
	assert(bucket != NULL);

	CBucketStruct* buckets = bucket->getOwner();
	assert(buckets != NULL);

	CProblem* prob = buckets->getOwner();
	assert(prob != NULL);

	// Update the variables count.
	int* argv = c->getArgv();
	int argc = c->getArgc();
	
	for (int v = 0; v < argc; ++v)
	{
		int var = argv[v];

		if (prob->isAssigned(var))
			continue;	// skip evidence.

		if (m_variables.find(var) == m_variables.end())
			m_variables.insert(var);
	}
}

void CMiniBucket::remove(CFunction* c)
{
	// Safety checks.
	assert(c != NULL);

	function_v::iterator it;
	it = find(m_functions.begin(), m_functions.end(), c);
	if (it != m_functions.end())
		m_functions.erase(it);
}

// This function checks whether the number of different
// variables within the mini-bucket exceeds the "ibound",
// after inserting the new scope.
bool CMiniBucket::isScopeInBound(int ibound, int argc, int* argv)
{
	int cnt = (int) m_variables.size();

	// Get problem instance.
	CBucket* bucket = m_owner;
	assert(bucket != NULL);

	CBucketStruct* buckets = bucket->getOwner();
	assert(buckets != NULL);

	CProblem* prob = buckets->getOwner();
	assert(prob != NULL);

	// Check each new variable
	for (int v = 0; v < argc; ++v)
	{
		int var = argv[v];

		if (prob->isAssigned(var))
			continue;	// skip evidence.

		if (m_variables.find(var) == m_variables.end())
			cnt++;	// new variable.

		if (cnt > ibound)
			return false;
	}

	return true;
}

// This function checks if the mini-bucket is empty.
bool CMiniBucket::isEmpty()
{
	return m_functions.empty();
}

// This function processes the mini-bucket by eliminating the
// bucket variable. The MPE operations are PROD/MAX.
bool CMiniBucket::process(CFunction*& fun)
{
	// Safety checks.
	assert(size() > 0);
	
	CBucket* bucket = getOwner();
	assert(bucket != NULL);

	CBucketStruct* buckets = bucket->getOwner();
	assert(buckets != NULL);

	CProblem* prob = buckets->getOwner();
	assert(prob != NULL);

	int elim = bucket->getVariable();
	return prodSum(elim, prob, fun);	
}

bool CMiniBucket::prodMax(int elim, CProblem* prob, CFunction*& fun)
{
	// Safety checks.
	assert(prob != NULL);

	// Initialize.
	fun = NULL;

	// We need to devise the separator; that is the scope of the
	// new function that is to be recorded.

	variable_v sep;
	variable_s::iterator it = m_variables.begin();
	for (; it != m_variables.end(); ++it)
	{
		int var = (*it);
		if (var != elim)
			sep.push_back(var);
	}

	// If there is no separator, then a constant must be generated.
	if (sep.empty())
	{
		// A constant is recorded.
		CProbabilityTable* nc = new CProbabilityTable(NOID, CPT, 0, NULL);
		nc->setOwner(prob);
		nc->setOriginal(false);
		nc->setConstant(true);
		nc->init();		// special init if constant.

		// Backup current assignment.
		prob->backupAssignment();

		// The maximum value that must be recorded.
		double maxVal = -1;

		for (int k = 0; k < prob->getStaticDomainSize(elim); ++k)
		{
			prob->setValue(elim, k);

			// Do the join of all functions on current scope combination.
			double crtVal, prodVal = 1.0;
			function_v::iterator it = m_functions.begin();
			for (; it != m_functions.end(); ++it)
			{
				CFunction* c = (*it);
				c->getCurrentValue(&crtVal);

				prodVal *= crtVal;
			}

			if (prodVal > maxVal)
				maxVal = prodVal;
		}

		// Save current entry.
		nc->setCurrentValue(&maxVal);

		// Restore current assignment.
		prob->restoreAssignment();

		// Record the new function.
		fun = nc;
	}
	else
	{
		// Some variables.
		int i, j, k, s_pos = 0;

		// "Separator" temporary buffers.
		int  s_argc = sep.size();
		int* s_argv = new int[s_argc];
		int* s_vals = new int[s_argc];

		// Create new scope.
		variable_v::iterator its;
		for (its = sep.begin(); its != sep.end(); ++its) s_argv[s_pos++] = (*its);
		sep.clear();

		// Create new table.
		CProbabilityTable* nc = new CProbabilityTable(NOID, CPT, s_argc, s_argv);
		nc->setOwner(prob);
		nc->setOriginal(false);
		nc->init();

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

			double maxVal = -1;

			// Enumerate "eliminator" values.
			for (k = 0; k < prob->getStaticDomainSize(elim); ++k)
			{
				prob->setValue(elim, k);

				// Do the join of all functions on current scope combination.
				double crtVal, prodVal = 1.0;
				function_v::iterator it = m_functions.begin();
				for (; it != m_functions.end(); ++it)
				{
					CFunction* c = (*it);
					c->getCurrentValue(&crtVal);

					prodVal *= crtVal;
				}

				if (prodVal > maxVal)
					maxVal = prodVal;
			}

			// Save current entry.
			nc->setCurrentValue(&maxVal);
		}

		// Clear temporary buffers.
		delete[] s_vals;

		// Restore current assignment.
		prob->restoreAssignment();

		// Record the new function.
		fun = nc;
	}

	return true;
}


bool CMiniBucket::prodSum(int elim, CProblem* prob, CFunction*& fun)
{
	// Safety checks.
	assert(prob != NULL);

	// Initialize.
	fun = NULL;

	// We need to devise the separator; that is the scope of the
	// new function that is to be recorded.

	variable_v sep;
	variable_s::iterator it = m_variables.begin();
	for (; it != m_variables.end(); ++it)
	{
		int var = (*it);
		if (var != elim)
			sep.push_back(var);
	}

	// If there is no separator, then a constant must be generated.
	if (sep.empty())
	{
		// A constant is recorded.
		CProbabilityTable* nc = new CProbabilityTable(NOID, CPT, 0, NULL);
		nc->setOwner(prob);
		nc->setOriginal(false);
		nc->setConstant(true);
		nc->init();		// special init if constant.

		// Backup current assignment.
		prob->backupAssignment();

		// The sum of the values
		double sumVal = 0;

		for (int k = 0; k < prob->getStaticDomainSize(elim); ++k)
		{
			prob->setValue(elim, k);

			// Do the join of all functions on current scope combination.
			double crtVal, prodVal = 1.0;
			function_v::iterator it = m_functions.begin();
			for (; it != m_functions.end(); ++it)
			{
				CFunction* c = (*it);
				c->getCurrentValue(&crtVal);

				prodVal *= crtVal;
				if ( prodVal == 0 )
					break;
			}

			sumVal += prodVal;
		}

		// Save current entry.
		nc->setCurrentValue(&sumVal);

		// Restore current assignment.
		prob->restoreAssignment();

		// Record the new function.
		fun = nc;
	}
	else
	{
		// Some variables.
		int i, j, k, s_pos = 0;

		// "Separator" temporary buffers.
		int  s_argc = sep.size();
		int* s_argv = new int[s_argc];
		int* s_vals = new int[s_argc];

		// Create new scope.
		variable_v::iterator its;
		for (its = sep.begin(); its != sep.end(); ++its) s_argv[s_pos++] = (*its);
		sep.clear();

		// Create new table.
		CProbabilityTable* nc = new CProbabilityTable(NOID, CPT, s_argc, s_argv);
		nc->setOwner(prob);
		nc->setOriginal(false);
		nc->init();

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

			double sumVal = 0;

			// Enumerate "eliminator" values.
			bool flagPruneZeros = false;
			for (k = 0; k < prob->getStaticDomainSize(elim); ++k)
			{
				prob->setValue(elim, k);

				// Do the join of all functions on current scope combination.
				double crtVal, prodVal = 1.0;
				function_v::iterator it = m_functions.begin();
				for (; it != m_functions.end(); ++it)
				{
					CFunction* c = (*it);
					c->getCurrentValue(&crtVal);

					if ( (crtVal == 0) && !c->isMemberOf(elim) )
						flagPruneZeros = true;

					prodVal *= crtVal;
					if ( prodVal == 0 )
					{
						break;
					}
				}

				sumVal += prodVal;
			}

			// Save current entry.
			nc->setCurrentValue(&sumVal);

//			if (sumVal == 0)
//				cout << "cucu" << endl;

/*			// robert code (to check for 0);
			if ( (sumVal == 0) && flagPruneZeros )
			{
				cout << "lala" << endl;
				int var = s_argv[i];
				int lastVal = prob->getStaticDomainSize(var) - 1;
				if (s_vals[i] < lastVal)
					++s_vals[i];
				else
				{
					for ( k = i ; k <=s_argc-1 ; k++)
					{
						int var1 = s_argv[k];
						int lastVal1 = prob->getStaticDomainSize(var1) - 1;
						s_vals[k] = lastVal1;
					}
				}
			}
*/

		
		}

		// Clear temporary buffers.
		delete[] s_vals;

		// Restore current assignment.
		prob->restoreAssignment();

		// Record the new function.
		fun = nc;
	}

	return true;
}

