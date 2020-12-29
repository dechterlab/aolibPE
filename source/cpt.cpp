// cpt.cpp -- Bayesian Conditional Probability Table (aka CPT).

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


#include "cpt.h"
#include "problem.h"

CProbabilityTable::CProbabilityTable() : CFunction()
{
}

CProbabilityTable::CProbabilityTable(int id, CTYPE type, int argc, int* argv)
	: CFunction(id, type, argc, argv)
{
	m_tableSize = 0;
	m_table = NULL;
}

CProbabilityTable::CProbabilityTable(int id, CTYPE type, int argc, int* argv, int tab_size, double *tab)
	: CFunction(id, type, argc, argv), m_tableSize(tab_size), m_table(tab)
{
}

CProbabilityTable::~CProbabilityTable()
{
	if (m_table)
		delete[] m_table;
}

double CProbabilityTable::getValueAt(int addr)
{
	// Safety checks.
	assert(m_table != NULL);
	assert(addr < m_tableSize);

	return m_table[addr];
}

void CProbabilityTable::setValueAt(int addr, double val)
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
int CProbabilityTable::getCurrentValue(void* val)
{
	assert(m_table != NULL);

	int addr = (m_constant) ? 0 : getCurrentAddress();
	double v = m_table[addr];

	*((double*) val) = v;

	return addr;
}

// This function sets the table entry, based on
// the current scope assignment. Requires that 
// all variables in the constraint's scope have
// already been instantiated.
int CProbabilityTable::setCurrentValue(void* val)
{
	assert(m_table != NULL);
	
	int addr = (m_constant) ? 0 : getCurrentAddress();
	double v = *((double*) val);
		
	m_table[addr] = v;

	return addr;
}

// This function allocates/initializes an empty soft
// constraint. Requires the scope to have been already 
// initialized.
void CProbabilityTable::init()
{
	if (!m_constant)
	{
		assert(m_argv != NULL);
		assert(m_argc > 0);

		m_tableSize = computeTableSize();
		m_table = new double[m_tableSize];
		for(int i=0; i<m_tableSize; i++)
			m_table[i] = 0;	
	}
	else
	{
		m_tableSize = 1;
		m_table = new double[m_tableSize];
		for(int i=0; i<m_tableSize; i++)
			m_table[i] = 0;	
	}
}

void CProbabilityTable::destroy()
{
	if (m_table)
		delete[] m_table;

	m_table = NULL;
	m_tableSize = 0;
	m_constant = false;
}

int CProbabilityTable::getChild()
{
	// Safety checks.
	assert(m_argv != NULL);
	assert(m_argc > 0);

	return m_argv[m_argc - 1];
}


int CProbabilityTable::getParent(int idx)
{
	// Safety checks.
	assert(m_argv != NULL);
	assert(m_argc > 0);
	assert((idx >= 0) && (idx < (m_argc - 1)));

	return m_argv[idx];
}

int CProbabilityTable::getNumberOfParents()
{
	// Safety checks.
	assert(m_argv != NULL);
	assert(m_argc > 0);

	return (m_argc - 1);
}

// This function substitutes a variable in the scope with
// specified value. The result is a new (non-original) function.
void CProbabilityTable::substitute(int var, int val, CFunction*& fun)
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
		CProbabilityTable* pt = new CProbabilityTable(NOID, CPT, 0, NULL);
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
		CProbabilityTable* pt = new CProbabilityTable(NOID, CPT, s_argc, s_argv);
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

void CProbabilityTable::create(int tableSize, double* table)
{
	init();
	m_tableSize = tableSize;
	for (int i = 0; i < m_tableSize; ++i) m_table[i] = table[i];
}

// This function creates the probability table. The order of
// arguments is as follows: child on the last position preceded
// by the parents. Several types of probability tables are supported.
void CProbabilityTable::create(int ptType)
{
	// Safety checks.
	assert(m_argc > 0);
	assert(m_argv != NULL);

	// Allocate the probability table.
	init();

	// Get the problem instance.
	CProblem* prob = m_owner;
	assert(prob != NULL);

	// Enumerate all parent combinations.
	int i, j;
	int child = m_argv[m_argc - 1];
	int k = prob->getStaticDomainSize(child);
	int np = m_argc - 1;	// Number of parents.
	int npc = 1;			// Number of parent value combinations.
	for (i = 0; i < np; ++i) 
		npc *= prob->getStaticDomainSize(m_argv[i]);

	int* parents = (np > 0) ? (new int[np]) : NULL;
	if (np > 0)
	{
		for (i = 0; i < np; ++i) parents[i] = 0;
		parents[np - 1] = -1;
	}

	// Process all parent value combinations.
	for (j = 0 ; j < npc*k ; j += k) 
	{
		// Enumerate parent value combination. 
		// Find next combination from the current.
		for (int pcai = np - 1 ; pcai >= 0 ; --pcai) 
		{
			if (++parents[pcai] < k) break ;
			parents[pcai] = 0 ;
		}

		// Check probability table type	
		switch (ptType)
		{
		case PT_UNIFORM:
			{
				double norm = (double) 0.0;		// Normalization constant.
				double p;						// Probability of child's value.
				int l;
				for (l = 0 ; l < k ; ++l) 
				{
					p = 0.0001 + RandUniformDouble();
					norm += m_table[j + l] = (double) p;
				}

				// Normalize.
				for (l = 0 ; l < k ; l++) m_table[j + l] /= norm;

				break;
			}
		case PT_UNIFORM_BIASED:
			{
				double norm = (double) 0.0;		// Normalization constant.
				double p;						// Probability of child's value.

				// Biased decision.
				double bias = 0.2;
				double r = RandUniformDouble();
				if (r < bias)
				{
					int val = (int)RandUniform(k);
					for (int l = 0; l < k; ++l)
					{
						m_table[j + l] = (l == val) ? 1.0 : 0.0;
					}
				}
				else
				{
					int l;
					for (l = 0 ; l < k ; ++l) 
					{
						p = 0.0001 + RandUniformDouble();
						norm += m_table[j + l] = (double) p;
					}

					// Normalize.
					for (l = 0 ; l < k ; l++) m_table[j + l] /= norm;
				}

				break;
			}
		case PT_NOISYOR:
			{
				// TO BE DECIDED HOW ...
				break;
			}
		}
	}

	if (parents)
		delete[] parents;
}

bool CProbabilityTable::verify()
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
/*
void CProbabilityTable::print(bool scope, bool table)
{
	if (scope && !table)
	{
		cout << "(";
		for (int i = 0; i < m_argc - 1; ++i)
			cout << m_argv[i] << ",";
		cout << m_argv[i] << ")";
	}
}
*/

void CProbabilityTable::print(bool scope, bool table)
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