// function.cpp -- Abstract super class for defining functions.

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


#include "function.h"
#include "problem.h"

// Create an empty function.
CFunction::CFunction()
{
	m_id = NOID;
	m_type = unknown;
	m_argv = NULL;
	m_argc = 0;
	m_original = true;
	m_constant = false;
	m_used = false;
	m_reference = 0;
}

CFunction::CFunction(int id, CTYPE type, int argc, int* argv)
	: m_id(id), m_type(type), m_argc(argc), m_argv(argv)
{
	m_original = true;
	m_constant = false;
	m_used = false;
	m_reference = 0;
}

// Destructor.
CFunction::~CFunction()
{
	if (m_argv)
		delete[] m_argv;
}

/**
 * Test if a given variable is part of the constraint's scope
 * @param [var] variable to look for
 * @return 'true' if so, 'false' otherwise
 */
bool CFunction::isMemberOf(int var)
{
	if (m_argv)
	{
		assert(m_argc != 0);
		for (int i = 0; i < m_argc; ++i)
		{
			if (var == m_argv[i])
				return true;
		}
	}

	return false;
}

// This function checks if the scope is instantiated.
bool CFunction::isScopeBounded()
{
	// Safety checks.
	assert(m_owner != NULL);

	for (int i = 0; i < m_argc; ++i)
	{
		int var = m_argv[i];
		if (!m_owner->isAssigned(var))
			return false;
	}

	return true;
}

int CFunction::getAddress(int argc, int* vals)
{
	// Safety checks.
	assert(argc == m_argc);
	assert(m_argv != NULL);

	int i = argc - 1;
	int addr = vals[i] ;
	int j = m_owner->getStaticDomainSize(m_argv[i]) ;
	for (--i ; i >= 0 ; i--)
	{
		// here we should add to 'adr' the product of :
		//		#-of-values-of-variable-i * current-value-of-i
		int val = vals[i] ;
		addr += j*val ;
		j *= m_owner->getStaticDomainSize(m_argv[i]) ;
	}

	return addr;
}

int CFunction::getAddress(int argc, int* vars, int* vals)
{
	// Safety checks.

	int i = argc - 1;
	int addr = vals[i] ;
	int j = m_owner->getStaticDomainSize(vars[i]) ;
	for (--i ; i >= 0 ; i--)
	{
		// here we should add to 'adr' the product of :
		//		#-of-values-of-variable-i * current-value-of-i
		int val = vals[i] ;
		addr += j*val ;
		j *= m_owner->getStaticDomainSize(vars[i]) ;
	}

	return addr;
}

// This function computes the table address
// of the current scope assignment. Requires that
// all variables in the constraint's scope have
// already been instantiated.
int CFunction::getCurrentAddress()
{
	if (m_constant)
		return 0;

	// Safety checks.
	assert(m_argc > 0);
	assert(m_argv != NULL);
	assert(m_owner != NULL);

	int i = m_argc - 1;
	int addr = m_owner->getValue(m_argv[i]);
	assert(addr >= 0);

	int j = m_owner->getStaticDomainSize(m_argv[i]);
	for (--i ; i >= 0 ; i--)
	{
		// here we should add to 'adr' the product of :
		//		#-of-values-of-variable-i * current-value-of-i
		int val = m_owner->getValue(m_argv[i]);
		assert(val >= 0);

		addr += j*val ;
		j *= m_owner->getStaticDomainSize(m_argv[i]) ;
	}

	return addr;

}

// This function computes the table size. Requires the scope
// to have been already allocated.
int CFunction::computeTableSize()
{
	assert(m_owner != NULL);
	assert(m_argv != NULL);
	assert(m_argc > 0);

	int size = 1;
	for (int i = 0; i < m_argc; ++i)
		size *= m_owner->getStaticDomainSize(m_argv[i]);

	return size;
}

int CFunction::computeTableSize(int argc, int* argv)
{
	assert(m_owner != NULL);
	assert(argc > 0);

	int size = 1;
	for (int i = 0; i < argc; ++i)
		size *= m_owner->getStaticDomainSize(argv[i]);

	return size;
}

// This function returns the highest argument in the ordering.
// Requires that the problem instance has already been defined.
int CFunction::getHighestNode(int& highestVar, int& highestPos, bool ignoreAssigned)
{
	// Safety checks.
	assert(m_owner != NULL);
	assert(m_argv != NULL);
	assert(m_argc > 0);

	int* position = m_owner->getPosition();
	assert(position != NULL);

	int hPos = -1, hVar = -1;
	for (int i = 0; i < m_argc; ++i)
	{
		int arg = m_argv[i];

		int pos = position[arg];
		if (pos > hPos)
		{
			hPos = pos;
			hVar = arg;
		}
	}

	highestPos = hPos;
	highestVar = hVar;

	return 1;
}

int CFunction::getArity()
{
	// Safety checks.
	assert(m_owner != NULL);

	CProblem* prob = m_owner;
	int arity = 0;
	for (int i = 0; i < m_argc; ++i)
	{
		int arg = m_argv[i];
		if (prob->isAssigned(arg))
			continue;	// skip evidence.

		++arity;
	}

	return arity;
}

//////////////////////////////////////////////////////////////////////
// CFunctionList class implementation

CFunctionList::CFunctionList(CProblem* csp)
{
	m_owner = csp;
	m_root = NOID;
}

CFunctionList::~CFunctionList()
{
	m_functions.clear();
	m_variables.clear();
}

void CFunctionList::add(CFunction* fun)
{
	// Safety checks.
	assert(fun != NULL);

	m_functions.push_back(fun);

	int* argv = fun->getArgv();
	int argc = fun->getArgc();
	for (int i = 0; i < argc; ++i)
	{
		int var = argv[i];
		if (var != m_root)
		{
			if (m_variables.find(var) == m_variables.end())
				m_variables.insert(var);
		}
	}
}

int CFunctionList::size()
{
	return (int) m_functions.size();
}


// Local Variables:
// mode: C++
// End:
