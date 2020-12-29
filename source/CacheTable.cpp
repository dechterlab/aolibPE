// CacheTable.cpp: Cache Table

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


#include "CacheTable.h"
//#include "defs.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CCacheTable::CCacheTable() : CFunction()
{
}

CCacheTable::CCacheTable(int id, CTYPE type, int argc, int* argv) 
	: CFunction(id, type, argc, argv)
{
	m_tableSize = 0;
	m_table = NULL;
}

CCacheTable::~CCacheTable()
{
	if (m_table)
		delete[] m_table;
}

// This function allocates/initializes an empty soft
// constraint. Requires the scope to have been already 
// initialized. Initial value in cache is -1
void CCacheTable::init()
{
	if (!m_constant)
	{
		assert(m_argv != NULL);
		assert(m_argc > 0);

		m_tableSize = computeTableSize();
		m_table = new double[m_tableSize];
		for(int i=0; i<m_tableSize; i++)
			m_table[i] = -1;	
	}
	else
	{
		m_tableSize = 1;
		m_table = new double[m_tableSize];
		for(int i=0; i<m_tableSize; i++)
			m_table[i] = -1;	
	}
}

void CCacheTable::purge()
{
	for(int i=0; i<m_tableSize; i++)
		m_table[i] = -1;	
}

double CCacheTable::getValueAt(int addr)
{
	// Safety checks.	
	assert(m_table != NULL);
	assert(addr < m_tableSize);

	return m_table[addr];
}

void CCacheTable::setValueAt(int addr, double val)
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
int CCacheTable::getCurrentValue(void* val)
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
int CCacheTable::setCurrentValue(void* val)
{
	assert(m_table != NULL);
	
	int addr = (m_constant) ? 0 : getCurrentAddress();
	double v = *((double*) val);
		
	m_table[addr] = v;

	return addr;
}

void CCacheTable::destroy()
{
	if (m_table)
		delete[] m_table;

	m_table = NULL;
	m_tableSize = 0;
	m_constant = false;
}


void CCacheTable::print(bool scope, bool table)
{
	if (scope && !table)
	{
		cout << "(";
		int i;
		for (i = 0; i < m_argc - 1; ++i)
			cout << m_argv[i] << ",";
		cout << m_argv[i] << ")";
	}
}

////////////// This function substitutes a variable in the scope with
// specified value. The result is a new (non-original) function.
void CCacheTable::substitute(int var, int val, CFunction*& fun)
{

}



bool CCacheTable::verify()
{
	return true;
}

