// cpt.h -- Bayesian Conditional Probability Table (aka CPT).

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

#pragma warning (disable : 4786)

#ifndef _MPELIB_CPT_H
#define _MPELIB_CPT_H

#include "function.h"

#define PT_UNIFORM			0
#define PT_UNIFORM_BIASED	1
#define PT_NOISYOR			2


//////////////////////////////////////////////////////////////////////////
// CProbabilityTable class interface.
class CProbabilityTable : public CFunction
{
protected:
	
	double* m_table;		// Probability table.
	int m_tableSize;		// Probability table size;

public:
	CProbabilityTable();
	CProbabilityTable(int id, CTYPE type, int argc, int* argv);
	CProbabilityTable(int id, CTYPE type, int argc, int* argv, int, double*);
	virtual ~CProbabilityTable();

public:
	void init();

	void create(int ptType = PT_UNIFORM);
	void create(int tableSize, double* table);


	double getValueAt(int addr);
	void   setValueAt(int addr, double val);

	int getCurrentValue(void* val);
	int setCurrentValue(void* val);

	void substitute(int var, int val, CFunction*& fun);

	int getChild();
	int getParent(int idx);
	int getNumberOfParents();

	void destroy();
	void print(bool scope, bool table);
	bool verify();

	double* getTable() { return m_table; };
	int getTableSize() { return m_tableSize; };
};


#endif	// _DFSLIB_CPT_H

// Local Variables:
// mode: C++
// End:
