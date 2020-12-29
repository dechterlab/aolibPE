// ConstraintTable.h: Constraint Tables.
//
//////////////////////////////////////////////////////////////////////

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

/*
 * NOTE: This is an internal header file.
 * You should not attempt to use it directly.
 */

#pragma warning (disable : 4786)

#ifndef _MPELIB_CONSTRAINT_H
#define _MPELIB_CONSTRAINT_H

#include "function.h"



//////////////////////////////////////////////////////////////////////////
// CConstraintTable class interface.
class CConstraintTable : public CFunction  
{
protected:
	
	int* m_table;		// Constraint table.
	int m_tableSize;		// Constraint table size;
	int m_highestVarInScope; // Highest var in scope

	// implement a Shannon tree to represent the tuples; this is to be used for constraint propagation (unit resolution)
	// a somewhat better idea would be a BDD

	CConstraintNode* m_rootShannonTree; // root of the Shannon tree of the constraint
	int* m_argvShannonTree; // the scope ordered according to problem ordering, for the Shannon Tree


public:
	CConstraintTable();
	CConstraintTable(int id, CTYPE type, int argc, int* argv);

	virtual ~CConstraintTable();

public:
	void init();

	void create(double tightness); //tightness = percentage of allowed tuples
	void create(int tableSize, int* table);

	int get_highestVarInScope()			{ return m_highestVarInScope; }; 
	void set_highestVarInScope(int val) { m_highestVarInScope = val; };
	
	int getValueAt(int addr);
	void setValueAt(int addr, int val);

	int getCurrentValue(void* val);
	int setCurrentValue(void* val);

	void substitute(int var, int val, CFunction*& fun);

//	int getChild();
	int getParent(int idx);
	int getNumberOfParents();

	void destroy();
	void print(bool scope, bool table);
	bool verify();

	int* getTable() { return m_table; };
	int getTableSize() { return m_tableSize; };

	CConstraintNode* getRootShannonTree()			{ return m_rootShannonTree; };
	int* getArgvShannonTree()						{ return m_argvShannonTree; };
	void setRootShannonTree(CConstraintNode* root)	{ m_rootShannonTree = root; };
	void createShannonTree(); 
	void destroyShannonTree();
	void printShannonTree();
};

#endif // _MPELIB_CONSTRAINT_H

// Local Variables:
// mode: C++
// End:
