// function.h -- Abstract super class for defining functions.

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

#ifndef _MPELIB_FUNCTION_H
#define _MPELIB_FUNCTION_H

#include "defs.h"

typedef enum {unknown, CPT, CONSTRAINT, CACHE, COST} CTYPE;

class CProblem;

//////////////////////////////////////////////////////////////////////////
// CFunction abstract class interface.
class CFunction
{
protected:
	// The reference counter.
	long m_reference;

	// The constraint's id.
	int m_id;

	// The constraint's type.
	CTYPE m_type;

	// The constraint's scope.
	int* m_argv;
	int  m_argc;

	// Problem instance.
	CProblem* m_owner;

	// Original constraint flag.
	bool m_original;
	bool m_constant;
	bool m_used;

public:
	CFunction();
	CFunction(int id, CTYPE type, int argc, int* argv);

	virtual ~CFunction();

public:
	void _addRef()					{ ++m_reference; };
	long _release()					{ return --m_reference; };

	virtual bool isMemberOf(int var);

	virtual int getID()				{ return m_id; };
	virtual void setID(int id)		{ m_id = id; };

	virtual CTYPE type()			{ return m_type; };
	virtual void setType(CTYPE type) { m_type = type; };

	virtual int* getArgv()			{ return m_argv; };
	virtual int  getArgc()			{ return m_argc; };
	virtual int  getArity();

	virtual int getAddress(int argc, int* vals);
	virtual int getAddress(int argc, int* vars, int* vals);

	bool isOriginal()				{ return m_original; };
	bool isConstant()				{ return m_constant; };

	void setOriginal(bool flg)		{ m_original = flg; };
	void setConstant(bool flg)		{ m_constant = flg; };

	CProblem* getOwner()			{ return m_owner; };
	void setOwner(CProblem* owner)	{ m_owner = owner; };

	bool isUsed()					{ return m_used; };
	void setUsed(bool flag)			{ m_used = flag; };

	virtual int getCurrentValue(void* val) = 0;
	virtual int setCurrentValue(void* val) = 0;

	virtual void substitute(int var, int val, CFunction*& fun) = 0;
	virtual bool verify() = 0;
protected:
	// Initialize an empty constraint, based on the scope.
	virtual void init() = 0;

public:
	virtual void print(bool scope, bool table) = 0;

	int getCurrentAddress();

	int computeTableSize();
	int computeTableSize(int argc, int* argv);

	// Checks if all variables in scope are instantiated.
	bool isScopeBounded();

	// Returns the highest node in the ordering.
	int getHighestNode(int& highestVar, int& highestPos, bool ignoreAssigned = false);
};

//////////////////////////////////////////////////////////////////////////
// Class implementing a list of constraints.
class CFunctionList
{
private:
	CFunctionList(const CFunctionList&);
	CFunctionList& operator=(const CFunctionList&);

protected:
	CProblem* m_owner;
	int m_root;					// Root variable.
	function_v m_functions;
	variable_s m_variables;		// List of variables connected to.

public:
	void add(CFunction* fun);
	int size();

	void setOwner(CProblem* csp)				{ m_owner = csp; };
	CProblem* getOwner()						{ return m_owner; };

	void setRoot(int root)						{ m_root = root; };
	int getRoot()								{ return m_root; };

	variable_s& getConnections()				{ return m_variables; };

public:
	CFunctionList(CProblem* csp = NULL);
	virtual ~CFunctionList();
};


inline bool varity_asc(CFunction* a, CFunction* b)
{
	return (a->getArity() < b->getArity());
}

inline bool varity_desc(CFunction* a, CFunction* b)
{
	return (a->getArity() > b->getArity());
}


#endif	// _DFSLIB_FUNCTION_H

// Local Variables:
// mode: C++
// End:
