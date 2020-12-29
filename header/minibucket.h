
// minibucket.h -- Mini-Bucket structure.

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

#pragma warning (disable : 4786)

#ifndef _MPELIB_MINIBUCKET_H
#define _MPELIB_MINIBUCKET_H

#include "defs.h"
#include "function.h"
#include "bucket.h"

/////////////////////////////////////////////////////////////////////////
// "Mini-Bucket" structure.
class CMiniBucket
{
	friend class CBucket;

private:
	// Copy constructor not supported.
	CMiniBucket(const CMiniBucket&);
	CMiniBucket& operator=(const CMiniBucket&);

protected:
	int m_id;						// Mini-Bucket's ID.
	CBucket* m_owner;				// Mini-Bucket's owner (it is a CBucket).
	function_v m_functions;			// Mini-Bucket's functions.
	variable_s m_variables;			// Mini-Bucket's variables (ordered set).

public:
	int getID()						{ return m_id; };
	void setID(int id)				{ m_id = id; };

	CBucket* getOwner()				{ return m_owner; };
	void setOwner(CBucket* owner)	{ m_owner = owner; };

	int size();						// Mini-Bucket's size (as number of functions).
	bool isEmpty();					// Is the mini-bucket empty?

	bool isScopeInBound(int ibound, int argc, int* argv);
	
	bool process(CFunction*& fun);

	void add(CFunction* c);		// Add a new function.
	void remove(CFunction* c);	// Remove a function.

protected:

	// Processing methods.
	bool prodMax(int elim, CProblem* csp, CFunction*& fun);
	bool prodSum(int elim, CProblem* prob, CFunction*& fun);

public:
	CMiniBucket(int id = NOID, CBucket* owner = NULL);
	virtual ~CMiniBucket();
};

#endif	// _CSPLIBR_MINIBUCKET_H

// Local Variables:
// mode: C++
// End:
