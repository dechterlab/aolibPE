// bucket.h -- Bucket structure.

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

#ifndef _MPELIB_BUCKET_H
#define _MPELIB_BUCKET_H

#include "defs.h"
#include "function.h"

class CBucketStruct;

//////////////////////////////////////////////////////////////////////////
// CBucket class interface.
class CBucket
{
private:
	// Copy constructor not supported.
	CBucket(const CBucket&);
	CBucket& operator=(const CBucket&);

protected:
	int m_index;					// Bucket's index.
	int m_variable;					// Bucket's variable.
	function_v m_functions;			// Functions in the bucket.
	function_v m_augment;			// Used for heuristic pre-computation.
	minibucket_v m_minibuckets;		// List of mini-buckets (partition).

	CBucketStruct* m_owner;			// Bucket structure of the problem.
	CFunction* m_output;			// The output function of the bucket.
	
public:
	// Constructor/Destructor
	CBucket(int var = UNKNOWN);
	virtual ~CBucket();

protected:
	void destroy();					// Destroy the content of the bucket.
	void createPartition(int ib);	// Create mini-buckets partition.
	void destroyPartition();		// Destroy mini-buckets partition.

public:
	void clean();					// Clean the bucket of non-original functions.

public:
	int size();						// Bucket's size.
	void print();
	function_v& functions()			{ return m_functions; };

	int getVariable()				{ return m_variable; };
	void setVariable(int var)		{ m_variable = var; };

	void add(CFunction* c);			// Add a new funcition.
	void remove(CFunction* c);		// Remove a function.

	void setOwner(CBucketStruct* owner)	{ m_owner = owner; };
	CBucketStruct* getOwner()			{ return m_owner; };

	CFunction* getOutput()			{ return m_output; };

	void setIndex(int idx)			{ m_index = idx; };
	int  getIndex()					{ return m_index; };

public:
	void augment(CFunction* c);							// Augment the bucket.
	bool process(int type, int ibound, int bottom);		// Process the bucket.
	bool solve(int k, double*& bounds);					// Solve first bucket in the ordering.

	double multiply();									// Multiply all functions in the bucket.

	bool evalStaticHeuristic(int k, double*& costs);
	bool evalDynamicHeuristic(int ibound, int k, double*& costs);
};

#endif	// _CSPLIBR_BUCKET_H

// Local Variables:
// mode: C++
// End:
