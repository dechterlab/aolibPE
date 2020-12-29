// bucketstruct.h -- Bucket structure associated with a (sub)problem.

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

#ifndef _MPELIB_BUCKETSTRUCT_H
#define _MPELIB_BUCKETSTRUCT_H

#include "bucket.h"
#include "minibucket.h"
#include "problem.h"

class CBucketStruct
{
	friend class CProblem;
	friend class CBucket;
	friend class CMiniBucket;

private:
	// Copy constructor not supported.
	CBucketStruct(const CBucketStruct&);
	CBucketStruct& operator=(const CBucketStruct&);

protected:
	CProblem* m_owner;					// Problem instance.
	CBucket** m_buckets;				// Complete bucket structure.
	int m_size;							// Size of the structure.

	int* m_ordering;					// Local variable ordering.
	int* m_position;					// Local variable position.
	bool m_locals;						// Flag indicating local ordering/position allocated.

protected:
	void destroy();						// Destroy the bucket structure.

public:	
	void init(bool dfsOrdering);		// Initialize the bucket structure.
	void init(int n, variable_v& ordering);
	void clean();						// Clean structure (preserve originals only).

	void print();

	CBucket* getBucketAt(int pos);		// Lookup a bucket at specified position.
	void getHighestNode(CFunction* fun, int& highestVar, int& highestPos);

	CProblem* getOwner()				{ return m_owner; };

	bool process(int type, int ibound, int& k, double*& result, int level = 0);
public:
	CBucketStruct(CProblem* owner);
	virtual ~CBucketStruct();
};

#endif	// _MPELIB_BUCKETSTRUCT_H

// Local Variables:
// mode: C++
// End:
