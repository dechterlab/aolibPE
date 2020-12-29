// constraintnode.cpp -- Node to represent the constraint table as a Shannon Tree

/*
 * Copyright (C) 2006 Robert Mateescu
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

#ifndef _MPELIB_CONSTRAINTNODE_H
#define _MPELIB_CONSTRAINTNODE_H

#include "defs.h"
#include "problem.h"
#include "constraintnode.h"

////////////////////////////////////////////////////////////////////////
// CConstraintNode class interface.

class CConstraintNode
{
	// Copy constructor not supported.
	CConstraintNode(const CConstraintNode&);
	CConstraintNode& operator=(const CConstraintNode&);

protected:
	int m_variable;						// Variable of the Node
	int m_value;						// Value of the variable

	double m_gvalue;					// g-value of the node.

	CConstraintNode* m_parent;					// Parent link.
	CConstraintNode** m_children;				// Children links.

public:
	void destroy();
	void setVariable(int var)			{ m_variable = var; };
	int variable()						{ return m_variable; };
	void setValue(int val)				{ m_value = val; };
	int value()							{ return m_value; };

	CConstraintNode* parent()			{ return m_parent; };
	CConstraintNode** children()		{ return m_children; };

	void createChildren(int var, int size);
	void setParent(CConstraintNode* n);
	
	double getG()						{ return m_gvalue; };
	void setG(double gvalue)			{ m_gvalue = gvalue; };
	void addToG(double val)				{ m_gvalue += val; }
	

public:
	CConstraintNode::CConstraintNode();
	CConstraintNode(int var, int val);
	virtual ~CConstraintNode();
};


#endif	// _MPELIB_CONSTRAINTNODE_H

// Local Variables:
// mode: C++
// End:

