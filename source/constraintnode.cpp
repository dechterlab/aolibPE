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

#include "constraintnode.h"


////////////////////////////////////////////////////////////////////////
// CConstraintNode class implementation.

CConstraintNode::CConstraintNode()
{
	m_parent = NULL;
	m_gvalue = 0;
	m_children = NULL;
}
CConstraintNode::CConstraintNode(int var, int val)
	: m_variable(var), m_value(val)
{
	m_parent = NULL;
	m_gvalue = 0;
	m_children = NULL;
}

CConstraintNode::~CConstraintNode()
{
	destroy();
}

void CConstraintNode::destroy()
{
	m_parent = NULL;
	if(m_children)
		delete[] m_children;
}

void CConstraintNode::createChildren(int var, int size)
{
	assert(size > 0);
	m_children = new CConstraintNode*[size];
	for(int i=0; i < size; i++)
	{
		m_children[i] = new CConstraintNode(var, i);		
		m_children[i]->setParent(this);
	}
}


void CConstraintNode::setParent(CConstraintNode* n)
{
	m_parent = n;
}


// Local Variables:
// mode: C++
// End:
