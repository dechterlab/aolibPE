// aonode.h -- AND-OR search tree node.

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

#ifndef _MPELIB_AONODE_H
#define _MPELIB_AONODE_H

#include "defs.h"
#include "problem.h"

typedef enum{AND, OR} NTYPE;

////////////////////////////////////////////////////////////////////////
// CAONode class interface.

class CAONode
{
	// Copy constructor not supported.
	CAONode(const CAONode&);
	CAONode& operator=(const CAONode&);

protected:
	NTYPE m_type;						// Type of the node (i.e. AND, OR).
	int m_label;						// Label of the node (variable - OR, value - AND).
	bool m_solved;						// Solvable state of the node.

	bool m_updated;						// The node's g-value has been updated from below.

	double m_gvalue;					// g-value of the node.
	double m_gcountvalue;				// g-countvalue of the node - to count solutions of CSP
	double m_lvalue;					// l-value keeps the bucket-multiply value

	CAONode* m_parent;					// Parent link.
	aonode_v m_children;				// Children links.

	double m_cost;						// Value cost estimate (upper bound) for AND nodes.

	vector<VALUECOST*> m_estimates;

	variable_multimap m_invalidValues;

public:
	void clear();
	void destroy();
	NTYPE type()						{ return m_type; };
	int label()							{ return m_label; };
	double cost()						{ return m_cost; };
	CAONode* parent()					{ return m_parent; };
	aonode_v& children()				{ return m_children; };
	int childrenSize()					{ return (int) m_children.size(); };

	void addChild(CAONode* n);
	void setParent(CAONode* n);
	
	double getG()						{ return m_gvalue; };
	void setG(double gvalue)			{ m_gvalue = gvalue; };
	void setCost(double bound)			{ m_cost = bound; };

	double getGCount()					{ return m_gcountvalue; };
	void setGCount(double gcountvalue)	{ m_gcountvalue = gcountvalue; };

	double getL()						{ return m_lvalue; };
	void setL(double lvalue)			{ m_lvalue = lvalue; };


	void remove(CAONode* node);
	void prune(CAONode* node, double g, vector<CAONode*>& siblings);

	bool isUpdated()					{ return m_updated; };
	void setUpdated(bool flag)			{ m_updated = flag; };

	void buildValueCosts(CProblem* prob, int ibound, double& maxCost);
	vector<VALUECOST*>& getValueCosts()	{ return m_estimates; };

	void computeExactSubtree(CProblem* prob, int ibound, double& sumCost);

	variable_multimap& invalidValues()		{ return m_invalidValues; };
	void addInvalidValue(int var, int val);		// adds an invalid pair (var, val) to list of modifications 
	void cleanInvalidValues(bool** table, int* count, int* assignment);


public:
	CAONode(NTYPE type, int label);
	virtual ~CAONode();
};

#endif	// _MPELIB_AONODE_H

// Local Variables:
// mode: C++
// End:

