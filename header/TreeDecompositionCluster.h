// TreeDecompositionNode.h: interface for the CTreeDecompositionNode class.
//
//////////////////////////////////////////////////////////////////////

/*
 * Copyright (C) 2005 Robert Mateescu
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


#ifndef _MPELIB_TREEDECOMPOSITIONCLUSTER_H
#define _MPELIB_TREEDECOMPOSITIONCLUSTER_H


#include "defs.h"
#include "problem.h"

class CTreeDecomposition;

////////////////////////////////////////////////////////////////////////
// CProblemMixed class implementation.
// defines a cluster of a tree decomposition

class CTreeDecompositionCluster  
{

protected:
	int		m_id;						// the id of the cluster
	variable_s	m_variables;			// ascending list of the variables contained in the cluster
	CTreeDecomposition* m_owner;		// owner of the cluster (tree decomposition)
	CProblem* m_owner_prob;
	CTreeDecompositionCluster* m_parent;	// parent in the tree decomposition
	treeDecompositionCluster_v m_children;	// children in the tree decomposition
	int m_visited;						// used in traversing the tree decomposition to assemble the cutset

public:
	CTreeDecompositionCluster();
	CTreeDecompositionCluster(int id, variable_s& variables);
	virtual ~CTreeDecompositionCluster();
	void destroy();


	void setOwner(CTreeDecomposition* owner)	{ m_owner = owner; };
	CTreeDecomposition* getOwner()			{ return m_owner; };

	void setOwnerProb(CProblem* owner)	{ m_owner_prob = owner; };
	CProblem* getOwnerProb()			{ return m_owner_prob; };

	int getVisited()					{ return m_visited; };
	void setVisited(int visited)		{ m_visited = visited; };

	bool isMemberOf(int var);	// checks if var belongs to cluster
	variable_s& getVariables()			{ return m_variables; };
	int getSize()			{ return (int) m_variables.size(); };
	int getID()				{ return m_id; };
	void setID(int id)		{ m_id = id; };
	void print();
	void addVar(int var);									// adds one more variable to the scope
	void deleteVar(int var);									// deletes variable var from the scope

	CTreeDecompositionCluster* parent()					{ return m_parent; };
	treeDecompositionCluster_v& children()				{ return m_children; };
	int childrenSize()					{ return (int) m_children.size(); };
	void addChild(CTreeDecompositionCluster* cluster);
	void setParent(CTreeDecompositionCluster* cluster)	{ m_parent = cluster; };
	int highestVar();

	void eliminateVar(int var);		// eliminates var from cluster, if it is contained in it
	

};

#endif // _MPELIB_TREEDECOMPOSITIONCLUSTER_H


// Local Variables:
// mode: C++
// End:
