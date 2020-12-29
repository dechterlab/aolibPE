// legaltree.h -- Rooted Tree Arangement.

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

#ifndef _MPELIB_LEGALTREE_H
#define _MPELIB_LEGALTREE_H

#include "defs.h"
#include "graph.h"

class CLegalTreeNode;

//////////////////////////////////////////////////////////////////////////
// CLegalTree class interface.
class CLegalTree
{
	friend class CLegalTreeNode;

protected:
	CLegalTreeNode* m_root;				// Root of the tree.
	CLegalTreeNode** m_nodes;			// Nodes of the tree.
	int m_height;						// Height of the tree.
	int* m_height_w_cutset;				// Vector of heights of AO-w-cutsets; m_height_w_cutset[0] = m_height;
	int m_size;							// Size of the tree as number of nodes.
	int* m_ordering;
	int* m_position;

public:
	void destroy();						// Destroy the tree.
	int height() { return m_height; };	// Return the height of the tree.
	int size() { return m_size; };		// Return the size of the tree.

	void create(CGraph* graph);			// Create rooted tree arrangement (induced graph).
	void createDFS(CGraph* graph);		// Create DFS tree of the moral graph.
	void createChain(int ordering[]);
	
	void eraseLinks();					// This prepares the legal tree to be made into a chain

	CLegalTreeNode** getNodes()			{ return m_nodes; };
	void getChildren(int var, variable_v& children);

	int* getWCutsetHeight()		{ return m_height_w_cutset; };

	void order(int*& ordering, int*& position);		// DFS ordering of the legal tree.
	int descendants(int node, variable_v& desc);

	void setRoot (CLegalTreeNode* root)		{ m_root = root; };
	CLegalTreeNode* getRoot()				{ return m_root; };

	// add a chain of variables to node; used for w-cutset pseudo-tree; variables is and ordered vector, 
	//.begin() renders smallest width of remaining graph, etc.
	// w is an array inidcating the width of the remaining graph after assigning the corresponding variable
	CLegalTreeNode* addVariables(CLegalTreeNode* parent, variable_v &variables, int *&w); 

	virtual int* getOrdering();
	virtual int* getPosition();

	void setOrdering(int* ordering);
	void setPosition(int* position);


	void print();
public:
	CLegalTree(int size, int* ordering, int* position);
	virtual ~CLegalTree();
};

//////////////////////////////////////////////////////////////////////////
// CLegalTreeNode class interface.
class CLegalTreeNode
{
	friend class CLegalTree;

protected:
	CLegalTree* m_owner;				// Tree instance.

	CLegalTreeNode* m_parent;			// Parent of node.
	legaltreenode_v m_children;			// Children on node.

	int m_variable;						// Variable assigned to node.
	int m_height;						// Relative height assigned to node.
	int m_w_cutset;						// The max w for which (node belongs to minimal w-cutset)


public:
	void destroy();						// Destroy the node.
	legaltreenode_v& children()			{ return m_children; };
	CLegalTreeNode* parent()			{ return m_parent; };
	int variable()						{ return m_variable; };
	void setHeight(int height)			{ m_height = height; };
	void setWCutset(int w)				{ m_w_cutset = w; };
	void setParent(CLegalTreeNode* parent)	{ m_parent = parent; };
	int getHeight()						{ return m_height; };
	int getWCutsetHeight()				{ return m_w_cutset; };

public:
	CLegalTreeNode(int var);
	virtual ~CLegalTreeNode();
};



#endif	// _DFSLIB_LEGALTREE_H

// Local Variables:
// mode: C++
// End:
