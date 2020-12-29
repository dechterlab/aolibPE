// legaltree.cpp -- Rooted Tree Arangement.

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


#include "legaltree.h"

//////////////////////////////////////////////////////////////////////////
// CLegalTree class implementation.
CLegalTree::CLegalTree(int size, int* ordering, int* position)
{
	m_height = 0;
	m_size = size;

	m_nodes = new CLegalTreeNode*[size];
	for (int i = 0; i < size; ++i)
	{
		int var = ordering[i];
		m_nodes[i] = new CLegalTreeNode(var);
	}

	m_ordering = new int[size];
	memcpy(m_ordering, ordering, size * sizeof(int));

	m_position = new int[size];
	memcpy(m_position, position, size * sizeof(int));

	m_height_w_cutset = new int[size];
	for(int i=0; i<size; i++)
		m_height_w_cutset[i] = 0;	

	m_root = m_nodes[0];
}

CLegalTree::~CLegalTree()
{
	destroy();
}

void CLegalTree::destroy()
{
	if (m_ordering)
		delete[] m_ordering;

	if (m_position)
		delete[] m_position;

	if (m_height_w_cutset)
		delete (m_height_w_cutset);

	if (m_nodes)
	{
		for (int i = 0; i < m_size; ++i)
			delete m_nodes[i];
		delete[] m_nodes;
	}

	m_root = NULL;
}

// This function returns the list of child nodes of a given node.
// The nodes in the tree are arranged according to the variable
// ordering of the induced graph.
void CLegalTree::getChildren(int var, variable_v& children)
{
	// Safety checks.
	assert(m_position != NULL);
	assert(m_ordering != NULL);
	assert((var >= 0) && (var < m_size));

	// Get the tree node associated with the variable.
	CLegalTreeNode* node = m_nodes[m_position[var]];
	assert(node != NULL);

	// Get the list of children.
	children.clear();
	legaltreenode_v& ch = node->children();
	legaltreenode_v::iterator it = ch.begin();
	for (; it != ch.end(); ++it)
	{
		CLegalTreeNode* child = (*it);
		children.push_back(child->variable());
	}
}

// This function creates the rooted tree arrangement from the induced graph.
void CLegalTree::create(CGraph* graph)
{
	// Safety checks.
	assert(m_nodes != NULL);
	assert(m_root != NULL);
	assert(m_size > 0);

	// Visited nodes.
	bool* visited = new bool[m_size];
	for(int i=0; i<m_size; i++)
		visited[i] = false;	

	// Mark root as visited.
	visited[m_root->m_variable] = true;
	stack<int> dfsStack;

	dfsStack.push(m_root->m_variable);
	while (!dfsStack.empty())
	{
		int var = dfsStack.top();

		// Look for the earliest successor.
		int succ = -1;
		int start = m_position[var];
		for (int pos = start + 1; pos < m_size; ++pos)
		{
			int cand = m_ordering[pos];
			if (graph->isConnected(var, cand))
			{
				if (!visited[cand])
				{
					succ = cand;
					break;
				}
			}
		}

		if (-1 == succ)
		{
			dfsStack.pop();
		}
		else
		{
			CLegalTreeNode* currNode = m_nodes[m_position[var]];
			CLegalTreeNode* succNode = m_nodes[m_position[succ]];

			currNode->m_children.push_back(succNode);
			succNode->m_parent = currNode;

			succNode->m_height = currNode->m_height + 1;

			if (succNode->m_height > m_height)
				m_height = succNode->m_height;

			visited[succ] = true;
			dfsStack.push(succ);
		}
	}

	delete[] visited;
}

// This function creates a legal tree that is a chain, such that 
// the AND/OR search space would be identical to the OR space.
void CLegalTree::createChain(int ordering[])
{
	// Safety checks.
	assert(m_nodes != NULL);
	assert(m_root != NULL);
	assert(m_size > 0);

	CLegalTreeNode* currNode = NULL;
	CLegalTreeNode* nextNode = NULL;

	for (int pos = 0; pos < m_size ; ++pos)
	{
		CLegalTreeNode* node = m_nodes[pos];
		node->m_children.clear();
		node->m_parent = NULL;
	}
	
	for (int i = 0; i<m_size-1; i++)
	{
		int currVar = ordering[i];
		int nextVar = ordering[i+1];
		//find node whose var is ordering[i]
		int j;
		for (j = 0; j < m_size; j++)
			if (m_nodes[j]->variable() == currVar) 
			{
				currNode = m_nodes[j];
				break;
			}
					
		for (j = 0; j < m_size; j++)
			if (m_nodes[j]->variable() == nextVar) 
			{
				nextNode = m_nodes[j];
				break;
			}

		currNode->m_children.push_back(nextNode);
		nextNode->m_parent = currNode;

		nextNode->m_height = currNode->m_height + 1;

		if (nextNode->m_height > m_height)
			m_height = nextNode->m_height;
	}
}

void CLegalTree::eraseLinks()
{
	// Safety checks.
	assert(m_nodes != NULL);
	assert(m_root != NULL);
	assert(m_size > 0);

	for (int pos = 0; pos < m_size ; ++pos)
	{
		CLegalTreeNode* currNode = m_nodes[pos];
		
		currNode->m_children.clear();
		currNode->m_parent = NULL;
	}

	m_height = 0;
}

// This function derives a DFS ordering of the legal tree.
void CLegalTree::order(int*& ordering, int*& position)
{
	// Safety checks.
	assert(m_nodes != NULL);
	assert(m_root != NULL);
	assert(m_size > 0);

	// Allocate DFS ordering.
	if (ordering) delete[] ordering;
	if (position) delete[] position;

	ordering = new int[m_size];
	assert(ordering != NULL);
	for(int i=0; i<m_size; i++)
		ordering[i] = -1;	

	position = new int[m_size];
	assert(position != NULL);
	for(int i=0; i<m_size; i++)
		position[i] = -1;	

	// DFS ordering of the legal tree.
	variable_v succ;
	stack<int> dfsStack;
	int root = m_root->variable();
	dfsStack.push(root);
	int pos = 0;
	
	while (!dfsStack.empty())
	{
		int node = dfsStack.top();
		dfsStack.pop();

		// Add current node to ordering.
		ordering[pos] = node;
		position[node] = pos;
		++pos;

		// Put children of node onto stack (in reverse order).
		getChildren(node, succ);

		variable_v::reverse_iterator it = succ.rbegin();
		for (; it != succ.rend(); ++it)
		{
			int child = (*it);
			dfsStack.push(child);
		}

		succ.clear();
	}
}

// This function creates a descendant list of a specific node in 
// the legal tree. The order of the descendants is a DFS ordering.
int CLegalTree::descendants(int node, variable_v& desc)
{
	// Safety checks.
	assert(m_nodes != NULL);
	assert(m_root != NULL);
	assert(m_size > 0);

	desc.clear();

	// DFS ordering of the legal sub-tree.
	variable_v succ;
	stack<int> dfsStack;
	int root = node;
	dfsStack.push(root);
	int count = 0;
	
	while (!dfsStack.empty())
	{
		int n = dfsStack.top();
		dfsStack.pop();

		// Add current node to ordering.
		desc.push_back(n);
		++count;

		// Put children of node onto stack (in reverse order).
		getChildren(n, succ);

		variable_v::reverse_iterator it = succ.rbegin();
		for (; it != succ.rend(); ++it)
		{
			int child = (*it);
			dfsStack.push(child);
		}

		succ.clear();
	}

	return count;
}

void CLegalTree::print()
{
	cout << "\n Legal Tree rooted at: " << m_root->m_variable;
	for (int i = 0; i < m_size; ++i)
	{
		CLegalTreeNode* node = m_nodes[i];
		cout << "\n  - node: " << node->m_variable;

		if (node->m_parent)
			cout << "[" << node->m_parent->m_variable << "] - ";
		else
			cout << "[] - ";

		legaltreenode_v::iterator it = node->m_children.begin();
		for (; it != node->m_children.end(); ++it)
		{
			CLegalTreeNode* child = (*it);
			cout << child->m_variable << " ";
		}
		//if (i%150 == 149)
		//	break;
	}
	cout << endl;

	cout << "tree height:" << m_height << endl;
}

// This function creates the DFS tree of the moral graph.
void CLegalTree::createDFS(CGraph* graph)
{
	// Safety checks.
	assert(m_nodes != NULL);
	assert(m_root != NULL);
	assert(m_size > 0);

	// Visited nodes.
	bool* visited = new bool[m_size];
	for(int i=0; i<m_size; i++)
		visited[i] = false;	

	// Mark root as visited.
	stack<int> dfsStack;

	CLegalTreeNode* lastNode = NULL;
	dfsStack.push(m_root->m_variable);
	while (!dfsStack.empty())
	{
		int var = dfsStack.top();
		dfsStack.pop();

		// Get current tree node.
		CLegalTreeNode* currNode = m_nodes[m_position[var]];

		if (!visited[var])
		{
			visited[var] = true;

			// Get neighbors in moral graph.
			variable_s neighbors;
			graph->getMoralNeighbors(var, neighbors);

			int cnt = 0;
			variable_s::iterator it = neighbors.begin();
			for (; it != neighbors.end(); ++it)
			{
				int n = (*it);
				if (!visited[n])
				{
					dfsStack.push(n);
					++cnt;
				}
			}

			neighbors.clear();

			// Make connections to parent variable.
			if (lastNode != NULL)
			{
				lastNode->m_children.push_back(currNode);
				currNode->m_parent = lastNode;
				currNode->m_height = lastNode->m_height + 1;

				if (currNode->m_height > m_height)
					m_height = currNode->m_height;

				// Update last node.
				lastNode = (cnt > 0) ? currNode : currNode->m_parent;
			}
			else
			{
				lastNode = m_root;
			}
		}
	}

	delete[] visited;
}


CLegalTreeNode* CLegalTree::addVariables(CLegalTreeNode* parent, variable_v &variables, int *&w)
{
	CLegalTreeNode* parentNode = parent;
	
	variable_v::iterator it = variables.begin();
	int posCount = 0;

	for ( ; it != variables.end() ; it++)
	{
		// attention: m_ordering and m_position of the legal tree are the old ones, from min-fill;
		// but it's fine, because it's inside the legal tree, and they are consistent.
		// retrieve the node  from m_nodes
		CLegalTreeNode* currNode = m_nodes[m_position[(*it)]];
		
		currNode->m_variable = (*it);
		currNode->m_parent = parentNode;
		parentNode->m_children.push_back(currNode);
				
		currNode->m_height = parentNode->m_height + 1;
		if ( currNode->m_height > m_height )
			m_height = currNode->m_height;

		currNode->m_w_cutset = w[posCount];
		if ( m_height_w_cutset[ w[posCount] ] < currNode->m_height )
			m_height_w_cutset[ w[posCount] ] = currNode->m_height;

		parentNode = currNode;
		posCount++;
	}

	return parentNode;	
}

int* CLegalTree::getOrdering()
{
	return (m_ordering);
}

int* CLegalTree::getPosition()
{
	return (m_position);
}

void CLegalTree::setOrdering(int* ordering)
{
	memcpy(m_ordering, ordering, sizeof(int) * m_size);
}

void CLegalTree::setPosition(int* position)
{
	memcpy(m_position, position, sizeof(int) * m_size);
}


//////////////////////////////////////////////////////////////////////////
// CLegalTreeNode class implementation.
CLegalTreeNode::CLegalTreeNode(int var)
{
	m_variable = var;
	m_height = 0;
	m_w_cutset = 0;
	m_parent = NULL;
}

CLegalTreeNode::~CLegalTreeNode()
{
	destroy();
}

void CLegalTreeNode::destroy()
{
	m_children.clear();
	m_parent = NULL;
}


// Local Variables:
// mode: C++
// End:
