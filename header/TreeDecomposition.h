// TreeDecomposition.h: interface for the CTreeDecomposition class.
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


#ifndef _MPELIB_TREEDECOMPOSITION_H
#define _MPELIB_TREEDECOMPOSITION_H

#include "defs.h"
#include "TreeDecompositionCluster.h"
#include "TreeDecompositionSeparator.h"
#include "problem.h"
#include "ProblemMixed.h"


////////////////////////////////////////////////////////////////////////
// CProblemMixed class implementation.

class CTreeDecomposition  
{	
	friend class CProblem;
	friend class CTreeDecompositionCluster;
	friend class CTreeDecompositionSeparator;

protected: 
	int m_size;				// number of clusters in the tree decomposition
	int m_width;			// max cluster size of the tree decomposition
	CProblem* m_owner;
	treeDecompositionCluster_v m_clusters;				// Complete cluster structure.
	CTreeDecompositionCluster* m_root;				// m_clusters[m_root] is the root of the tree decomposition
							// need safety checks here m_root < m_clusters.size();
							// and the links in each cluster should be correct (m_parent, m_children)
	treeDecompositionSeparator_v m_separators;			// Separators between clusters

	CLegalTreeNode*			m_cutsetTreeNode;			// used in building the andor cutset tree

protected:
	void destroy();						// Destroy the cluster structure.



public:
	CTreeDecomposition(CProblem *owner);
	virtual ~CTreeDecomposition();

	CProblem* getOwner()				{ return m_owner; };	
	void setOwner(CProblem* owner)		{ m_owner =  owner; };	

	CTreeDecompositionCluster* getRootCluster()		{ return m_root; };
	void setRootCluster(CTreeDecompositionCluster* cluster)		{ m_root = cluster; };

	treeDecompositionCluster_v& getClusters()		{ return m_clusters; };
	void setClusters(treeDecompositionCluster_v &clusters);	

	treeDecompositionSeparator_v& getSeparators()		{ return m_separators; };
	void setSeparators(treeDecompositionSeparator_v &separators);

	CLegalTreeNode* getCutsetTreeNode()				{ return m_cutsetTreeNode; };
	void setCutsetTreeNode(CLegalTreeNode* node)	{ m_cutsetTreeNode = node; }; // points to the node in the pseudo tree where the next variable will be attached


	int getWidth()						{ return m_width; };
	int setWidth(int width)				{ m_width = width; };

	int getSize()						{ return m_size; };
	void setSize(int size)				{ m_size = size; };

	void print();

	void eliminateVar(int var);				// eliminates var from the tree, if it is contained

	int computeSeparatorWeights(int scopeSize);
	// get the weight of subtree rooted at cluster, but which does not include subtree rooted at cutoff_cluster
	// the clusters are thinned by variables (we eliminate the 'variables' from each cluster)
	double computeSubtreeWeight(CTreeDecompositionCluster* root, CTreeDecompositionCluster* cutoff, variable_s variables, int scopeSize);  
	
	int computeRemainingWidth(variable_s& variables);
	int computeRemainingWidth(int var);

	CTreeDecompositionSeparator* highestWeightSeparator();
	int connectedComponents();		// given the variables to be eliminated (of a separator), identify connected components and return their number
	void getComponent(int component, treeDecompositionCluster_v &clusters );

	CTreeDecompositionSeparator* getSeparator(CTreeDecompositionCluster* parent, CTreeDecompositionCluster* child); // get the separator between parent and child

	int computeBestVariableGWC(int w);		// computes the best variable to be eliminated by GWC heuristic
	
	bool cyclic(); // tests treeness

};

#endif // _MPELIB_TREEDECOMPOSITION_H

// Local Variables:
// mode: C++
// End:
