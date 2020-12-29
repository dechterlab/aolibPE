// TreeDecompositionSeparator.h: interface for the CTreeDecompositionSeparator class.
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


#ifndef _MPELIB_TREEDECOMPOSITIONSEPARATOR_H
#define _MPELIB_TREEDECOMPOSITIONSEPARATOR_H


#include "defs.h"

class CTreeDecomposition;


class CTreeDecompositionSeparator  
{
protected:
	variable_s	m_variables;			// ascending list of the variables contained in the separator
	CTreeDecompositionCluster* m_clusterParent;	// the cluster up towards the root
	CTreeDecompositionCluster* m_clusterChild;	// the cluster down towards the root
	CTreeDecomposition* m_owner;		// owner of the cluster (tree decomposition)

	int m_clusterIDParent;					// same as above, only ID of cluster
	int m_clusterIDChild;
	double m_weight;						//weight of the separator
	//double m_weightParent;					// weight up = (width of upper portion of tree-dec) * log(number of clusters up)
	//double m_weightChild;					// similar, down
										// these help find the most balanced separator to be broken for a balanced w-cutset
public:
	CTreeDecompositionSeparator();
	CTreeDecompositionSeparator(variable_s& variables);

	virtual ~CTreeDecompositionSeparator();
	void destroy();

	void setOwner(CTreeDecomposition* owner)				{ m_owner = owner; };
	CTreeDecomposition* getOwner()							{ return m_owner; };


	bool isMemberOf(int var);			// checks if var belongs to cluster
	variable_s& getVariables()								{ return m_variables; };
	int getSize()											{ return (int) m_variables.size(); };
	void print();
	void addVar(int var);				// adds one more variable to the scope
	void deleteVar(int var);			// deletes variable var from the scope

	void setClusterIDParent(int var)						{ m_clusterIDParent = var; };
	void setClusterIDChild(int var)							{ m_clusterIDChild = var; };
	int getClusterIDParent()								{ return m_clusterIDParent; };
	int getClusterIDChild()									{ return m_clusterIDChild; };

	void setClusterParent(CTreeDecompositionCluster* cluster)	{ m_clusterParent = cluster; };
	void setClusterChild(CTreeDecompositionCluster* cluster)	{ m_clusterChild = cluster; };
	CTreeDecompositionCluster* getClusterParent()				{ return m_clusterParent; };
	CTreeDecompositionCluster* getClusterChild()				{ return m_clusterChild; };

	double computeWeight(int scopeSize);
	double weight()												{ return m_weight; };

	void eliminateVar(int var);				// eliminates var from the separator, if it is contained


};


#endif // _MPELIB_TREEDECOMPOSITIONSEPARATOR_H

// Local Variables:
// mode: C++
// End:

