// graph.h -- A graph structure.

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


#ifndef _MPELIB_GRAPH_H
#define _MPELIB_GRAPH_H

#include "defs.h"
#include "function.h"


typedef enum{mindegree, minwidth, minfill, maxcard} OTYPE;

////////////////////////////////////////////////////////////////////////
// CGraph class interface.

class CGraph
{
protected:
	int m_size;							// Size of the graph (nodes).
	char** m_graph;						// The actual graph (adjacency matrix).
	char** m_backup;
	int m_width;						// Width of the induced graph.
	bool m_induced;

public:
	void init(function_map& functions);
	// This function creates the mixed network graph.
	void init(function_map& functions, function_map& constraints);
	void backup();
	void restore();

	void connect(variable_s& nodes);
	void connect(int var1, int var2);
	void disconnect(int var1, int var2);

	virtual void order(int type, int*& ordering, int*& position);
	void createInduced(int* ordering, int* position);
	CGraph* clone();
	void extractCliques(int* ordering, vector<set<int> >& cliques);

	bool isConnected();
	bool isConnected(int var1, int var2);

	int getMoralNeighbors(int var, variable_s& neighbors);

	int getWidth()						{ return m_width; };
	void setWidth(int w)				{ m_width = w;};

	bool isEdge(int var1, int var2, bool original = false);

	// was protected, made public by Robert
	int width(int* ordering, int* position);
	void cleanConnections();
	int neighborhood(int var, variable_s& neighbors);

protected:
	int fillset(int var, bool* used);

protected:
	int select(OTYPE type, bool* used);

private:
	void orderingDefault(int*& ordering, int*& position);
	void orderingMinDegree(int*& ordering, int*& position);
	void orderingMinWidth(int*& ordering, int*& position);
	void orderingMinFill(int*& ordering, int*& position);
	void orderingMaxCard(int*& ordering, int*& position);
	void orderingRandomized(int*& ordering, int*& position);

public:
	CGraph(int n);
	virtual ~CGraph();
};


////////////////////////////////////////////////////////////////////////
// from Radu
class CGraphNode
{
public:
	int m_variable;

public:
	CGraphNode(int v) : m_variable(v) {};
	~CGraphNode() {};
};


////////////////////////////////////////////////////////////////////////
// CGraphHash class interface.
class CGraphHash
{
protected:
	set<int> m_vertices;
	map<int,variable_s*> m_neighbors;

public:
	variable_s* neighbors(int var);

public:
	void addNode(int i);
	void addEdge(int i, int j);
	void addDirectedEdge(int i, int j);
	void addClique(vector<int>& clique);
	void addClique(int n, int* clique);
	void addCliqueSet(set<int>& clique);
	
	bool containsEdge(int i, int j);

	void removeNode(int i);
	void removeEdge(int i, int j);

	int degree(int i);

	void removeAndConnect(int i);
	void removeAndConnectAncestral(int i, variable_s& parents); // add only ancestral edges

	void induce(vector<int>& order, int& width);
	void eliminate(int etype, vector<int>& order, int* domains, bool* deterministic, variable_s* parents);
	void eliminate(int root, int etype, vector<int>& order, int* domains);

	double scoreFunction(int etype, int var, int* domains);

	void clear();

	bool isConnected();
	void init(function_v& funs, map<int,int>& idxVarToKey, set<int>& ancestors);
	void addAncestralEdges(variable_s* parents);

protected:
	double scoreMinDegree(int var);
	double scoreMinWidth(int var);
	double scoreMinFill(int var);
	double scoreMinWeight(int var, int* domains);
	double scoreMaxDomainMinFill(int var, int* domains);
	double scoreWeightedMinFill(int var, int* domains);

	

public:
	CGraphHash(CGraphHash& g);
	CGraphHash& operator=(const CGraphHash& g);

	CGraphHash();
	~CGraphHash();
};



#endif	// _DFSLIB_GRAPH_H

// Local Variables:
// mode: C++
// End:
