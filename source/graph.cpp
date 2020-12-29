// graph.cpp -- A graph structure implementation.

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

#include "graph.h"
#include "heap.h"


CGraph::CGraph(int n)
{
	m_size = n;
	m_graph = new char*[n];
	m_backup = new char*[n];
	for (int i = 0; i < n; ++i)
	{
		m_graph[i] = new char[n];
		m_backup[i] = new char[n];
		for(int j=0; j<n; j++)
		{
			m_graph[i][j] = 0;	
			m_backup[i][j] = 0;
		}
	}
}

CGraph::~CGraph()
{
	if (m_graph)
	{
		for (int i = 0; i < m_size; ++i)
			delete[] m_graph[i];
		delete[] m_graph;
	}

	if (m_backup)
	{
		for (int i = 0; i < m_size; ++i)
			delete[] m_backup[i];
		delete[] m_backup;
	}
}

// This function connects two nodes in the graph.
void CGraph::connect(int var1, int var2)
{
	// Safety checks.
	assert(m_graph != NULL);
	assert((var1 >= 0) && (var1 < m_size));
	assert((var1 >= 0) && (var1 < m_size));

	m_graph[var1][var2] = 1;
	m_graph[var2][var1] = 1;
}

// This funtion connects all nodes in the set.
void CGraph::connect(variable_s& nodes)
{
	// Safety checks.
	assert(m_graph != NULL);

	int* temp = new int[nodes.size()];
	assert(temp != NULL);

	int pos = 0;
	variable_s::iterator it = nodes.begin();
	for (; it != nodes.end(); ++it)
	{
		int node = (*it);
		temp[pos++] = node;
	}

	int N = nodes.size();
	for (int i = 0; i < N - 1; ++i)
	{
		for (int j = i + 1; j < N; ++j)
		{
			int var1 = temp[i];
			int var2 = temp[j];

			if (!isConnected(var1, var2))
			{
				m_graph[var1][var2] = 3;
				m_graph[var2][var1] = 3;
			}
		}
	}

	delete[] temp;
}

// this function removes the old induced edges from the graph
void CGraph::cleanConnections()
{
	for (int i = 0; i < m_size; i++)
		for (int j = 0; j < m_size; j++)
			if (m_graph[i][j] == 3) 
				m_graph[i][j] = 0;		
}

// This function disconnects two nodes in the graph.
void CGraph::disconnect(int var1, int var2)
{
	// Safety checks.
	assert(m_graph != NULL);
	assert((var1 >= 0) && (var1 < m_size));
	assert((var1 >= 0) && (var1 < m_size));

	m_graph[var1][var2] = 2;
	m_graph[var2][var1] = 2;
}

// This function tests the connection between two nodes.
bool CGraph::isConnected(int var1, int var2)
{
	// Safety checks.
	assert(m_graph != NULL);
	assert((var1 >= 0) && (var1 < m_size));
	assert((var1 >= 0) && (var1 < m_size));

	return ((1 == m_graph[var1][var2]) || 
			(3 == m_graph[var1][var2]));
}

// This function checks if the graph is connected. It performs
// a DFS step on the original graph. If not all nodes were marked
// as visited upon DFS completion, then the graph is disconnected.
bool CGraph::isConnected()
{
	// Safety checks.
	assert(m_graph != NULL);
	assert(m_size > 0);

	bool* visited = new bool[m_size];
	for(int i=0; i<m_size; i++)
		visited[i] = false;	

	// Init search.
	stack<int> dfsStack;
	dfsStack.push(0);

	while (!dfsStack.empty())
	{
		int var = dfsStack.top();
		dfsStack.pop();

		visited[var] = true;

		// Get the not yet visited neighbors.
		variable_s neighbors;
		neighborhood(var, neighbors);

		variable_s::iterator it = neighbors.begin();
		for (; it != neighbors.end(); ++it)
		{
			int n = (*it);
			if (!visited[n])
				dfsStack.push(n);
		}

		neighbors.clear();
	}

	bool conn = true;
	for (int i = 0; i < m_size; ++i)
	{
		if (!visited[i])
		{
			conn = false;
			break;
		}
	}

	delete[] visited;

	return conn;
}

// This function creates the network graph.
void CGraph::init(function_map& functions)
{
	function_map::iterator it = functions.begin();
	for (; it != functions.end(); ++it)
	{
		CFunction* fun = (*it).second;
		if (fun->getArgc() > 1)
		{
			int argc = fun->getArgc();
			int* argv = fun->getArgv();

			for (int i = 0; i < argc - 1; ++i)
			{
				for (int j = i + 1; j < argc; ++j)
				{
					int v1 = argv[i];
					int v2 = argv[j];

					// Connect the two variables in the graph.
					connect(v1, v2);
				}
			}
		}
	}
}

// This function creates the mixed network graph.
void CGraph::init(function_map& functions, function_map& constraints)
{
	function_map::iterator it = functions.begin();
	for (; it != functions.end(); ++it)
	{
		CFunction* fun = (*it).second;
		if (fun->getArgc() > 1)
		{
			int argc = fun->getArgc();
			int* argv = fun->getArgv();

			for (int i = 0; i < argc - 1; ++i)
			{
				for (int j = i + 1; j < argc; ++j)
				{
					int v1 = argv[i];
					int v2 = argv[j];

					// Connect the two variables in the graph.
					connect(v1, v2);
				}
			}
		}
	}
	function_map::iterator iter = constraints.begin();
	for (; iter != constraints.end(); ++iter)
	{
		CFunction* fun = (*iter).second;
		if (fun->getArgc() > 1)
		{
			int argc = fun->getArgc();
			int* argv = fun->getArgv();

			for (int i = 0; i < argc - 1; ++i)
			{
				for (int j = i + 1; j < argc; ++j)
				{
					int v1 = argv[i];
					int v2 = argv[j];

					// Connect the two variables in the graph.
					connect(v1, v2);
				}
			}
		}
	}
}

void CGraph::backup()
{
	// Safety checks.
	assert(m_graph != NULL);
	assert(m_size > 0);

	for (int i = 0; i < m_size; ++i)
	{
		memcpy(m_backup[i], m_graph[i], m_size * sizeof(char));
	}
}

void CGraph::restore()
{
	// Safety checks.
	assert(m_graph != NULL);
	assert(m_size > 0);

	for (int i = 0; i < m_size; ++i)
	{
		memcpy(m_graph[i], m_backup[i], m_size * sizeof(char));
	}
}

int CGraph::neighborhood(int var, variable_s& neighbors)
{
	// Safety checks.
	assert(m_graph != NULL);
	assert((var >= 0) && (var < m_size));

	neighbors.clear();
	for (int i = 0; i < m_size; ++i)
	{
		if (isConnected(var, i))
			neighbors.insert(i);
	}

	return neighbors.size();
}

int CGraph::getMoralNeighbors(int var, variable_s& neighbors)
{
	assert(m_graph != NULL);
	assert((var >= 0) && (var < m_size));

	neighbors.clear();
	for (int i = 0; i < m_size; ++i)
	{
		// Do not check induced edges (3).
		if ((1 == m_graph[var][i]) ||
			(1 == m_graph[i][var]))
			neighbors.insert(i);
	}

	return neighbors.size();
}

int CGraph::fillset(int var, bool* used)
{
	// Safety checks.
	assert(m_graph != NULL);
	assert((var >= 0) && (var < m_size));

	int edges = 0;

	variable_v parents;
	for (int v = 0; v < m_size; ++v)
	{
		if (!used[v])
		{
			if (isConnected(var, v))
			{
				// Potential parent in the induced graph.
				parents.push_back(v);
			}
		}
	}

	// Compute the fillset edges.
	int* temp = new int[parents.size()];
	assert(temp != NULL);

	int pos = 0;
	int N = parents.size();
	variable_v::iterator it = parents.begin();
	for (; it != parents.end(); ++it)
	{
		int node = (*it);
		temp[pos++] = node;
	}

	for (int i = 0; i < N - 1; ++i)
	{
		for (int j = i + 1; j < N; ++j)
		{
			int var1 = temp[i];
			int var2 = temp[j];

			if (!isConnected(var1, var2))
				edges++;
		}
	}

	parents.clear();
	delete[] temp;

	return edges;
}

int CGraph::select(OTYPE type, bool* used)
{
	// Safety checks.
	assert(m_graph != NULL);

	int result = -1;
	int N = m_size;

	switch (type)
	{
	case mindegree:	// Select a variable with min degree.
		{
			int mindegree = INT_MAX;
			int candidate = -1;
			variable_s neighbors;

			int var;
			for (var = 0; var < N; ++var)
			{
				neighbors.clear();
				if (used[var])
					continue;	// skip used nodes.

				int degree = neighborhood(var, neighbors);
				
				if (degree < mindegree)
				{
					mindegree = degree;
					candidate = var;
				}
			}

			// Check for valid candidate.
			if (candidate >= 0)
			{
				neighbors.clear();
				neighborhood(candidate, neighbors);

				variable_s::iterator it = neighbors.begin();
				for (; it != neighbors.end(); ++it)
				{
					var = (*it);
					disconnect(candidate, var);
				}
			}

			neighbors.clear();

			// Record the result.
			result = candidate;

			break;
		}
	case minwidth:
		{
			int mindegree = INT_MAX;
			int candidate = -1;
			variable_s neighbors;

			int var;
			for (var = 0; var < N; ++var)
			{
				neighbors.clear();
				if (used[var])
					continue;	// skip used nodes.

				int degree = neighborhood(var, neighbors);
				
				if (degree < mindegree)
				{
					mindegree = degree;
					candidate = var;
				}
			}

			// Check for valid candidate.
			if (candidate >= 0)
			{
				neighbors.clear();
				neighborhood(candidate, neighbors);

				// Connect neighborhood.
				connect(neighbors);

				// Remove candidate from graph.
				variable_s::iterator it = neighbors.begin();
				for (; it != neighbors.end(); ++it)
				{
					var = (*it);
					disconnect(candidate, var);
				}
			}

			neighbors.clear();

			// Record the result.
			result = candidate;

			break;
		}
	case minfill:
		{
			int minfillset = INT_MAX;
			int candidate = -1;
			variable_s neighbors;

			int var;
			for (var = 0; var < N; ++var)
			{
				if (used[var])
					continue;	// skip used nodes.

				int fs = fillset(var, used);
				
				if (fs < minfillset)
				{
					minfillset = fs;
					candidate = var;
				}
			}

			// Check for valid candidate.
			if (candidate >= 0)
			{
				neighborhood(candidate, neighbors);

				// Connect neighborhood.
				connect(neighbors);

				// Remove candidate from graph.
				variable_s::iterator it = neighbors.begin();
				for (; it != neighbors.end(); ++it)
				{
					var = (*it);
					disconnect(candidate, var);
				}
			}

			neighbors.clear();

			// Record the result.
			result = candidate;

			break;
		}
	case maxcard:
		{
			int maxEdges = -1;
			int candidate = -1;

			for (int var = 0; var < N; ++var)
			{
				if (used[var])
					continue;	// skip used nodes.

				int edges = 0;
				for (int cand = 0; cand < N; ++cand)
				{
					if (!used[cand])
						continue;	// skip not used nodes.

					if (isConnected(var, cand))
						edges++;
				}

				if (edges > maxEdges)
				{
					maxEdges = edges;
					candidate = var;
				}
			}

			// Record the result
			result = candidate;

			break;
		}
	}

	return result;
}

// This function computes the induced width of the ordering.
int CGraph::width(int* ordering, int* position)
{
	// Safety checks.
	assert(m_graph != NULL);
	assert(ordering != NULL);
	assert(position != NULL);
	assert(m_size > 0);

	int maxdegree = -1;
	variable_s parents;
	int N = m_size;
	for (int i = N - 1; i >= 1; --i)
	{
		int var = ordering[i];
		parents.clear();

		for (int j = i - 1; j >= 0; --j)
		{
			int parent = ordering[j];
			if (isConnected(var, parent))
			{
				parents.insert(parent);
			}
		}

		int degree = parents.size();
		if (degree > maxdegree)
		{
			maxdegree = degree;
		}

		// Connect parents.
		connect(parents);
	}

	parents.clear();

	return maxdegree;
}

// This function computes the "DEFAULT" variable ordering.
void CGraph::orderingDefault(int*& ordering, int*& position)
{
	// Safety checks.
	assert(m_graph != NULL);
	assert(m_size > 0);

	for (int i = 0; i < m_size; ++i)
	{
		ordering[i] = i;
		position[i] = i;
	}

}

// This function computes the "MIN-DEGREE" variable ordering.
void CGraph::orderingMinDegree(int*& ordering, int*& position)
{
	// Safety checks.
	assert(m_graph != NULL);
	assert(m_size > 0);

	// Declare variables.
	int N = m_size;
	bool* used = new bool[N];
	for(int i=0; i<N; i++)
		used[i] = false;	

	int pos = N - 1;
	int cnt = 1;
	while (cnt <= N)
	{
		// Look for next variable.
		int var = select(mindegree, used);
		assert(var >= 0);

		// Mark variable as used.
		used[var] = true;

		// Update position/order.
		ordering[pos] = var;
		position[var] = pos;

		--pos;
		++cnt;
	}

	delete[] used;
}

// This function computes the "MIN-INDUCED-WIDTH" variable ordering.
void CGraph::orderingMinWidth(int*& ordering, int*& position)
{
	// Safety checks.
	assert(m_graph != NULL);
	assert(m_size > 0);

	// Declare variables.
	int N = m_size;
	bool* used = new bool[N];
	for(int i=0; i<N; i++)
		used[i] = false;	

	int pos = N - 1;
	int cnt = 1;

	while (cnt <= N)
	{
		// Look for next variable.
		int var = select(minwidth, used);
		assert(var >= 0);

		// Mark variable as used.
		used[var] = true;

		// Update position/order.
		ordering[pos] = var;
		position[var] = pos;

		--pos;
		++cnt;
	}

	delete[] used;
}

// This function computes the "MIN-FILL" variable ordering.
void CGraph::orderingMinFill(int*& ordering, int*& position)
{
	// Safety checks.
	assert(m_graph != NULL);
	assert(m_size > 0);

	// Declare variables.
	int N = m_size;
	bool* used = new bool[N];
	for(int i=0; i<N; i++)
		used[i] = false;	

	int pos = N - 1;
	int cnt = 1;

	while (cnt <= N)
	{
		// Look for next variable.
		int var = select(minfill, used);
		assert(var >= 0);

		// Mark variable as used.
		used[var] = true;

		// Update position/order.
		ordering[pos] = var;
		position[var] = pos;

		--pos;
		++cnt;

		//debug
//		cout << "variable: " << cnt << "\n";
	}

	delete[] used;
}

// This function computes the "MAX-CARD" variable ordering.
void CGraph::orderingMaxCard(int*& ordering, int*& position)
{
	// Safety checks.
	assert(m_graph != NULL);
	assert(m_size > 0);

	// Declare variables.
	int N = m_size;
	bool* used = new bool[N];
	for(int i=0; i<N; i++)
		used[i] = false;	

	// Initialize the ordering with max-degree variable.
	int maxDegree = -1;
	int firstVar = -1;
	variable_s neighbors;
	for (int var = 0; var < N; ++var)
	{
		int degree = neighborhood(var, neighbors);
		if (degree > maxDegree)
		{
			maxDegree = degree;
			firstVar = var;
		}
	}

	neighbors.clear();
	ordering[0] = firstVar;
	position[firstVar] = 0;
	used[firstVar] = true;

	int pos = 1;
	int cnt = 2;
	while (cnt <= N)
	{
		// Look for next variable.
		int var = select(maxcard, used);
		assert(var >= 0);

		// Mark variable as used.
		used[var] = true;

		// Update position/order.
		ordering[pos] = var;
		position[var] = pos;

		++pos;
		++cnt;
	}

	delete[] used;
}


// This function computes the "RANDOMIZED" variable ordering (Superlink - style).
void CGraph::orderingRandomized(int*& ordering, int*& position)
{
	// Safety checks.
	assert(m_graph != NULL);
	assert(m_size > 0);

	// Declare variables.
	int N = m_size;
	bool* used = new bool[N];
	for(int i=0; i<N; i++)
		used[i] = false;	
}

// This function computes a variable ordering.
void CGraph::order(int type, int*& ordering, int*& position)
{
	// Safety checks.
	assert(m_graph != NULL);
	assert(m_size > 0);

	if (ordering)
	{
		delete[] ordering;
		ordering = NULL;
	}

	if (position)
	{
		delete[] position;
		position = NULL;
	}

	// Allocate mmory.
	ordering = new int[m_size];
	position = new int[m_size];
	assert(ordering != NULL);
	assert(position != NULL);

	// Backup the original graph connections.
	backup();

	// Compute and record the ordering.
	switch (type)
	{
	case VO_DEFAULT:
		{
			orderingDefault(ordering, position);
			break;
		}
	case VO_MINDEGREE:
		{
			orderingMinDegree(ordering, position);
			break;
		}
	case VO_MINWIDTH:
		{
			orderingMinWidth(ordering, position);
			break;
		}
	case VO_MINFILL:
		{
			orderingMinFill(ordering, position);
			break;
		}
	case VO_MAXCARD:
		{
			orderingMaxCard(ordering, position);
			break;
		}
	case VO_RANDOMIZED:
		{
			orderingRandomized(ordering, position);
			break;
		}
	}

	// Restore graph connections to original values.
	restore();

	// Create the induced graph/width of the ordering.
	createInduced(ordering, position);
}

// This function computes the induced graph/width of the ordering.
void CGraph::createInduced(int* ordering, int* position)
{
	// Safety checks.
	assert(m_graph != NULL);
	assert(ordering != NULL);
	assert(position != NULL);
	assert(m_size > 0);

	int maxdegree = -1;
	variable_s parents;
	int N = m_size;
	for (int i = N - 1; i >= 1; --i)
	{
		int var = ordering[i];
		parents.clear();

		for (int j = i - 1; j >= 0; --j)
		{
			int parent = ordering[j];
			if (isEdge(var, parent))
			{
				parents.insert(parent);
			}
		}

		int degree = (int)parents.size();
		if (degree > maxdegree)
		{
			maxdegree = degree;
		}

		// Connect parents.
		connect(parents);
	}

	parents.clear();

	// Update the induced width.
	m_width = maxdegree;
	m_induced = true;
}

CGraph* CGraph::clone()
{
	CGraph* cl = new CGraph(m_size);
	cl->m_induced = m_induced;
	cl->m_width = m_width;
	for (int i = 0; i < m_size; ++i)
	{
		memcpy(cl->m_graph[i], m_graph[i], m_size * sizeof(char));
		memcpy(cl->m_backup[i], m_backup[i], m_size * sizeof(char));
	}

	return cl;
}

// Cliques are extracted following the elimination order.
void CGraph::extractCliques(int* ordering, vector<set<int> >& cliques)
{
	// Allocate space for cliques.
	cliques.resize(m_size);

	int i, j, var, parent;
	for (i = m_size - 1; i >= 0; --i)
	{
		var = ordering[i];
		cliques[i].insert(var);

		for (j = i - 1; j >= 0; --j)
		{
			parent = ordering[j];
			if (isEdge(var, parent))
				cliques[i].insert(parent);
		}
	}
}

// This function tests the connection between two nodes.
bool CGraph::isEdge(int var1, int var2, bool original)
{
	// Safety checks.
	assert(m_graph != NULL);
	assert((var1 >= 0) && (var1 < m_size));
	assert((var1 >= 0) && (var1 < m_size));

	if (!original)
		return ((1 == m_graph[var1][var2]) || 
				(3 == m_graph[var1][var2]));
	else
		return (1 == m_graph[var1][var2]);
}


////////////////////////////////////////////////////////////////
// from Radu

////////////////////////////////////////////////////////////////////
// CGraphHash class implementation.

CGraphHash::CGraphHash()
{
}

// This is the copy constructor.
CGraphHash::CGraphHash(CGraphHash& g)
{
	set<int>& vertices = g.m_vertices;
	map<int,variable_s*>& neighbors = g.m_neighbors;

	// Copy the vertices.
	copy(vertices.begin(), vertices.end(), 
		inserter(m_vertices, m_vertices.begin()));

	// Copy the neighbors.
	map<int,variable_s*>::iterator it = neighbors.begin();
	for (; it != neighbors.end(); ++it)
	{
		int i = (*it).first;
		variable_s* s = (*it).second;

		m_neighbors.insert(make_pair(i, new variable_s(*s)));
	}
}

// This is the = operator.
CGraphHash& CGraphHash::operator=(const CGraphHash& g)
{
	set<int>& vertices = ((CGraphHash&)g).m_vertices;
	map<int,variable_s*>& neighbors = ((CGraphHash&)g).m_neighbors;

	// Copy the vertices.
	copy(vertices.begin(), vertices.end(), 
		inserter(m_vertices, m_vertices.begin()));

	// Copy the neighbors.
	map<int,variable_s*>::iterator it = neighbors.begin();
	for (; it != neighbors.end(); ++it)
	{
		int i = (*it).first;
		variable_s* s = (*it).second;

		m_neighbors.insert(make_pair(i, new variable_s(*s)));
	}

	return (*this);
}

CGraphHash::~CGraphHash()
{
	map<int,variable_s*>::iterator it = m_neighbors.begin();
	for (; it != m_neighbors.end(); ++it)
	{
		variable_s* s = (*it).second;
		s->clear();
		delete s;
	}

	m_vertices.clear();
	m_neighbors.clear();
}

void CGraphHash::clear()
{
	map<int,variable_s*>::iterator it = m_neighbors.begin();
	for (; it != m_neighbors.end(); ++it)
	{
		variable_s* s = (*it).second;
		s->clear();
		delete s;
	}

	m_vertices.clear();
	m_neighbors.clear();
}

variable_s* CGraphHash::neighbors(int var)
{
	map<int,variable_s*>::iterator it = m_neighbors.find(var);
	if (it != m_neighbors.end())
		return (*it).second;

	return NULL;
}

// This function adds a new node to the graph.
void CGraphHash::addNode(int i)
{
	if (m_vertices.find(i) == m_vertices.end())
	{
		m_vertices.insert(i);
		m_neighbors.insert(make_pair(i, new variable_s()));
	}
}

// This function adds the edge (i,j) to the graph. Also adds the reversed edge.
void CGraphHash::addEdge(int i, int j)
{
	map<int,variable_s* >::iterator it;
	
	// Add edge (i,j) to the graph.
	it = m_neighbors.find(i);
	if (it != m_neighbors.end())
	{
		variable_s* s = (*it).second;
		if (s->find(j) == s->end())
			s->insert(j);
	}

	// Add edge (j,i) to the graph.
	it = m_neighbors.find(j);
	if (it != m_neighbors.end())
	{
		variable_s* s = (*it).second;
		if (s->find(i) == s->end())
			s->insert(i);
	}
}

// This function adds the edge (i,j) to the graph. Does not add the reversed edge.
void CGraphHash::addDirectedEdge(int i, int j)
{
	map<int,variable_s* >::iterator it;
	
	// Add edge (i,j) to the graph.
	it = m_neighbors.find(i);
	if (it != m_neighbors.end())
	{
		variable_s* s = (*it).second;
		if (s->find(j) == s->end())
			s->insert(j);
	}
}

bool CGraphHash::containsEdge(int i, int j)
{
	variable_s* s = neighbors(i);
	if (s->find(j) != s->end())
		return true;
	else
		return false;
}

// This function removes node i from the graph.
void CGraphHash::removeNode(int i)
{
	map<int,variable_s*>::iterator ita, itb;
	ita = m_neighbors.find(i);
	if (ita != m_neighbors.end())
	{
		variable_s* neighborset = (*ita).second;
		for (set<int>::iterator its = neighborset->begin();
			its != neighborset->end(); ++its)
		{
			int j = (*its);
			itb = m_neighbors.find(j);
			if (itb != m_neighbors.end())
			{
				variable_s* s = (*itb).second;
				s->erase(i);
			}
		}

		neighborset->clear();
		delete neighborset;
		m_neighbors.erase(i);
	}

	m_vertices.erase(i);
}

// This function removes the edge (i,j) from the graph.
void CGraphHash::removeEdge(int i, int j)
{
	map<int,variable_s* >::iterator it;
	variable_s* s;

	// Remove node "j" from i's neighbors.
	it = m_neighbors.find(i);
	s = (*it).second;
	s->erase(j);

	// Remove node "i" from j's neighbors.
	it = m_neighbors.find(j);
	s = (*it).second;
	s->erase(i);
}

// This function adds a new clique to the graph.
void CGraphHash::addClique(vector<int>& clique)
{
	// Add vertices first.
	int i, j, n = (int)clique.size();
	for (i = 0; i < n; ++i)
		addNode(clique[i]);

	// Add edges next.
	for (i = 0; i < n - 1; ++i)
		for (j = i + 1; j < n; ++j)
			addEdge(clique[i], clique[j]);
}

void CGraphHash::addClique(int n, int* clique)
{
	// Add vertices first.
	int i, j;
	for (i = 0; i < n; ++i)
		addNode(clique[i]);

	// Add edges next.
	for (i = 0; i < n - 1; ++i)
		for (j = i + 1; j < n; ++j)
			addEdge(clique[i], clique[j]);
}

void CGraphHash::addCliqueSet(set<int>& clique)
{
	vector<int> temp;
	copy(clique.begin(), clique.end(), back_inserter(temp));

	addClique(temp);
	temp.clear();
}

// This function returns the degree of a node in the graph.
int CGraphHash::degree(int i)
{
	map<int,variable_s*>::iterator it;
	it = m_neighbors.find(i);
	if (it != m_neighbors.end())
	{
		variable_s* s = (*it).second;
		return (int)s->size();
	}

	return 0;
}

double CGraphHash::scoreFunction(int etype, int var, int* domains)
{
	double score;

	switch (etype)
	{
	case VO_MINDEGREE:
		score = scoreMinDegree(var);
		break;
	case VO_MINFILL:
	case VO_DETERMINISM:
		score = scoreMinFill(var);
		break;
	case VO_MINWIDTH:
		score = scoreMinWidth(var);
		break;
	case VO_MINWEIGHT:
		score = scoreMinWeight(var, domains);
		break;
	case VO_MAXDOMAINMINFILL:
		score = scoreMaxDomainMinFill(var, domains);
		break;
	case VO_WEIGHTEDMINFILL:
		score = scoreWeightedMinFill(var, domains);
		break;

	default:
		score = 0;
		break;
	};

	return score;
}

double CGraphHash::scoreMinDegree(int var)
{
	return degree(var);
}

double CGraphHash::scoreMinWidth(int var)
{
	return degree(var);
}

double CGraphHash::scoreMinFill(int var)
{
	variable_s* s = neighbors(var);
	variable_s::iterator it1, it2;
	double total = 0;

	for (it1 = s->begin(); it1 != s->end(); ++it1)
	{
		int var1 = (*it1);
		for (it2 = s->begin(); it2 != s->end(); ++it2)
		{
			int var2 = (*it2);
			if ((var1 != var2) && !containsEdge(var1, var2))
				total ++;
		}
	}

	return total;
}

double CGraphHash::scoreMinWeight(int var, int* domains)
{
	variable_s* s = neighbors(var);
	variable_s::iterator it;
	double total = 0;

	for (it = s->begin(); it != s->end(); ++it)
	{
		int var1 = (*it);
		total *= domains[var1];
	}

	return total;
}

double CGraphHash::scoreMaxDomainMinFill(int var, int*domains)
{
	variable_s* s = neighbors(var);
	variable_s::iterator it1, it2;
	double total = 0;

	for (it1 = s->begin(); it1 != s->end(); ++it1)
	{
		int var1 = (*it1);
		for (it2 = s->begin(); it2 != s->end(); ++it2)
		{
			int var2 = (*it2);
			if ((var1 != var2) && !containsEdge(var1, var2))
				total++;
		}
	}

	// now give huge weight to big domains 
	return total + (1000000.0 / (double)domains[var]);
}

double CGraphHash::scoreWeightedMinFill(int var, int *domains)
{
	variable_s* s = neighbors(var);
	variable_s::iterator it1, it2;
	double total = 0;

	for (it1 = s->begin(); it1 != s->end(); ++it1)
	{
		int var1 = (*it1);
		for (it2 = s->begin(); it2 != s->end(); ++it2)
		{
			int var2 = (*it2);
			if ((var1 != var2) && !containsEdge(var1, var2))
				total++;
		}
	}

	return total;
}

// This function removes a node from the graph and connects its neighbors.
void CGraphHash::removeAndConnect(int i)
{
	map<int,variable_s*>::iterator it;
	it = m_neighbors.find(i);

	// Connect its neihbors.
	variable_s* s = (*it).second;
	addCliqueSet(*s);

	// Remove the node from the graph.
	removeNode(i);
}
// Removes a node from the graph and adds only ancestral edges (parent to other neighbor) - for deterministic variables
void CGraphHash::removeAndConnectAncestral(int i, variable_s &parents)
{
	map<int,variable_s*>::iterator it;
	it = m_neighbors.find(i);

	// Get set of neighbors
	variable_s* s = (*it).second;

	//adjust parents to current graph (intersect with current neighborhood)
	variable_s tempParents;
	variable_s::iterator it1 = parents.begin();
	for ( ; it1 != parents.end() ; it1++)
	{
		if(s->find(*it1) != s->end())
			tempParents.insert(*it1);
	}

	//// safety check
	//variable_s::iterator ittt = parents.begin();
	//for ( ; ittt != parents.end() ; ittt++)
	//	if (s->find(*ittt) == s->end() )
	//		cout << "uuuuu  ";

	// add only ancestral edges
	variable_s::iterator itv = s->begin();
	for ( ; itv != s->end(); ++itv )
	{
		int var = (*itv);
		if ( tempParents.find(var) == tempParents.end() )
		{
			variable_v temp;
			copy(tempParents.begin(), tempParents.end(), back_inserter(temp));
			temp.push_back(var);
			addClique(temp);
			temp.clear();			
		}
	}

	tempParents.clear();

	// Remove the node from the graph.
	removeNode(i);
}

// This function finds the elimination order of the graph.
void CGraphHash::eliminate(int etype, vector<int>& order, int* domains, bool* deterministic, variable_s* parents)
{
	// Create a copy of the current graph.
	CGraphHash g(*this);
	
	int i, v, var, n = (int)m_vertices.size();
	vector<double> scores;
	vector<void*> nodes;
	map<int,CGraphNode*> mapNodes;
	map<int,CGraphNode*>::iterator itn;
	CGraphNode* node = NULL;

	set<int>::iterator it = m_vertices.begin();
	for (; it != m_vertices.end(); ++it)
	{
		var = (*it);
		node = new CGraphNode(var);
		nodes.push_back(node);
		mapNodes.insert(make_pair(var, node));
	/*	if(etype == VO_DETERMINISM && deterministic[var])
			scores.push_back(-10000000000);
		else*/
			scores.push_back(-scoreFunction(etype, var, domains));
	}

	CHeap heap(nodes, scores);
	for (i = 0; i < n; ++i)
	{
		CHeapElement* he = heap.top();
		heap.pop();

		var = ((CGraphNode*)(he->m_object))->m_variable;
		order.push_back(var);
		variable_s s(*(g.neighbors(var)));

		// debug
		// cout << "\n var = " << var ;

		switch (etype)
		{
		case VO_MINFILL:
		case VO_MINWIDTH:
		case VO_MAXDOMAINMINFILL:
		case VO_WEIGHTEDMINFILL:
            g.removeAndConnect(var);
			break;
		case VO_DETERMINISM:
			if (deterministic[var])				
				g.removeAndConnectAncestral(var, parents[var]);
			else
				g.removeAndConnect(var);				
			break;
		case VO_MINDEGREE:
		case VO_MINWEIGHT:
			g.removeNode(var);
			break;
		}

		// Update the scores of its neighbors.
		for (it = s.begin(); it != s.end(); ++it)
		{
			v = (*it);
			if (etype == VO_DETERMINISM && deterministic[v])
				continue;

			itn = mapNodes.find(v);
			node = (*itn).second;
			heap.updateScore(node, -(g.scoreFunction(etype, v, domains)));
		}

		s.clear();
	}

	// Free memory.
	scores.clear();
	mapNodes.clear();
	for (i = 0; i < (int)nodes.size(); ++i)
		delete nodes[i];
	nodes.clear();
}

// This function finds the elimination order of the graph, starting with root.
void CGraphHash::eliminate(int root, int etype, vector<int>& order, int* domains)
{
	// Create a copy of the current graph.
	CGraphHash g(*this);
	
	int i, v, var, n = (int)m_vertices.size();
	vector<double> scores;
	vector<void*> nodes;
	map<int,CGraphNode*> mapNodes;
	map<int,CGraphNode*>::iterator itn;
	CGraphNode* node = NULL;

	set<int>::iterator it = m_vertices.begin();
	for (; it != m_vertices.end(); ++it)
	{
		var = (*it);
		node = new CGraphNode(var);
		nodes.push_back(node);
		mapNodes.insert(make_pair(var, node));
		
		// Make sure that the root starts the ordering.
		if (var != root)
			scores.push_back(-scoreFunction(etype, var, domains));
		else
			scores.push_back(MINUS_INFINITY);
	}

	CHeap heap(nodes, scores);
	for (i = 0; i < n; ++i)
	{
		CHeapElement* he = heap.top();
		heap.pop();

		var = ((CGraphNode*)(he->m_object))->m_variable;
		order.push_back(var);
		variable_s s(*(g.neighbors(var)));

		switch (etype)
		{
		case VO_MINFILL:
		case VO_MINWIDTH:
            g.removeAndConnect(var);
			break;
		case VO_MINDEGREE:
			g.removeNode(var);
			break;
		}

		// Update the scores of its neighbors.
		for (it = s.begin(); it != s.end(); ++it)
		{
			v = (*it);
			itn = mapNodes.find(v);
			node = (*itn).second;

			if (v != root)
				heap.updateScore(node, -(g.scoreFunction(etype, v, domains)));
			else
				heap.updateScore(node, MINUS_INFINITY);
		}

		s.clear();
		delete he;
	}

	// Check for consistecy.
	if (find(order.begin(), order.end(), root) != order.end())
	{
		int& last = order.back();
		assert(last == root);
	}

	// Free memory.
	scores.clear();
	mapNodes.clear();
	for (i = 0; i < (int)nodes.size(); ++i)
		delete nodes[i];
	nodes.clear();
}

void CGraphHash::induce(vector<int>& order, int& width)
{
	width = 0;
	CGraphHash g(*this);
	vector<int>::iterator it = order.begin();
	for (; it != order.end(); ++it)
	{
		int var = (*it);
		int n = g.degree(var);
		if (n > width)
			width = n;

		variable_s s(*(g.neighbors(var)));
		g.removeNode(var);
		g.addCliqueSet(s);

		s.clear();
	}
}

// This function checks if the graph is connected. It performs
// a DFS step on the original graph. If not all nodes were marked
// as visited upon DFS completion, then the graph is disconnected.
bool CGraphHash::isConnected()
{
	int n = (int)m_vertices.size();
	bool* visited = new bool[n];
	for(int i=0; i<n; i++)
		visited[i] = false;	

	// Init search.
	stack<int> dfsStack;
	dfsStack.push(0);

	while (!dfsStack.empty())
	{
		int var = dfsStack.top();
		dfsStack.pop();

		visited[var] = true;

		// Get the not yet visited neighbors.
		variable_s* s = neighbors(var);
		variable_s::iterator it = s->begin();
		for (; it != s->end(); ++it)
		{
			int n = (*it);
			if (!visited[n])
				dfsStack.push(n);
		}
	}

	bool conn = true;
	for (int i = 0; i < n; ++i)
	{
		if (!visited[i])
		{
			conn = false;
			break;
		}
	}

	delete[] visited;

	return conn;
}

// This function creates the graph links from the subproblem's functions.
void CGraphHash::init(function_v& funs, map<int,int>& idxVarToKey, set<int>& ancestors)
{
	int argc, *argv, i;
	variable_v vars;

	function_v::iterator it = funs.begin();
	for (; it != funs.end(); ++it)
	{
		CFunction* fun = (*it);
		argc = fun->getArgc();
		argv = fun->getArgv();
		
		// Identify non-ancestor variables in the scope.
		vars.clear();
		for (i = 0; i < argc; ++i)
		{
			if (ancestors.find(argv[i]) == ancestors.end())
			{
				map<int,int>::iterator itv = idxVarToKey.find(argv[i]);
				vars.push_back((*itv).second);
			}
		}

		// Connect the variables.
		addClique(vars);
	}

	vars.clear();
}

void CGraphHash::addAncestralEdges(variable_s* parents)
{
	int i, count = 0;
	
	queue<int> Queue;
	bool* inQueue = new bool[m_vertices.size()];

	for (int i=0; i<m_vertices.size(); i++)
	{
		Queue.push(i);
		inQueue[i] = true;
	}

	while (!Queue.empty())
	{
		i = Queue.front();
		Queue.pop();
		inQueue[i] = false;

		// non-deterministic vars don't have parents assigned
		if ( parents[i].empty() )
			continue;

		map<int,variable_s*>::iterator it;
		it = m_neighbors.find(i);

		// Get set of neighbors
		variable_s* s = (*it).second;

		// add only ancestral edges
		variable_s::iterator itv = s->begin();
		for ( ; itv != s->end(); ++itv )
		{
			int var = (*itv);
			if ( parents[i].find(var) == parents[i].end() )
			{
				// check all possible edges between parents and var
				variable_s::iterator itp = parents[i].begin();
				for ( ; itp != parents[i].end(); itp++ )
				{
					int p = (*itp);
					if (!containsEdge(var, p))
					{
						addEdge(var, p);
						if (!inQueue[var])
						{
							Queue.push(var);
							inQueue[var] = true;
						}
						if (!inQueue[p])
						{
							Queue.push(p);
							inQueue[p] = true;
						}
						cout << count++ << endl;
					}
				}
			}
		}

	}// end while

	delete[] inQueue;
}

// Local Variables:
// mode: C++
// End:
