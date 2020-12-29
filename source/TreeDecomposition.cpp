// TreeDecomposition.cpp: implementation of the CTreeDecomposition class.
//
//////////////////////////////////////////////////////////////////////

#include "TreeDecomposition.h"
#include "defs.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CTreeDecomposition::CTreeDecomposition(CProblem* owner)
	: m_owner(owner)
{
	m_size = 0;
	m_cutsetTreeNode = NULL;
}

CTreeDecomposition::~CTreeDecomposition()
{
	destroy();
}

// This function destroys the cluster structure.
void CTreeDecomposition::destroy()
{
	treeDecompositionCluster_v::iterator it = m_clusters.begin();
	for ( ; it != m_clusters.end() ; it++)
		delete (*it);
	m_clusters.clear();

	treeDecompositionSeparator_v::iterator its = m_separators.begin();
	for ( ; its != m_separators.end() ; its++)
		delete (*its);
	m_separators.clear();
	m_root = NULL;
	m_size = 0;
}

// print the Cluster structure
void CTreeDecomposition::print()
{
	cout << "\n Tree Decomposition:\n";
	cout << "\n Clusters:\n" ;
	treeDecompositionCluster_v::reverse_iterator it = m_clusters.rbegin();
	for ( ; it != m_clusters.rend() ; it++)
	{
		(*it)->print();
	}
	//(*it)->print();

	cout << "\n Separators:\n" ;
	treeDecompositionSeparator_v::reverse_iterator its = m_separators.rbegin();
	for ( ; its != m_separators.rend() ; its++)
	{
		(*its)->print();
	}
	//(*its)->print();
}

// function to compute the separtor weights, to decide where to break the tree decomposition
// can only be called when the tree decomposition is well defined
// for each separator, it computes a < 1 weight, w_left/w_right or w_right/w_left, and picks the largest
// w_right = Sum_{i=1}^{log(number_of_clusters)}  [higest_w(cluster)],
// where cluster is in the right portion after the separator is eliminated
int CTreeDecomposition::computeSeparatorWeights(int scopeSize)
{
	//int count = 0;
	if (m_separators.empty())
		return 0;

	treeDecompositionSeparator_v::iterator it = m_separators.begin();
	for ( ; it != m_separators.end() ; it++ )
	{
		(*it)->computeWeight(scopeSize);
	//	cout << "sep: " << count++ << "\n";
	}

	return 1;
}


double CTreeDecomposition::computeSubtreeWeight(CTreeDecompositionCluster* root, CTreeDecompositionCluster* cutoff, variable_s variables, int scopeSize)
{
	// count the number of clusters in the subtree
	// and the max {cluster - variables} size
	int nClust = 0;
	int maxClustWidth = -1; // cluster.size() - variables.size()


	treeDecompositionCluster_v clustQueue;
	clustQueue.push_back(root);

	// maintain an ordered set of clusters sizes (when separator vars are removed) and then pick the first
	// log( # clusters ) from it
	multiset<int, greater<int> > orderedValues;

	int count = 0; //debug
	while(!clustQueue.empty())
	{
		//cout << clustQueue.size() << "\n";
		CTreeDecompositionCluster* curCluster = clustQueue.back();
		clustQueue.pop_back();
		if (curCluster == cutoff)
			continue;

		variable_s tempVariables;
		set_difference(	curCluster->getVariables().begin(), curCluster->getVariables().end(),
						variables.begin(), variables.end(), inserter(tempVariables, tempVariables.begin())
					   );

		int tempSize = tempVariables.size();
		orderedValues.insert(tempSize);

//		if (tempSize > maxClustWidth)
//			maxClustWidth = tempSize;

		//
		if (tempSize > scopeSize)
			nClust++;

		treeDecompositionCluster_v::iterator it = curCluster->children().begin();
		for ( ; it != curCluster->children().end() ; it++ )
			clustQueue.push_back(*it);

		tempVariables.clear();

		//cout << "      " << count++ << "\n";
	}

	double weight = 0;
	for (int i = 0 ; i < nClust; i++) // ( log((double) nClust) / log(2.0) ); i++) // 
	{
		multiset<int, greater<int> >::iterator it = orderedValues.begin();
		int val = (*it);
		if( val > scopeSize )
			weight += val - scopeSize;
		//else
		//	weight += 1;
		orderedValues.erase(it);
	}

	//test
	//
//	if (maxClustWidth < w)
//		 return 1;

//	if (nClust == 0)
//		return 1;

//	weight = maxClustWidth * log((double) nClust);
//	if (weight == 0)
//		return 1;
//	else

	orderedValues.clear();

	if (weight == 0)
		weight = 1;
	return weight;
}

CTreeDecompositionSeparator* CTreeDecomposition::highestWeightSeparator()
{
	double maxWeight = -1;
	CTreeDecompositionSeparator* bestSep = NULL;

	treeDecompositionSeparator_v::iterator it = m_separators.begin();
	for ( ; it != m_separators.end(); it++ )
	{
		if ( (*it)->weight() > maxWeight )
		{
			bestSep = (*it);
			maxWeight = (*it)->weight();
		}
	}

	assert (bestSep != NULL);


	return bestSep;
}

int CTreeDecomposition::computeRemainingWidth(variable_s& variables)
{
	int maxSize = -1;

	treeDecompositionCluster_v::iterator it = m_clusters.begin();
	for ( ; it != m_clusters.end() ; it++ )
	{
		variable_s tempVariables;
		set_difference(	(*it)->getVariables().begin(), (*it)->getVariables().end(),
						variables.begin(), variables.end(), inserter(tempVariables, tempVariables.begin())
					   );

		int tempSize = tempVariables.size();
		if (tempSize > maxSize)
			maxSize = tempSize;

		tempVariables.clear();
	}
	
	return maxSize;
}

int CTreeDecomposition::computeRemainingWidth(int var)
{
	int maxSize = -1;

	treeDecompositionCluster_v::iterator it = m_clusters.begin();
	for ( ; it != m_clusters.end() ; it++ )
	{	
		if ( (*it)->getVariables().empty() )
			continue;

		if ( (*it)->isMemberOf(var) )
		{
			int tempSize = (*it)->getVariables().size() - 1;
			if (tempSize > maxSize)
				maxSize = tempSize;
		}
		else
		{
			if ( (*it)->getVariables().size() > maxSize )
				maxSize = (*it)->getVariables().size();
		}
	}
	
	return maxSize;
}


int CTreeDecomposition::computeBestVariableGWC(int w)
{
	// set two arrays of m_N variables, to keep track of clusters > w+1, and total number of clusters
	int N = m_owner->getN();

	int* clustersGreaterThanWplusOne = new int[N];
	assert(clustersGreaterThanWplusOne != NULL);
	memset(clustersGreaterThanWplusOne, 0, N * sizeof(int));

	int* totalNumberOfClusters = new int[N];
	assert(totalNumberOfClusters != NULL);
	memset(totalNumberOfClusters, 0, N * sizeof(int));

	treeDecompositionCluster_v::iterator it = m_clusters.begin();
	for ( ; it != m_clusters.end() ; it++ )
	{
		CTreeDecompositionCluster* cluster = (*it);
		int size = cluster->getVariables().size();

		variable_s::iterator its = cluster->getVariables().begin();
		for ( ; its != cluster->getVariables().end() ; its++ )
		{
			if ( size > w )
				clustersGreaterThanWplusOne[(*its)] += 1;
			totalNumberOfClusters[(*its)] += 1;
		}
	}

	int var = -1;
	int mostBigClusters = -1;
	int mostClusters = -1;

	for ( int i = 0; i < N ; i++ )
	{
		if ( clustersGreaterThanWplusOne[i] > mostBigClusters )
		{
			var = i;
			mostBigClusters = clustersGreaterThanWplusOne[i];
			mostClusters = totalNumberOfClusters[i];
		}
		else
		{
			if ( (clustersGreaterThanWplusOne[i] == mostBigClusters) && (totalNumberOfClusters[i] > mostClusters) )
			{
				var = i;				
				mostClusters = totalNumberOfClusters[i];
			}
		}
	}

	delete[] clustersGreaterThanWplusOne;
	delete[] totalNumberOfClusters;

	assert(var != -1);

	return var;
}

int CTreeDecomposition::connectedComponents()
{
	// mark initial state
	treeDecompositionCluster_v::iterator it = m_clusters.begin();
	for ( ; it != m_clusters.end() ; it++ )
		(*it)->setVisited(INITIAL_STATE);

	int components = 0; // the number of the component

	// start marking components
	it = m_clusters.begin();
	for ( ; it != m_clusters.end() ; it++ )
	{
		if ( (*it)->getVisited() != INITIAL_STATE )
			continue;

		if ( (*it)->getVariables().empty() )
			continue;

		(*it)->setVisited(components);
		treeDecompositionCluster_v Stack;
		Stack.push_back( (*it) );
		
		while ( !Stack.empty() )
		{
			CTreeDecompositionCluster* currClust = Stack.back();
			Stack.pop_back();

			treeDecompositionCluster_v::iterator it_ch = currClust->children().begin();
			for ( ; it_ch != currClust->children().end() ; it_ch ++)
			{
				if ( !getSeparator(currClust, (*it_ch))->getVariables().empty() )
				{
					(*it_ch)->setVisited(components) ;
					Stack.push_back(*it_ch);
				}
			}
		}

		components++;
	}

	return components;
}


void CTreeDecomposition::getComponent(int component, treeDecompositionCluster_v &clusters)
{
	treeDecompositionCluster_v::iterator it = m_clusters.begin();
	for ( ; it != m_clusters.end() ; it++ )
	{
		if ( (*it)->getVisited() != component )
			continue;

		CTreeDecompositionCluster* tempClust = new CTreeDecompositionCluster( (*it)->getID(), (*it)->getVariables() );
		clusters.push_back(tempClust);		
	}
}

void CTreeDecomposition::eliminateVar(int var)
{
	treeDecompositionCluster_v::iterator it = m_clusters.begin();
	for ( ; it != m_clusters.end() ; it++ )
		(*it)->eliminateVar(var);

	treeDecompositionSeparator_v::iterator its = m_separators.begin();
	for ( ; its != m_separators.end() ; its++ )
		(*its)->eliminateVar(var);
}

CTreeDecompositionSeparator* CTreeDecomposition::getSeparator(CTreeDecompositionCluster* parent, CTreeDecompositionCluster* child)
{
	treeDecompositionSeparator_v::iterator it = m_separators.begin();
	for ( ; it != m_separators.end() ; it++ )
	{
		CTreeDecompositionCluster* par = (*it)->getClusterParent();
		CTreeDecompositionCluster* chi = (*it)->getClusterChild();
		if ( ( par == parent ) && ( chi == child ) )
			return (*it);
	}

	return NULL;
}

bool CTreeDecomposition::cyclic()
{
	treeDecompositionCluster_v::iterator it = m_clusters.begin();
	for ( ; it != m_clusters.end() ; it++ )
		(*it)->setVisited(INITIAL_STATE);

	treeDecompositionCluster_v Stack;
	it = m_clusters.begin();
	Stack.push_back( (*it) );
		
	while( !Stack.empty() )
	{
		CTreeDecompositionCluster* currClust = Stack.back();
		Stack.pop_back();

		if (currClust->getVisited() == 1)
			return true;

		currClust->setVisited(1);

		treeDecompositionCluster_v::iterator it_ch = currClust->children().begin();
		for ( ; it_ch != currClust->children().end() ; it_ch ++)
			Stack.push_back(*it_ch);		
	}

	return false;
}

void CTreeDecomposition::setClusters(treeDecompositionCluster_v &clusters)
{
	copy(clusters.begin(), clusters.end(), inserter(m_clusters, m_clusters.begin()));
}

void CTreeDecomposition::setSeparators(treeDecompositionSeparator_v &separators)
{
	copy(separators.begin(), separators.end(), inserter(m_separators, m_separators.begin()));
}	