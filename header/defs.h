/// defs.h -- Global definitions.

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

#ifndef _MPELIB_DEFS_H
#define _MPELIB_DEFS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <memory.h>
#include <sys/types.h>
#include <sys/timeb.h>

// STL kernel
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <functional>
#include <algorithm>
#include <deque>
#include <list>
#include <queue>
#include <set>
#include <stack>

using namespace std;

//==============================


//#define DEBUG_TIME	
//#define WIN32
#define LINUX



//==============================



#define UNKNOWN			-1
#define NOID			-1
#define INFINITY		-1

#define MINUS_INFINITY	-INT_MAX

// Execution output codes.
#define S_SUCCESS		0
#define S_FAILURE		1
#define S_TIMEOUT		2
#define S_BACKTRACK		3
#define S_END			4
#define S_SENDALL		5
#define S_EMPTY			6

#define E_SUCCESS		0
#define E_FAILURE		1

#define SV_SUCCESS		0
#define SV_NOVALUE		1
#define SV_FAILURE		2

#define MSG_LAMBDA		0
#define MSG_PI			1
#define MSG_NONE		2

#define MC_HARD			0
#define MC_SOFT			1
#define MC_BOTH			2

#define MB_HARD			0
#define MB_SOFT			1
#define MB_BOTH			2

#define SL_LOOKAHEAD	0
#define SL_FORWARDCHECK	1
#define SL_UNKNOWN		2

#define VO_DEFAULT		0
#define VO_MINDEGREE	1
#define VO_MINWIDTH		2
#define VO_MINFILL		3
#define VO_MAXCARD		4
#define VO_RANDOMIZED	5
#define VO_MINWEIGHT	6
#define VO_MAXDOMAINMINFILL	7
#define VO_DETERMINISM	8
#define VO_TOPOLOGICAL	9
#define VO_WEIGHTEDMINFILL	10

#define STATIC			0
#define DYNAMIC			1

#define SA_BACKTRACK	0
#define SA_BACKJUMP		1
#define SA_LOOKAHEAD	2
#define SA_FORWARDCHECK	3

#define GRAPH_RANDOM		0		// Adnan Darwiche
#define GRAPH_RANDOM_FIXED	1		// Rina Dechter
#define GRAPH_GRID			2
#define GRAPH_CODING		3

#define MBE_SIMPLE		0
#define MBE_AUGMENT		1
#define MBE_PARTIAL		2

// algorithms
#define PRUNING_NO							0
#define PRUNING_CONSTRAINTS_ONLY			1
#define PRUNING_FORWARD_CHECKING			2
#define PRUNING_FC_PROJECTION				3 // this checks all constraints by projecting on the current assignment
#define PRUNING_FULL_LOOK_AHEAD				4
#define AUXILIARY							5
#define BE_AUXILIARY						6 // bucket elimination on the auxiliary network
#define OR_PRUNING_CONSTRAINTS_ONLY			7
#define OR_PRUNING_FORWARD_CHECKING			8
#define OR_PRUNING_FC_PROJECTION			9 // this checks all constraints by projecting on the current assignment
#define OR_PRUNING_NO						10
#define AO_LINEAR_ADVANCED_ID				11
#define AO_LINEAR_BE_ID						12
#define AO_ADVANCED_ADVANCED_ID				13
#define AO_ADVANCED_BE_ID					14


// define the algorithms
#define AO_NO_PRUNING						1
#define	AO_CONSTRAINTS						2
#define AO_FC								4
#define	AO_FC_PROJ							8
//-----------------------------				
#define	AO_AUX								16	// this only prunes when it finds zeros, constraints are incorporated into CPTs
#define BE_AUX								32
//-----------------------------				
#define OR_CONSTRAINTS						64	// checking constraints only
#define OR_FC								128	
#define OR_FC_PROJ							256
#define OR_NO_PRUNING						512
//-----------------------------
// w-cutset algorithms
#define AO_LINEAR_ADVANCED					1024
#define AO_LINEAR_BE						2048
#define AO_ADVANCED_ADVANCED				4096
#define AO_ADVANCED_BE						8192

#define NUMBER_OF_ALGORITHMS	17
#define NUMBER_OF_MEASURES		20


#define INP_SIMPLE			1
#define INP_ORDER			2

#define MSB					1
#define LSB					2

// for building andor cutsets; when breaking the tree decomposition
#define INITIAL_STATE		-10
#define	EMPTY_CLUSTER		-5
#define PARENT_SEP_BURNT	0



// Min, Max macros
#define min(a,b)  (((a) < (b)) ? (a) : (b))
#define max(a,b)  (((a) > (b)) ? (a) : (b))

// Divide and either round up or down 
#define divRoundDown(n, s)  ((n) / (s))
#define divRoundUp(n, s)    (((n) / (s)) + ((((n) % (s)) > 0) ? 1 : 0))

class CFunction;
class CMiniBucket;
class CCluster;
class CMiniCluster;
class CMessage;
class CLegalTreeNode;
class CAONode;
class CConstraintNode;
class CTreeDecomposition;
class CTreeDecompositionCluster;		
class CTreeDecompositionSeparator;		


typedef map<unsigned long, CFunction*, less<unsigned long> > function_map;
typedef multimap<int, int, less<int> > variable_multimap;
typedef vector<CFunction*> function_v;
typedef vector<CMiniBucket*> minibucket_v;
typedef vector<int> variable_v;
typedef set<int, less<int> > variable_s;
typedef vector<CCluster*> cluster_v;
typedef vector<CMiniCluster*> minicluster_v;
typedef vector<CMessage*> message_v;
typedef vector<CLegalTreeNode*> legaltreenode_v;
typedef vector<CAONode*> aonode_v;
typedef vector<CConstraintNode*> constraintnode_v;
typedef vector<CTreeDecompositionCluster*> treeDecompositionCluster_v;
typedef vector<CTreeDecompositionSeparator*> treeDecompositionSeparator_v;
typedef vector<CTreeDecomposition*> treeDecomposition_v;



typedef enum {LAMBDA, PI} MTYPE;

// Defines the "search level" structure.
typedef struct tagSLINFO
{
	int level;						// Search tree level.
	int status;						// Search level status.
	
	int* domainSize;				// Current domain sizes.
	double* domainCost;				// Current domain costs (per value).
	bool* domainState;				// Current domain states (per value).

	struct tagSLINFO* prev;			// Previous link.

	void destroy()
	{
		delete[] domainSize;
		delete[] domainCost;
		delete[] domainState;
	}

} SLINFO;

typedef struct tagVALUECOST
{
	int value;
	double cost;

	tagVALUECOST(int v, double c) : value(v), cost(c) {};

} VALUECOST;

typedef struct tagVARIABLE
{
	int index;
	int name;
	int domain;

	tagVARIABLE(int i, int n, int d) : index(i) , name(n) , domain(d) {};

} VARIABLE;


inline bool vc_asc(VALUECOST* a, VALUECOST* b)
{
	return (a->cost < b->cost);
}

inline bool vc_desc(VALUECOST* a, VALUECOST* b)
{
	return (a->cost > b->cost);
}

long RandUniform(long);
double RandUniformDouble();

typedef map<int, VARIABLE* > variable_m;


#endif	// _DFSLIB_DEFS_H

// Local Variables:
// mode: C++
// End:
