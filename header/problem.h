// problem.h -- The super class of all the problems.

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


#ifndef _MPELIB_PROBLEM_H
#define _MPELIB_PROBLEM_H

#include "function.h"

///////////////////////////////////////////////////////////////////////
// CProblem abstract class interface.
class CProblem
{
protected:
	const string m_name;						// Generic name of the problem instance.

private:
	// Copy constructor not supported.
	CProblem(CProblem&);
	CProblem& operator=(const CProblem&);

protected:
	unsigned long		m_timePassed;			// Time used.
	time_t				m_time;					// Starting stop-watch time.

protected:
	int m_N;					// Number of variables.
	int m_K;					// Maximum limit of the domain size.
	int m_C;					// Number of constraints (CPTs)

	int* m_domains;				// Domain sizes for variables.

	int* m_ordering;			// Ordering of variables.
	int* m_position;			// Position of variables in the ordering.

	int* m_assignment;
	bool* m_evidence;
	
	int* m_backupAssignment;
	bool* m_backupEvidence;

protected:
    // The constraints in this problem.  We store them as pointers
    // we want to be able to have variables and constraints of
    // different types. We also prefer hash tables to lists and
    // vectors since they allow O(1) search & removals w/o any
    // need for address readjustment.

	function_map m_functions;		// Constraints of the problem.

	CFunctionList** m_adjacencies;	// The adjacency matrix.

public:
	function_map& functions()		{ return m_functions; };

	// Save current complete assignment/evidence.
	void backupAssignment();

	// Restore current complete assignment/evidence.
	void restoreAssignment();

	// Get the domain size for a variable.
	virtual int getStaticDomainSize(int var);

	// Get the "degree" of a variable.
	virtual int getVariableConflicts(int var);

	// Get functions associated with a variable.
	virtual function_v& getBucket(int var, bool dfsOrd = true) = 0;
	virtual variable_v& getTreeDescendants(int var) = 0;

public:
	// Get maximum domain size.
	virtual int getMaxDomainSize() { return m_K; };

	// Get a function from the repository.
	CFunction* getFunctionByID(int id);
	
	// Get variable ordering/position.
	virtual int* getOrdering(bool dfsOrd = false) = 0;
	virtual int* getPosition(bool dfsOrd = false) = 0;

	// Return (N, K, C) parameters.
	virtual int getN()		{ return m_N; };
	virtual int getK()		{ return m_K; };
	virtual int getC()		{ return m_C; };

	// Add a new function to the repository.
	virtual void addFunction(CFunction* c, bool w_id = true);

	// Single variable value maniputation routines.
	virtual void setValue(int var, int val);
	virtual int  getValue(int var);
	virtual void resetValue(int var);
	virtual bool isAssigned(int var);

	// Evidence flag manipulation routines.
	virtual bool isEvidence(int var);
	virtual void setEvidence(int var);
	virtual void setEvidence(int var, int val);
	virtual void resetEvidence(int var);
	
	// Return highest ordered argument from a list of arguments (argv).
	virtual void getHighestNode(int argc, int* argv, int& highestVar, int& highestPos);

	// Internal clock manipulation routines.
	/*void startClock(time_t *time_counter);
	unsigned long stopClock(int add_to_passed_time, time_t *time_counter);*/

public:
	// Constructor/Destructor.
	CProblem();
	CProblem(int n, int k, int p);
	virtual ~CProblem();
};


#endif	// _MPELIB_PROBLEM_H

// Local Variables:
// mode: C++
// End:
