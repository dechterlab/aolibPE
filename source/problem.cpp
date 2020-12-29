// Problem.cpp -- The super class of all the CSPs.

/*
 * Copyright (C) 2003 Radu Marinescu
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

#include "problem.h"
#include "graph.h"

CProblem::CProblem()
{
	m_N = 0;
	m_K = 0;
	m_C = 0;

	m_ordering = NULL;
	m_position = NULL;
	m_assignment = NULL;
	m_evidence = NULL;
	m_domains = NULL;
	m_timePassed = 0;
	m_adjacencies = NULL;
	m_backupAssignment = NULL;
	m_backupEvidence = NULL;
}

CProblem::CProblem(int n, int k, int p)
	: m_N(n), m_K(k), m_C(p)
{
	// m_C = 0;

	m_ordering = NULL;
	m_position = NULL;
	m_assignment = NULL;
	m_evidence = NULL;
	m_domains = NULL;
	m_timePassed = 0;
	m_adjacencies = NULL;
	m_backupAssignment = NULL;
	m_backupEvidence = NULL;
}

CProblem::~CProblem()
{
	// Deallocate memory.
	if (m_domains)
		delete[] m_domains;

	if (m_ordering)
		delete[] m_ordering;

	if (m_position)
		delete[] m_position;

	if (m_assignment)
		delete[] m_assignment;

	if (m_evidence)
		delete[] m_evidence;

	if (m_backupAssignment)
		delete[] m_backupAssignment;

	if (m_backupEvidence)
		delete[] m_backupEvidence;


	function_map::iterator it = m_functions.begin();
	for ( ; it != m_functions.end(); ++it)
		delete (*it).second;
	m_functions.clear();

	// Clear the adjacency matrix
	if (m_adjacencies)
	{
		for (int i = 0; i < m_N; ++i)
			delete m_adjacencies[i];
		delete[] m_adjacencies;
	}
}

int CProblem::getStaticDomainSize(int var)
{
	// Safety checks.
	assert(m_domains != NULL);
	assert(var < m_N);

	return m_domains[var];
}

CFunction* CProblem::getFunctionByID(int id)
{
	function_map::iterator it = m_functions.find(id);

	return (it != m_functions.end()) ?
		(*it).second : NULL;
}

void CProblem::addFunction(CFunction* c, bool w_id)
{
	unsigned long id = (w_id) ? ++m_C : c->getID();
	if (w_id) c->setID(id);

	m_functions.insert(make_pair(id, c));
}

void CProblem::setValue(int var, int val)
{
	// Safety checks.
	assert(m_assignment != NULL);
	assert((var >= 0) && (var < m_N));

	m_assignment[var] = val;
}

void CProblem::resetValue(int var)
{
	// Safety checks.
	assert(m_assignment != NULL);
	assert((var >= 0) && (var < m_N));

	m_assignment[var] = -1;
}

int CProblem::getValue(int var)
{
	// Safety checks.
	assert(m_assignment != NULL);
	assert((var >= 0) && (var < m_N));

	return m_assignment[var];
}

bool CProblem::isEvidence(int var)
{
	// Safety checks.
	assert(m_assignment != NULL);
	assert((var >= 0) && (var < m_N));

	return m_evidence[var];
}

bool CProblem::isAssigned(int var)
{
	// Safety checks.
	assert(m_assignment != NULL);
	assert((var >= 0) && (var < m_N));

	return (m_assignment[var] >= 0);
}

void CProblem::setEvidence(int var)
{
	// Safety checks.
	assert(m_assignment != NULL);
	assert((var >= 0) && (var < m_N));

	m_evidence[var] = true;
}

void CProblem::setEvidence(int var, int val)
{
	// Safety checks.
	assert(m_assignment != NULL);
	assert((var >= 0) && (var < m_N));

	m_evidence[var] = true;
	m_assignment[var] = val;
}

void CProblem::resetEvidence(int var)
{
	// Safety checks.
	assert(m_assignment != NULL);
	assert((var >= 0) && (var < m_N));

	m_evidence[var] = false;
}

//void CProblem::startClock(time_t *time_counter)
//{
//	// get and store the current timer value.
//	// note, this is really the time elapsed since some reference point.
//	struct _timeb timebuffer ;
//	_ftime(&timebuffer) ;
//	if (time_counter) 
//	{
//		*time_counter = 1000*timebuffer.time + timebuffer.millitm ;
//	}
//	else 
//	{
//		m_time = 1000*timebuffer.time + timebuffer.millitm ;
//	}
//}
//
//unsigned long CProblem::stopClock(int add_to_passed_time, time_t *time_counter)
//{
//	// get the current timer value.
//	// note, this is really the time elapsed since some reference point.
//	struct _timeb timebuffer ;
//	_ftime(&timebuffer) ;
//	time_t endTime = 1000*timebuffer.time + timebuffer.millitm ;
//
//	// compute the time passed since the clock was started
//	unsigned long t ;
//	if (time_counter) 
//	{
//		t = endTime - (*time_counter) ;
//	}
//	else 
//	{
//		t = endTime - m_time ;
//	}
//
//	if (add_to_passed_time) m_timePassed += t ;
//
//	return t ;
//}

// This function returns the highest position / variable in "argv"
// according to the current variable ordering.
void CProblem::getHighestNode(int argc, int* argv, int& highestVar, int& highestPos)
{
	// Safety checks.
	assert(m_ordering != NULL);
	assert(m_position != NULL);
	assert(argv != NULL);
	assert(argc > 0);

	int hPos = -1, hVar = -1;
	for (int i = 0; i < argc; ++i)
	{
		int pos = m_position[argv[i]];
		if (pos > hPos)
		{
			hPos = pos;
			hVar = argv[i];
		}
	}

	highestPos = hPos;
	highestVar = hVar;
}

// This function returns the number of constraints (hard)
// a variable is involved in.
int CProblem::getVariableConflicts(int var)
{
	assert(m_adjacencies != NULL);
	assert((var >= 0) && (var < m_N));

	return m_adjacencies[var]->size();
}

void CProblem::backupAssignment()
{
	// Safety checks.
	assert(m_assignment != NULL);
	assert(m_evidence != NULL);
	assert(m_backupAssignment != NULL);
	assert(m_backupEvidence != NULL);
	assert(m_N > 0);

	memcpy(m_backupAssignment, m_assignment, m_N * sizeof(int));
	memcpy(m_backupEvidence, m_evidence, m_N * sizeof(bool));
}

void CProblem::restoreAssignment()
{
	// Safety checks.
	assert(m_assignment != NULL);
	assert(m_evidence != NULL);
	assert(m_backupAssignment != NULL);
	assert(m_backupEvidence != NULL);
	assert(m_N > 0);

	memcpy(m_assignment, m_backupAssignment, m_N * sizeof(int));
	memcpy(m_evidence, m_backupEvidence, m_N * sizeof(bool));
}
