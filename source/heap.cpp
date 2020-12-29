// heap.cpp -- Heap data structure.

/*
 * Copyright (C) 2005 Radu Marinescu
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

#include "heap.h"


CHeap::CHeap()
{
	m_size = 0;
}

CHeap::CHeap(int n, void** objects, double* scores)
{
	int i;

	m_elements.resize(n);
	for (i = 0; i < n; ++i)
	{
		m_elements[i]->m_object = objects[i];
		m_elements[i]->m_score  = scores[i];
	}

	m_size = n;
	for (i = n - 1; i >= 0; --i)
		heapify(i);
}

CHeap::CHeap(vector<void*>& objects, vector<double>& scores)
{
	int i, n = (int)objects.size();

	m_elements.resize(n);
	for (i = 0; i < n; ++i)
	{
		m_elements[i] = new CHeapElement(objects[i], scores[i]);
	}

	m_size = n;
	for (i = n - 1; i >= 0; --i)
		heapify(i);
}

CHeap::~CHeap()
{
	for (int i = 0; i < (int)m_elements.size(); ++i)
		delete m_elements[i];
	m_elements.clear();
}

int CHeap::index(void* obj)
{
	for (int i = 0; i < (int)m_elements.size(); ++i)
		if (obj == m_elements[i]->m_object)
			return i;

	return UNKNOWN;
}

// Push an new element in the heap.
void CHeap::push(void* obj, double sc)
{
	CHeapElement* he = new CHeapElement(obj, sc);
	int i = m_size;
	m_elements.push_back(he);
	m_size++;
	propagateUp(i);
}

// Update the priority of an element in the heap.
void CHeap::updateScore(void* obj, double sc)
{
	int i = index(obj);
	assert((i != UNKNOWN) && (i < m_size));
	CHeapElement* he = m_elements[i];
	bool increased = (sc > he->m_score);

	he->m_score = sc;
	if (increased)
		propagateUp(i);
	else
		heapify(i);
}

// Extracts the element with highest priority.
CHeapElement* CHeap::top()
{
	return m_elements[0];
}

// Remove the element from top of the heap.
void CHeap::pop()
{
	m_size--;
	swap(m_elements[0], m_elements[m_size]);
	m_elements.erase((m_elements.begin() + m_size));

	heapify(0);
}

// Arrange the elements in the list as heap.
void CHeap::heapify(int i) 
{
	int bestind;
	if (left(i) >= m_size) 
	{
		return;
	} else if (score(left(i)) > score(i)) 
	{
		bestind = left(i);
	} else 
	{
		bestind = i;
	}
	
	if ((right(i) < m_size) && (score(right(i)) > score(bestind))) 
	{
		bestind = right(i);
	}
	
	if (bestind != i) 
	{
		swap(m_elements[bestind], m_elements[i]);
		heapify(bestind);
	}
}

void CHeap::propagateUp(int i) 
{
	double sc = score(i);
	while (sc > score(parent(i))) 
	{
		swap(m_elements[i], m_elements[parent(i)]);
		i = parent(i);
	}
}

double CHeap::score(int i) 
{
	CHeapElement* he = (CHeapElement*) m_elements[i];
	
	return he->m_score;
}

int CHeap::parent(int i)
{
	return (i - 1) / 2;
}

int CHeap::left(int i)
{
	return 2 * i + 1;
}

int CHeap::right(int i)
{
	return 2 * i + 2;
}
