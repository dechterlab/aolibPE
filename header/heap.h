// heap.h -- Heap data structure.

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

#pragma warning (disable : 4786)


#ifndef _LIB_HEAP_H
#define _LIB_HEAP_H

#include "defs.h"

class CHeapElement
{
public:
	void* m_object;
	double m_score;


	CHeapElement() : m_object(NULL), m_score(0) {};
	CHeapElement(void* obj, double score) : m_object(obj), m_score(score) {};
	~CHeapElement() {};
};

class CHeap
{
protected:
	vector<CHeapElement*> m_elements;
	int m_size;

public:
	int size()				{ return m_size; };
	bool empty()			{ return (m_size <= 0); };

public:
	void push(void* obj, double sc);
	void updateScore(void* obj, double sc);

	CHeapElement* top();
	void pop();

protected:
	void heapify(int i);
	void propagateUp(int i);
	int parent(int i);
	int left(int i);
	int right(int i);
	
	int index(void* obj);
	double score(int i);
public:
	CHeap();
	CHeap(int n, void** objects, double* scores);
	CHeap(vector<void*>& objects, vector<double>& scores);
	~CHeap();

};

#endif	// _LIB_HEAP_H
