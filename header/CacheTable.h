// CacheTable.h: Cache Tables.
//
//////////////////////////////////////////////////////////////////////

/*
 * Copyright (C) 2004 Robert Mateescu
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

#ifndef _MPELIB_CACHE_H
#define _MPELIB_CACHE_H

#include "function.h"


class CCacheTable : public CFunction  
{

protected:

	double* m_table;		// Cache table.
	int m_tableSize;		// Cache table size;

public:
	CCacheTable();
	CCacheTable(int id, CTYPE type, int argc, int* argv);


	virtual ~CCacheTable();

public:
	void init();
	void purge(); // erases all values in cache and sets them to -1

	double getValueAt(int addr);
	void setValueAt(int addr, double val);

	int getCurrentValue(void* val);
	int setCurrentValue(void* val);


	void destroy();
	void print(bool scope, bool table);

	int getTableSize() { return m_tableSize; };

	void substitute(int var, int val, CFunction*& fun);
	bool verify();
};

#endif // _MPELIB_CACHE_H
