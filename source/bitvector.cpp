// bitvector.h -- Bit vector class.

/*
 * Copyright (C) 2003-2006 Radu Marinescu
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

#include "bitvector.h"

CBitVector::CBitVector()
{
	m_numBits = 0;
	m_numWords = 0;
	storage = NULL;
}

CBitVector::CBitVector(int n) 
{ 
    m_numBits = n;
    m_numWords = divRoundUp(m_numBits, BITSINWORD);

    storage = new unsigned int[m_numWords];
	memset(storage, 0, m_numWords * sizeof(unsigned int));
}

CBitVector::CBitVector(int n, bool val) 
{ 
    m_numBits = n;
    m_numWords = divRoundUp(m_numBits, BITSINWORD);

    storage = new unsigned int[m_numWords];

	if(val)
	{
		for (int i=0; i<n; i++)
			setBit(i);
	}
	else
	{
		for (int i=0; i<n; i++)
			freeBit(i);
	}
}

CBitVector::~CBitVector()
{ 
	if (storage)
		delete storage;
}

void CBitVector::setBit(int i) 
{ 
    assert(i >= 0 && i < m_numBits);
    storage[i / BITSINWORD] |= 1 << (i % BITSINWORD);
}
    
void CBitVector::freeBit(int i) 
{
    assert(i >= 0 && i < m_numBits);
    storage[i / BITSINWORD] &= ~(1 << (i % BITSINWORD));
}

bool CBitVector::testBit(int i)
{
    assert(i >= 0 && i < m_numBits);
	return (storage[i / BITSINWORD] & (1 << (i % BITSINWORD))) ? true : false;
}

void CBitVector::reset()
{
	memset(storage, 0, m_numWords * sizeof(unsigned int));
}

int CBitVector::size()
{
	return m_numBits;
}

void CBitVector::dump(ostream& os)
{
	for (int i = 0; i < size(); ++i)
		os << " " << (testBit(i) ? "1" : "0");
}

