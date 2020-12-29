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

/*
 * NOTE: This is an interal header file.
 * You should not attempt to use it directly.
 */

#ifndef AOLIB_BITVECTOR_H
#define AOLIB_BITVECTOR_H

#include "defs.h"

#define BITSINBYTE 	8
#define BITSINWORD 	32

// Divide and either round up or down 
#define divRoundDown(n, s)  ((n) / (s))
#define divRoundUp(n, s)    (((n) / (s)) + ((((n) % (s)) > 0) ? 1 : 0))


////////////////////////////////////////////////////////////////
// CBitVector class interface.

class CBitVector
{
public:
	CBitVector();
    CBitVector(int n);						// Initialize a bitvector, with "n" bits.
											// initially, all bits are cleared.
	CBitVector::CBitVector(int n, bool val); // initialize "n" bits to val

    virtual ~CBitVector();					// De-allocate bitvector.
    
    void setBit(int i);   					// Set the "i-th" bit.
    void freeBit(int i);  					// Clear the "i-th" bit.
    bool testBit(int i);   					// Is the "i-th" bit set?
    
	int size();								// Returns number of bits.
	void reset();

	void dump(ostream& os);					// Dump the content to output stream.

private:
    int m_numBits;							// number of bits in the bitvector
    int m_numWords;							// number of words of bitvector storage
											// (rounded up if m_numBits is not a
											//  multiple of the number of bits in
											//  a word)

    unsigned int *storage;					// bit storage

public:

	unsigned int* getBitVector()			{ return storage; };
	int			  getNumWords()				{ return m_numWords; };

};

#endif // AOLIB_BITVECTOR_H

// Local Variables:
// mode: C++
// End:
