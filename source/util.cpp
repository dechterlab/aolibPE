// util.cpp -- Utility functions (random generators).

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

#include "defs.h"

/* CJM used in Uniform() */

long Seed_1 = 7683, Seed_2 = 4673 ;  /* signed 32-bit integers */


/*
	given n returns an integer k, 0 <= k < n.
*/

long RandUniform(long n)
{
	double z;
	long k;

	k = Seed_1 / 53668;
	Seed_1 = 40014 * (Seed_1 - k * 53668) - k * 12211;
	if (Seed_1 < 0) Seed_1 = Seed_1 + 2147483563;

	k = Seed_2 / 52774;
	Seed_2 = 40692 * (Seed_2 - k * 52774) - k * 3791;
	if (Seed_2 < 0) Seed_2 = Seed_2 + 2147483399;

	z = (double) (Seed_1 - Seed_2);

	if (z < 1.0) z = z + 2147483562.0;

/*
	return(z * 4.65661E-10);  // (z & 0xF00) >> 16 
*/

	long random_value = (long) (((double) n)*(z * 4.65661E-10)) ;
	return random_value ;
}


/*
	This function returns a random real number between 0 and 1.
*/

double RandUniformDouble(void)
{
	long z, k ;

    k = Seed_1 / 53668;
    Seed_1 = 40014 * (Seed_1 - k * 53668) - k * 12211;
    if (Seed_1 < 0) 
        Seed_1 = Seed_1 + 2147483563;

    k = Seed_2 / 52774;
    Seed_2 = 40692 * (Seed_2 - k * 52774) - k * 3791;
    if (Seed_2 < 0) 
        Seed_2 = Seed_2 + 2147483399;

    z = Seed_1 - Seed_2;

    if (z < 1) 
        z = z + 2147483562;

	double random_value = (z * 4.65661E-10) ;
	return random_value ;
}

/*********************************************************************
  3. This random number generator is from William H. Press, et al.,
     _Numerical Recipes in C_, Second Ed. with corrections (1994), 
     p. 282.  This excellent book is available through the
     WWW at http://nr.harvard.edu/nr/bookc.html.
     The specific section concerning ran2, Section 7.1, is in
     http://cfatab.harvard.edu/nr/bookc/c7-1.ps
*********************************************************************/

#define IM1   2147483563
#define IM2   2147483399
#define AM    (1.0/IM1)
#define IMM1  (IM1-1)
#define IA1   40014
#define IA2   40692
#define IQ1   53668
#define IQ2   52774
#define IR1   12211
#define IR2   3791
#define NTAB  32
#define NDIV  (1+IMM1/NTAB)
#define EPS   1.2e-7
#define RNMX  (1.0 - EPS)

/* ran2() - Return a random floating point value between 0.0 and
   1.0 exclusive.  If idum is negative, a new series starts (and
   idum is made positive so that subsequent calls using an unchanged
   idum will continue in the same sequence). */

double ran2(long *idum)
{
	int j ;
	long k ;
	static long idum2 = 123456789 ;
	static long iy = 0 ;
	static long iv[NTAB] ;
	float temp ;

	if (*idum <= 0) {                             /* initialize */
		if (-(*idum) < 1)                           /* prevent idum == 0 */
			*idum = 1 ;
		else
			*idum = -(*idum) ;                         /* make idum positive */
		idum2 = (*idum) ;
		for (j = NTAB + 7; j >= 0; j--) {           /* load the shuffle table */
			k = (*idum) / IQ1 ;
			*idum = IA1 * (*idum - k*IQ1) - k*IR1 ;
			if (*idum < 0)
				*idum += IM1 ;
			if (j < NTAB)
				iv[j] = *idum ;
			}
		iy = iv[0] ;
		}
      
	k = (*idum) / IQ1 ;
	*idum = IA1 * (*idum - k*IQ1) - k*IR1 ;
	if (*idum < 0)
		*idum += IM1 ;
	k = idum2/IQ2 ;
	idum2 = IA2 * (idum2 - k*IQ2) - k*IR2 ;
	if (idum2 < 0)
		idum2 += IM2 ;
	j = iy / NDIV ;
	iy = iv[j] - idum2 ;
	iv[j] = *idum ;
	if (iy < 1)
		iy += IMM1 ;
	if ((temp = AM * iy) > RNMX)
		return RNMX ;                                /* avoid endpoint */
	else
		return temp ;
}



/* Generator of normally distributed deviate
   with zero mean and unit variance using ran2(idum)
   as the source of uniform deviates
   (from "Numerical recipes in C", page 289)
*/

static long idum = -1 ;

float gasdev(void)
{
	static int iset = 0 ;
	static float gset ;
	float fac, rsq, v1, v2 ;

	if (0 == iset) {          /* We don't have an extra deviate handy, so */
		do {
			v1 = 2.0*ran2(&idum) - 1.0 ;  /* pick two uniform numbers in the square */
			v2 = 2.0*ran2(&idum) - 1.0 ;  /* extending from -1 to +1 in each */
			rsq = v1*v1 + v2*v2 ;      /* direction, see if they are in the unit */
			} while (rsq >= 1.0 || rsq == 0.0) ;    /* circle, and if they are not, */
		/* try again. */
		fac = sqrt(-2.0*log(rsq)/rsq) ;
		/* Now make the Box-Muller transformation to get two normal deviates. */
		/* Return one and save the other for the next time. */
		gset = v1*fac ;
		iset = 1 ;                   /* set flag */
		return v2*fac ;
		}
	else {                      /* we have an extra deviate handy, so */
		iset = 0 ;                   /* unset the flag */
		return gset ;                /* and return it. */
		}
}
