// timers.cpp -- Timer management routines.

/*
 * Copyright (C) 2005-2006 Radu Marinescu
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


#include "timers.h"

/* --------------------------------------------------------------------
// Timer management functions
// -------------------------------------------------------------------- */

// WINDOWS IMPLEMENTATION
#ifdef WIN32

//double cpuTime()
//{
//	clock_t a = clock();
//	double tick = (double)a;
//
//	return (tick / (double)CLOCKS_PER_SEC);
//}

double cpuTime()
{      
   __int64 cycles;
   _asm rdtsc; // won't work on 486 or below - only pentium or above
   _asm lea ebx,cycles;
   _asm mov [ebx],eax;
   _asm mov [ebx+4],edx;
   return ((double)cycles) / 2258000000.0 ;
}

#endif

// LINUX IMPLEMENTATION
#ifdef LINUX

#include <unistd.h>
#include <sys/time.h>
#include <sys/times.h>

//// Return CPU time in seconds
//double cpuTime()
//{
//  static struct tms buf;
//
//  times(&buf);
//  return ((double) (buf.tms_utime+buf.tms_stime+buf.tms_cutime+buf.tms_cstime)) / ((double) sysconf(_SC_CLK_TCK));
//}


// Reads the TimeStampCounter on the Pentium and returns a 64 bit result.
//long long cpuTime()
//{  
//	long long t;  
//	asm (".byte 0x0f,0x31"  
//		: "=A"(t)                      
//		: : "%eax", "%edx"); 
//	return t;   
//}

double cpuTime()
{
	double t;
	struct timeval tv;
	struct timezone tz;
	gettimeofday(&tv, &tz);

	t = ((double)	tv.tv_sec) + ((double) tv.tv_usec)*0.000001;
	return t;
}

#endif
