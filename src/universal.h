/* -*-C-*-

Copyright (C) 2015 Massachusetts Institute of Technology

universal.h is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

universal.h is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with universal.c; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301,
USA.

*/

/* written by Jack Wisdom in June and July, 2015
   David M. Hernandez helped test and improve it.
   Further edits by David M. Hernandez, summer 2020.
    If you use universal.c please cite Wisdom and Hernandez (2015).
    */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef _STATE_
#define _STATE_
typedef struct {
  double x, y, z, xd, yd, zd, xdd, ydd, zdd;
} State;
#endif

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#define SUCCESS 1
#define FAILURE 0

