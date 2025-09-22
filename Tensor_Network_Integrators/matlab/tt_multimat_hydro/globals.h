/*****************************************************************************
*  Globals 
*  *
*  *  Global declarations common to all library code. Standard headers that are
*  *  always utilized are included here. Commonly used constants and Macro
*  *  functions are also declared.
*  *****************************************************************************/
#ifndef _Globals_
#define _Globals_

#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/time.h>
#include <time.h>

/* max num of dimension for which this tool will work. */

#define PI               3.14159265358979323846264338327

#define ABS(x)   (((x) >= 0) ? (x):(-(x)))  
#define SIGN(x)  (((x) >= 0) ? 1.0:(-1.0))
#define MAX(x,y) (((x) >= (y)) ? (x):(y))
#define MIN(x,y) (((x) <= (y)) ? (x):(y)) 

#ifndef SWAP
#define SWAP(a, b, t) {t h; h = (a); (a) = (b); (b) = h; }
#endif

#ifndef null
#define null NULL
#endif

#endif

