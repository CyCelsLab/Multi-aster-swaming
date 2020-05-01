//RCS: $Id: smath.cc,v 2.14 2005/03/30 12:29:00 nedelec Exp $

#include "smath.h"
#include <cstring>

//--------------------------------------------------------
//returns 1 if w[] is inside the
//ellipse of axis   (1,0,0), (0,1,0), (0,0,1)
//and of half width  d[0],    d[1],    d[2]
int insideEllipse(const real * w, const real * d)
{
  real a, dis = 0;
  for( int ii = 0; ii < 3; ++ii )
    if ( d[ii] > 0 ) {
      a = w[ii] / d[ii];
      dis += a * a;
    } else if ( w[ii] ) return 0;
  return ( dis < 1.0 );
}

//--------------------------------------------------------
//returns 1 if the segments [a,b] and [c,d] are crossing
int cross( const real a[], const real b[], const real c[], const real d[],
           real * x, real * y)
{
  real abx = b[0] - a[0];
  real aby = b[1] - a[1];
  real cdx = d[0] - c[0];
  real cdy = d[1] - c[1];

  real det = abx * cdy - aby * cdx;

  if ( fabs(det) == 1e-6 ) return -1;

  *x = ( ( c[0] - a[0] ) * cdy - ( c[1] - a[1] ) * cdx ) / det;
  *y = ( ( c[0] - a[0] ) * aby - ( c[1] - a[1] ) * abx ) / det;

  if ( ( *x < 0 ) || ( *x > 1.0 ) ) return 0;
  if ( ( *y < 0 ) || ( *y > 1.0 ) ) return 0;
  return 1;
}

//--------------------------------------------------------
void copyBytes( void * dest, const void * src, const int cnt)
{
  for( int ii=0; ii < cnt; ++ii )
    ((char*)dest)[ii] = ((char*)src)[ii];
}

//--------------------------------------------------------
void printBits(const void * v, const int size)
{
  for(int ii=0; ii < size; ++ii ) {
    char c = ((char*)v)[ii];
    for(int jj=7; jj >=0; --jj )
      printf( "%d", ( c >> jj ) & 1 );
    printf(".");
  }
  printf("\n");
}

//--------------------------------------------------------
/*int isnan_array( const real* ap, const int size )
{
  int nbNan = 0;
#ifndef __APPLE__
  for( int ii=0; ii<size; ii++ ) {
#ifdef __INTEL_COMPILER
    nbNan += isnan( ap[ii] );
#elif defined __CYGWIN__
    nbNan += isnan( ap[ii] );
#else
    nbNan += __gnu_cxx::isnan( ap[ii] );
#endif
  }
#endif
  return nbNan;
}


int isinf_array( const real* ap, const int size )
{
  int nbInf = 0;
#ifndef __APPLE__
  for( int ii=0; ii<size; ii++ ) {
#ifdef __INTEL_COMPILER
    nbInf += isinf( ap[ii] );
#elif defined __CYGWIN__
    nbInf += isinf( ap[ii] );
#else
    nbInf += __gnu_cxx::isinf( ap[ii] );
#endif
  }
#endif
  return nbInf;
}


int isnotfinite_array( const real* ap, const int size )
{
  int nbNonFin = 0;
#ifndef __APPLE__
  for( int ii=0; ii<size; ii++ ) {
#ifdef __INTEL_COMPILER
    nbNonFin += ! isfinite( ap[ii] );
#elif defined __CYGWIN__
    nbNonFin += ! finite( ap[ii] );
#else
    nbNonFin += ! __gnu_cxx::isfinite( ap[ii] );
#endif
  }
#endif
  return nbNonFin;
}
*/
//==========================================================================
//                    MATRIX & VECTOR  BLAS-STYLED DUMP
//==========================================================================

void vectPrint( int m, const real * X, FILE * f, int digits )
{
  char digitStr[10];

  if ( X == 0 ) { fprintf(f, "void\n"); return; }

  snprintf(digitStr, sizeof(digitStr), "%%%d.%df ", digits+3, digits);
  for(int ii = 0; ii < m; ++ii)
    if ( X[ ii ] == 0 )
      fprintf(f, "   .   ");
    else
      fprintf(f, digitStr, X[ ii ] );
  fprintf(f, "\n");
}

//--------------------------------------------------------
void vectPrint( int m, const real * X, int digits, FILE * f )
{
  vectPrint(m, X, f, digits);
}

//--------------------------------------------------------
void vectDump( int m, const real * X, FILE * f, int digits )
{
  char digitStr[10];

  if ( X == 0 ) { fprintf(f, "void\n"); return; }

  snprintf(digitStr, sizeof(digitStr), "%%.%de\n", digits);
  for(int ii = 0; ii < m; ++ii)
    fprintf(f, digitStr, X[ ii ] );
}

//--------------------------------------------------------
void vectDump( int m, const real * X, int digits, FILE * f )
{
  vectDump(m, X, f, digits);
}

//--------------------------------------------------------
void matPrint( int m, int n, const real * X, FILE * f, int digits )
{
  char digitStr[10];

  if ( X == 0 ) { fprintf(f, "void\n"); return; }

  snprintf(digitStr, sizeof(digitStr), "%%%d.%df ", digits+3, digits);
  for(int ii = 0; ii < m; ++ii) {
    for(int jj = 0; jj < n; ++jj )
      if ( fabs( X[ ii + m * jj ] ) < 1e-4 )
	fprintf(f, "   .   ");
      else
	fprintf(f, digitStr, X[ ii + m * jj ] );
    fprintf(f, "\n");
  }
}

//--------------------------------------------------------
void matPrint( int m, int n, const real * X, int digits, FILE * f )
{
  matPrint(m, n, X, f, digits);
}

//--------------------------------------------------------
void matDump( int m, int n, const real * X, FILE * f, int digits )
{
  char digitStr[10];

  if ( X == 0 ) { fprintf(f, "void\n"); return; }

  snprintf(digitStr, sizeof(digitStr), "%%i %%i  %%.%de\n", digits);
  for(int ii = 0; ii < m; ++ii)
    for(int jj = 0; jj < n; ++jj )
      fprintf(f, digitStr, ii, jj, X[ ii + m * jj ] );
}

//--------------------------------------------------------
void matDump( int m, int n, const real * X, int digits, FILE * f )
{
  matDump(m, n, X, f, digits);
}

//--------------------------------------------------------
void matDumpOffset( int m, int n, const real * X, int off, FILE * f, int digits )
{
  char digitStr[10];

  if ( X == 0 ) { fprintf(f, "void\n"); return; }

  snprintf(digitStr, sizeof(digitStr), "%%i %%i  %%.%de\n", digits);
  for(int ii = 0; ii < m; ++ii)
    for(int jj = 0; jj < n; ++jj )
      fprintf(f, digitStr, ii+off, jj+off, X[ ii + m * jj ] );
}

//--------------------------------------------------------
void matDumpOffset( int m, int n, const real * X, int off, int digits, FILE * f )
{
  matDumpOffset(m, n, X, off, f, digits);
}

