//RCS: $Id: matrix.cc,v 2.1 2005/01/10 15:43:17 foethke Exp $

#include "matrix.h"
#include "assert_macro.h"
#include "cblas.h"
#include "main.h"

//----------------------------------------------------------------------
real Matrix::value(const int x, const int y) const 
{
  real * v = addr( x, y );
  if ( v == 0 ) 
    return 0;
  else
    return *v;
}

//----------------------------------------------------------------------
real Matrix::maxNorm() const 
{
  real result = 0;
  for( int ii = 0; ii < size; ++ii )
    for( int jj = 0; jj < size; ++jj ) {
      real * v = addr( ii, jj );
      if ( v  &&  ( *v > result ) )
        result = *v;
    }
  return result;
}

//----------------------------------------------------------------------
bool Matrix::nonZero() const
{
  for( int ii = 0; ii < size; ++ii )
    for( int jj = 0; jj < size; ++jj )
      if ( 0 != value( ii, jj ) )
        return true;
  return false;
}

//----------------------------------------------------------------------
int  Matrix::nbNonZeroElements() const
{
  int result = 0;
  for( int ii = 0; ii < size; ++ii )
    for( int jj = 0; jj < size; ++jj )
      result += ( 0 != value( ii, jj ) );
  return result;
}

//----------------------------------------------------------------------
// M <- the block contained in [x, x+sx, y, y+sy]
void Matrix::copyBlock(real* M, int x, int sx, int y, int sy) const
{
  assert( x + sx <= size );
  assert( y + sy <= size );
  for(int ii = 0; ii < sx; ++ii)
    for(int jj = 0; jj < sy; ++jj)
      M[ ii + sx * jj ] = value( x + ii, y + jj );
}

//----------------------------------------------------------------------
// add the diagonal block ( x, x, x+sx, x+sx ) to M
void Matrix::addDiagonalBlock(real* M, const int x, const int sx) const
{
  for(int ii = 0; ii < sx; ++ii)
    for(int jj = 0; jj < sx; ++jj)
      M[ ii + sx * jj ] += value( x + ii, x + jj );  
}

//----------------------------------------------------------------------
// add the upper triagular block ( x, x, x+sx, x+sx ) to M
void Matrix::addTriangularBlock(real* M, const int x, const int sx) const
{
  for(int ii = 0; ii < sx; ++ii)
    for(int jj = ii; jj < sx; ++jj)
      M[ ii + sx * jj ] += value( x + ii, x + jj );    
}


//----------------------------------------------------------------------
// multiplication of a vector: Y = M * X, dim(X) = dim(M)
void Matrix::multiplyVector( const real * X, real * Y ) const 
{
  blas_xzero( size, Y, 1 ); 
  multiplyVectorAddIsotropic( X, Y );
}

//----------------------------------------------------------------------
void Matrix::multiplyVectorIsotropic( const real * X, real * Y ) const 
{
  blas_xzero( DIM*size, Y, 1 ); 
  multiplyVectorAddIsotropic( X, Y );
}


//----------------------------------------------------------------------
void Matrix::printFull(FILE * f) const
{
  //printf("%i %i\n", size, size);
  for(int ii = 0; ii < size; ++ii ) {
    for(int jj = 0; jj < size; ++jj )
      fprintf(f, "%6.2f ", value(ii, jj) );
    fprintf(f, "\n");
  }
}

//----------------------------------------------------------------------
void Matrix::printSparse(FILE * f) const
{
  for (int ii = 0; ii < size; ++ii )
    for(int jj = 0; jj < size; ++jj ) 
      if ( addr( ii, jj ) )
        fprintf(f, "%i %i %16.8e\n", ii, jj, *addr(ii, jj) );
}



