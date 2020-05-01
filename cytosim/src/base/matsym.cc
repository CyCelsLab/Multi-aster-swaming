//RCS: $Id: matsym.cc,v 2.2 2005/03/01 13:56:57 nedelec Exp $
//======================================================================
//=               mat.cc  -  a full symmetrix matrix                   =
//======================================================================

#include "matsym.h"
#include "main.h"
#include "cblas.h"


//----------------------------------------------------------------------
void MatrixSymmetric::allocate()
{
  size      = 0;
  allocated = 0;
  val       = 0;
  do_not_delete_array = 0;
}

//----------------------------------------------------------------------
void MatrixSymmetric::allocate( int sz )
{
  size = sz;
  if ( size > allocated ) {
    allocated = size;
    if ( val ) delete[] val;
    val = new real[ size * size ];
    
    if ( val == 0 ) {
      fprintf(stderr, "MatrixSymmetric::allocate(int):: memory allocation failed\n");
      exit(1);
    }
  }
}

//----------------------------------------------------------------------
void MatrixSymmetric::deallocate()
{
  if ( do_not_delete_array ) return;
  if ( val ) delete[] val;
  allocated = 0;
  val = 0;
}

//----------------------------------------------------------------------
void MatrixSymmetric::setToZero()
{
  for( int ii = 0; ii < size * size; ++ii )
    val[ ii ] = 0;
}

//----------------------------------------------------------------------
void MatrixSymmetric::scale( real a )
{
  for( int ii = 0; ii < size * size; ++ii )
    val[ ii ] *= a;
}

//----------------------------------------------------------------------
real & MatrixSymmetric::operator()(int x, int y)
{
  assert( ( x >= 0 ) && ( x < size ) );
  assert( ( y >= 0 ) && ( y < size ) );
  if ( x < y )
    return val[ x + size * y ];
  else
    return val[ y + size * x ];
}

//----------------------------------------------------------------------
real * MatrixSymmetric::addr(int x, int y) const
{
  assert( ( x >= 0 ) && ( x < size ) );
  assert( ( y >= 0 ) && ( y < size ) );
  if ( x < y )
    return &val[ x + size * y ];
  else
    return &val[ y + size * x ];
}

//----------------------------------------------------------------------
bool MatrixSymmetric::nonZero() const
{
  return true;
}

//----------------------------------------------------------------------
int MatrixSymmetric::nbNonZeroElements() const
{
  return size * size;
}

//----------------------------------------------------------------------
char * MatrixSymmetric::shortDescription() const
{
  static char msg[16];
  snprintf(msg, sizeof(msg), " fullSym");
  return msg;
}

//========================================================================
//========================================================================

void MatrixSymmetric::multiplyVectorAdd( const real * X, real * Y ) const 
{
  blas_xsymv( 'U', size, 1.0, val, size, X, 1, 1.0, Y, 1 );
}

//----------------------------------------------------------------------
void MatrixSymmetric::multiplyVectorAddIsotropic( const real * X, real * Y ) const
{
  for( int d = 0; d < DIM; ++d )
    blas_xsymv( 'U', size, 1.0, val, size, X+d, DIM, 1.0, Y+d, DIM );  
}

