//RCS: $Id: matsparse.cc,v 2.2 2005/03/01 13:56:27 nedelec Exp $
//======================================================================
//=                           sparse.cc                                =
//======================================================================
#include <cstdio>

#include "matsparse.h"
#include "cblas.h"
#include "main.h"
#include "smath.h"

#define FAST_VECTMULT

#define VAL_NOT_USED -1
#define END_OF_COLUMN -2

//----------------------------------------------------------------------
void MatrixSparse::allocate()
{ 
  size      = 0;
  allocated = 0;
  Vcol      = 0;
  Vrow      = 0;
}


//----------------------------------------------------------------------
void MatrixSparse::allocate( const int sz )
{
  size = sz;
  if ( size > allocated ) {

    real ** Vcol_new = new real * [ size ];	  
    int  ** Vrow_new = new int * [ size ];	 
    
    if ((Vcol_new == 0) || (Vrow_new == 0)) {
      fprintf(stderr, "MatrixSparse::allocate(int):: memory allocation failed\n");
      exit(1);
    }
    
    int ii = 0;
    if ( Vcol ) {
      for( ; ii < allocated; ++ii ) {
        Vcol_new[ ii ] =  Vcol[ ii ];
        Vrow_new[ ii ] =  Vrow[ ii ];
      }
      delete[] Vcol;
      delete[] Vrow;
    }
    
    for( ; ii < size; ++ii ) {
      Vcol_new[ ii ] = 0;
      Vrow_new[ ii ] = 0;
    }
  
    Vcol = Vcol_new;
    Vrow = Vrow_new;
    allocated = size;
  }
}

//----------------------------------------------------------------------
void MatrixSparse::deallocate()
{
  if ( Vcol ) {
    for (int ii = 0; ii < allocated; ++ii )
      if ( Vcol[ ii ] ) {
        delete[] Vcol[ ii ];
        delete[] Vrow[ ii ];
      };
    delete[] Vcol;
    delete[] Vrow;
  }
}

//----------------------------------------------------------------------
void MatrixSparse::allocateColumn( const int jj, const int sz )
{
  assert((jj >= 0) && (jj < size));
  assert( sz > 0 );
  //printf("new S-COL %i %i\n", jj, sz );

  real * Vcol_new  = new real[ sz ];
  int  * Vrow_new  = new  int[ sz ];
  
  if ((Vcol_new == 0) || (Vrow_new == 0)) {
    fprintf(stderr, "MatrixSparse::allocateColumn(int):: memory allocation failed\n");
    exit(1);
  }
  
  int ii = 0;
  if ( Vcol[ jj ] ) {	  
    for( ; Vrow[ jj ][ ii ] != END_OF_COLUMN ; ++ii ) {
      Vcol_new[ ii ] =  Vcol[ jj ][ ii ];
      Vrow_new[ ii ] =  Vrow[ jj ][ ii ];
    }
    
    delete[] Vcol[ jj ];
    delete[] Vrow[ jj ];
  }
  for( ; ii < sz-1 ; ++ii ) 
    Vrow_new[ ii ] = VAL_NOT_USED;
  Vrow_new[ sz-1 ] = END_OF_COLUMN;
  
  Vcol[ jj ]  = Vcol_new;
  Vrow[ jj ]  = Vrow_new;
}


//----------------------------------------------------------------------
real & MatrixSparse::operator()( int x, int y)
  //allocate the position if necessary:
{
  assert( ( x >= 0 ) && ( x < size ) );
  assert( ( y >= 0 ) && ( y < size ) );
 
  if ( Vrow[ y ] ) {
    
    int ii = 0;
    for( ; Vrow[ y ][ ii ] >= 0; ++ii )
      if ( Vrow[ y ][ ii ] == x )
        return Vcol[ y ][ ii ];
    
    if( Vrow[ y ][ ii ] == END_OF_COLUMN )
      allocateColumn( y, ii + 1 + SPCHUNCK );
    assert( Vrow[ y ][ ii ] == VAL_NOT_USED );
    Vrow[ y ][ ii ] = x;
    Vcol[ y ][ ii ] = 0;
    //printf("allo. %3i %3i\n", x, y );
    return Vcol[ y ][ ii ];
  }
  
  allocateColumn( y, SPCHUNCK );
  //printf("allo. %3i %3i\n", nx, ny );
  assert( Vrow[ y ][ 0 ] == VAL_NOT_USED );

  //put the diagonal term first:
  Vrow[ y ][ 0 ] = y;
  Vcol[ y ][ 0 ] = 0;
  if ( x == y )
    return Vcol[ y ][ 0 ];
  
  Vrow[ y ][ 1 ] = x;
  Vcol[ y ][ 1 ] = 0;
  return Vcol[ y ][ 1 ];
}


//----------------------------------------------------------------------
real * MatrixSparse::addr( int x, int y) const
  //does not allocate the position:
{
  int * row = Vrow[ y ];
  if ( row )
    for( ; *row >= 0 ; ++row )
      if ( *row == x )
        return & Vcol[ y ][ row - Vrow[ y ] ];
  return 0;
}
  

//----------------------------------------------------------------------
void MatrixSparse::setToZero()
{
  for(int ii = 0; ii < size; ++ii ) 
    if ( Vrow[ ii ] ) 
      for(int jj = 0; Vrow[ ii ][ jj ] != END_OF_COLUMN; ++jj )
        Vrow[ ii ][ jj ] = VAL_NOT_USED;
}


//----------------------------------------------------------------------
void MatrixSparse::scale( real a )
{
  for(int ii = 0; ii < size; ++ii ) if ( Vrow[ ii ] ) 
    for(int jj = 0; Vrow[ ii ][ jj ] >= 0; ++jj )
      Vcol[ ii ][ jj ] *= a;
}

//----------------------------------------------------------------------
  // M <- M + the upper block contained in [x, x+sx, x, x+sx]
void MatrixSparse::addTriangularBlock(real* M, const int x, const int sx ) const
{
  assert( x + sx <= size );

  for(int jj = 0; jj < sx; ++jj ) {
    int * row = Vrow[ jj + x ];
    if ( row != 0 ) {
      real * col = Vcol[ jj + x ];
      for( ; *row >= 0 ; ++row, ++col) {
        int ii = *row - x;
        if ( ( ii >= 0 ) && ( ii < sx ) ) {
          assert( ii <= jj );	
          M[ ii + sx * jj ] += * col;
          //printf("Sp %4i %4i %.4f\n", ii, jj, a );
        }
      }
    }
  }
}

//----------------------------------------------------------------------
// M <- M + the upper block contained in [x, x+sx, x, x+sx]
void MatrixSparse::addDiagonalBlock(real* M, const int x, const int sx ) const
{
  assert( x + sx <= size );

  for(int jj = 0; jj < sx; ++jj ) {
    
    int *row = Vrow[ jj + x ];
    if ( row ) {
      real * col = Vcol[ jj + x ];
      for( ; *row >= 0 ; ++row, ++col) {
        int ii = *row - x;
        if (( ii >= 0 ) && ( ii < sx )) {
          M[ ii + sx * jj ] += * col;
          //printf("Sp %4i %4i %.4f\n", ii, jj, a );
        }
      }
    }
  }
}


//----------------------------------------------------------------------
void MatrixSparse::checkDataConsistency() const
{
  assert( size > 0 );
  for(int jj = 0; jj < size; ++jj ) if ( Vrow[ jj ] )
    for (int ii = 0; Vrow[ jj ][ ii ] >= 0; ++ii )
      assert( Vrow[ jj ][ ii ] < size );
}


//----------------------------------------------------------------------
void MatrixSparse::printSparse(FILE * f) const
{
  for(int jj = 0; jj < size; ++jj ) 
    if ( Vrow[ jj ] )
      for (int ii = 0; Vrow[ jj ][ ii ] >= 0; ++ii )
        fprintf(f, "%i %i %16.8e\n", Vrow[ jj ][ ii ], jj, Vcol[ jj ][ ii ] );
}

//----------------------------------------------------------------------
bool MatrixSparse::nonZero() const
{
  for( int jj = 0; jj < size; ++jj ) if ( Vrow[ jj ] )
    for ( int ii = 0; Vrow[ jj ][ ii ] >= 0; ++ii )
      if ( Vcol[ jj ][ ii ] != 0 )
        return true;
  return false;
}

//----------------------------------------------------------------------
int MatrixSparse::nbNonZeroElements() const
{
  //all allocated elements are counted, even if the value is zero
  int cnt = 0;
  for( int jj = 0; jj < size; ++jj ) if ( Vrow[ jj ] )
    for ( int ii = 0; Vrow[ jj ][ ii ] >= 0; ++ii )
      ++cnt;
  return cnt;
}

//----------------------------------------------------------------------
char * MatrixSparse::shortDescription() const
{
  static char msg[16];
  snprintf(msg, sizeof(msg), "SP nnz %i", nbNonZeroElements() );
  return msg;
}

//========================================================================
//========================================================================

void MatrixSparse::multiplyVectorAdd( const real * X, real * Y ) const
{
  int ii, jj, kk;
    
  for( jj = 0; jj < size; ++jj )
    if ( Vrow[ jj ] )
      for ( ii = 0; ( kk = Vrow[ jj ][ ii ] ) >= 0; ++ii ) {
        Y[ kk ] += Vcol[ jj ][ ii ] * X[ jj ];
      }
}

#ifndef FAST_VECTMULT

//----------------------------------------------------------------------
void MatrixSparse::multiplyVectorAddIsotropic( const real * X, real * Y ) const
{
  int ii, jj, kk;
  real a;
    
  for( jj = 0; jj < size; ++jj ) 
    if ( Vrow[ jj ] )
      for ( ii = 0; ( kk = Vrow[ jj ][ ii ] ) >= 0; ++ii ) {
        a = Vcol[ jj ][ ii ];
        Y[ DIM*kk     ] += a * X[ DIM*kk     ];
#if ( DIM > 1 )
        Y[ DIM*kk + 1 ] += a * X[ DIM*kk + 1 ];
#endif
#if ( DIM > 2 )
        Y[ DIM*kk + 2 ] += a * X[ DIM*kk + 2 ]; 
#endif
      }
}


#else
//========================================================================
//                        FAST_VECTMULT
//========================================================================
#if ( DIM == 1 )

//----------------------------------------------------------------------
void MatrixSparse::multiplyVectorAddIsotropic( const real * X, real * Y ) const
{      
  multiplyVectorAdd( X, Y );
}
#endif

#if ( DIM == 2 )

//----------------------------------------------------------------------
void MatrixSparse::multiplyVectorAddIsotropic( const real * X, real * Y ) const
{      
 int kk, ll, jj, * row;

 for(jj = 0; jj < size; ++jj )
   if ( ( row = Vrow[ jj ] ) != 0 )
     {
       real * col = Vcol[ jj ];
       ll = DIM * jj;

       real X1 = X[ ll     ];
       real X2 = X[ ll + 1 ];
       
       while ( ( kk = DIM * ( * row ) ) >= 0 ) {
         Y[ kk   ] += (*col) * X1;
         Y[ kk+1 ] += (*col) * X2;
         
         ++row;
         ++col;
       }
     }
}
#endif

#if ( DIM == 3 )

//----------------------------------------------------------------------
void MatrixSparse::multiplyVectorAddIsotropic( const real * X, real * Y ) const
{      
 int kk, ll, jj, * row;

 for(jj = 0; jj < size; ++jj )
   if ( ( row = Vrow[ jj ] ) != 0 )
     {
       real * col = Vcol[ jj ];
       ll = DIM * jj;

       real X1 = X[ ll     ];
       real X2 = X[ ll + 1 ];
       real X3 = X[ ll + 2 ];

       while ( ( kk = DIM * ( * row ) ) >= 0 ) {
         Y[ kk   ] += (*col) * X1;
         Y[ kk+1 ] += (*col) * X2;
         Y[ kk+2 ] += (*col) * X3;
         
         ++row;
         ++col;
       }
     }
}
#endif

#endif

