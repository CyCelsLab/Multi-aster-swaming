//RCS: $Id: matsparsesym.cc,v 2.2 2005/03/01 13:56:48 nedelec Exp $
#include <cstdio>

#include "matsparsesym.h"
#include "cblas.h"
#include "main.h"
#include "smath.h"

#define FAST_VECTMULT

#define VAL_NOT_USED -1
#define END_OF_COLUMN -2

//----------------------------------------------------------------------
void MatrixSparseSymmetric::allocate()
{ 
  size      = 0;
  allocated = 0;
  Vcol      = 0;
  Vrow      = 0;
}


//----------------------------------------------------------------------
void MatrixSparseSymmetric::allocate( const int sz )
{
  size = sz;
  if ( size > allocated ) {

    real ** Vcol_new = new real * [ size ];	  
    int  ** Vrow_new = new int * [ size ];	 
    
    if ((Vcol_new == 0) || (Vrow_new == 0)) {
      fprintf(stderr, "MatrixSparseSymmetric::allocate(int):: memory allocation failed\n");
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
void MatrixSparseSymmetric::deallocate()
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
void MatrixSparseSymmetric::allocateColumn( const int jj, const int sz )
{
  assert((jj >= 0) && (jj < size));
  assert( sz > 0 );
  //printf("new S-COL %i %i\n", jj, sz );

  real * Vcol_new  = new real[ sz ];
  int  * Vrow_new  = new  int[ sz ];
  
  if ((Vcol_new == 0) || (Vrow_new == 0)) {
    fprintf(stderr, "MatrixSparseSymmetric::allocateColumn(int):: memory allocation failed\n");
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
real & MatrixSparseSymmetric::operator()( int x, int y )
  //allocate the position if necessary:
{
  assert( ( x >= 0 ) && ( x < size ) );
  assert( ( y >= 0 ) && ( y < size ) );

  //may swap to get the order right
  int nx, ny;
  if ( y < x ) {
    nx = y;
    ny = x;
  } else {
    nx = x;
    ny = y;
  }
  
  if ( Vrow[ ny ] ) {
    
    int ii = 0;
    for( ; Vrow[ ny ][ ii ] >= 0; ++ii )
      if ( Vrow[ ny ][ ii ] == nx )
        return Vcol[ ny ][ ii ];
    
    if( Vrow[ ny ][ ii ] == END_OF_COLUMN )
      allocateColumn( ny, ii + 1 + SPCHUNCK );
    assert( Vrow[ ny ][ ii ] == VAL_NOT_USED );
    Vrow[ ny ][ ii ] = nx;
    Vcol[ ny ][ ii ] = 0;
    //printf("allo. %3i %3i\n", nx, ny );
    return Vcol[ ny ][ ii ];
  }
  
  allocateColumn( ny, SPCHUNCK );
  //printf("allo. %3i %3i\n", nx, ny );
  assert( Vrow[ ny ][ 0 ] == VAL_NOT_USED );

  //put the diagonal term first:
  Vrow[ ny ][ 0 ] = ny;
  Vcol[ ny ][ 0 ] = 0;
  if ( nx == ny )
    return Vcol[ ny ][ 0 ];
  
  Vrow[ ny ][ 1 ] = nx;
  Vcol[ ny ][ 1 ] = 0;
  return Vcol[ ny ][ 1 ];
}


//----------------------------------------------------------------------
real * MatrixSparseSymmetric::addr( int x, int y) const
  //does not allocate the position:
{
  //may swap to get the order right
  int nx, ny;
  if ( y < x ) {
    nx = y;
    ny = x;
  } else {
    nx = x;
    ny = y;
  }
  
  int * row = Vrow[ ny ];
  if ( row )
    for( ; *row >= 0 ; ++row )
      if ( *row == nx )
	return &Vcol[ ny ][ row - Vrow[ ny ] ];
  return 0;
}
  

//----------------------------------------------------------------------
void MatrixSparseSymmetric::setToZero()
{
  for(int ii = 0; ii < size; ++ii ) 
    if ( Vrow[ ii ] ) 
      for(int jj = 0; Vrow[ ii ][ jj ] != END_OF_COLUMN; ++jj )
        Vrow[ ii ][ jj ] = VAL_NOT_USED;
}


//----------------------------------------------------------------------
void MatrixSparseSymmetric::scale( real a )
{
  for(int ii = 0; ii < size; ++ii ) if ( Vrow[ ii ] ) 
    for(int jj = 0; Vrow[ ii ][ jj ] >= 0; ++jj )
      Vcol[ ii ][ jj ] *= a;
}

//----------------------------------------------------------------------
void MatrixSparseSymmetric::addTriangularBlock(real* M, const int x, const int sx ) const
  // M <- M + the upper block contained in [x, x+sx, x, x+sx]
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
void MatrixSparseSymmetric::addDiagonalBlock(real* M, const int x, const int sx ) const
{
  assert( x + sx <= size );

  for(int jj = 0; jj < sx; ++jj ) {
    
    int *row = Vrow[ jj + x ];
    if ( row != 0 ) {
      
      real * col = Vcol[ jj + x ];
      for( ; *row >= 0 ; ++row, ++col) {
        int ii = *row - x;
        if ( ( ii >= 0 ) && ( ii < sx ) ) {
          
          assert( ii <= jj );	
          M[ ii + sx * jj ] += * col;
          if ( ii != jj )
            M[ jj + sx * ii ] += * col;
          //printf("Sp %4i %4i %.4f\n", ii, jj, a );
        }
      }
    }
  }
}


//----------------------------------------------------------------------
void MatrixSparseSymmetric::checkDataConsistency() const
{
  assert( size > 0 );
  for(int jj = 0; jj < size; ++jj ) if ( Vrow[ jj ] )
    for (int ii = 0; Vrow[ jj ][ ii ] >= 0; ++ii )
	assert( Vrow[ jj ][ ii ] < size );
}


//----------------------------------------------------------------------
void MatrixSparseSymmetric::printSparse(FILE * f) const
{
  for(int jj = 0; jj < size; ++jj ) if ( Vrow[ jj ] )
    for (int ii = 0; Vrow[ jj ][ ii ] >= 0; ++ii )
      fprintf(f, "%i %i %16.8e\n", Vrow[ jj ][ ii ], jj, Vcol[ jj ][ ii ] );
}


//----------------------------------------------------------------------
bool MatrixSparseSymmetric::nonZero() const
{
  //check for any non-zero sparse term:
  for( int jj = 0; jj < size; ++jj ) if ( Vrow[ jj ] )
    for ( int ii = 0; Vrow[ jj ][ ii ] >= 0; ++ii )
      if ( 0 != Vcol[ jj ][ ii ] )
        return true;
  
  //if here, the matrix is empty
  return false;
}

//----------------------------------------------------------------------
int MatrixSparseSymmetric::nbNonZeroElements() const
{
  //all allocated elements are counted, even if zero
  //the diagonal is not counted
  int cnt = 0;
  for( int jj = 0; jj < size; ++jj ) if ( Vrow[ jj ] )
    for ( int ii = 0; Vrow[ jj ][ ii ] >= 0; ++ii )
      ++cnt;
  return cnt;
}

//----------------------------------------------------------------------
char * MatrixSparseSymmetric::shortDescription() const
{
  static char msg[16];
  snprintf(msg, sizeof(msg), "SPS nnz %i", nbNonZeroElements() );
  return msg;
}

//========================================================================
//========================================================================


void MatrixSparseSymmetric::multiplyVectorAdd( const real * X, real * Y ) const
{
  int ii, jj, kk;
    
  for( jj = 0; jj < size; ++jj )
    if ( Vrow[ jj ] )
      for ( ii = 0; ( kk = Vrow[ jj ][ ii ] ) >= 0; ++ii ) {
        Y[ kk ] += Vcol[ jj ][ ii ] * X[ jj ];
        if ( kk != jj ) 
          Y[ jj ] += Vcol[ jj ][ ii ] * X[ kk ];
      }
}

#ifndef FAST_VECTMULT

//----------------------------------------------------------------------
void MatrixSparseSymmetric::multiplyVectorAddIsotropic( const real * X, real * Y ) const
{
  int ii, jj, kk;
    
  for( jj = 0; jj < size; ++jj ) 
    if ( Vrow[ jj ] )
      for ( ii = 0; ( kk = Vrow[ jj ][ ii ] ) >= 0; ++ii )
        if ( kk == jj ) 
          blas_xaxpy( DIM, Vcol[ jj ][ ii ], X + DIM*kk, 1, Y + DIM*kk, 1 );
        else {
          blas_xaxpy( DIM, Vcol[ jj ][ ii ], X + DIM*jj, 1, Y + DIM*kk, 1 );
          blas_xaxpy( DIM, Vcol[ jj ][ ii ], X + DIM*kk, 1, Y + DIM*jj, 1 );
        }
}

#else
//========================================================================
//                        FAST_VECTMULT
//========================================================================
#if ( DIM == 1 )

//----------------------------------------------------------------------
void MatrixSparseSymmetric::multiplyVectorAddIsotropic( const real * X, real * Y ) const
{      
  multiplyVectorAdd( X, Y );
}
#endif

#if ( DIM == 2 )

//----------------------------------------------------------------------
void MatrixSparseSymmetric::multiplyVectorAddIsotropic( const real * X, real * Y ) const
{      
 int kk, ll, jj, * row;

 for(jj = 0; jj < size; ++jj )
   if ( ( row = Vrow[ jj ] ) != 0 )
     {
       real * col = Vcol[ jj ];
       ll = DIM * jj;

       real X1 = X[ ll     ];
       real X2 = X[ ll + 1 ];

       real & Y1 = Y[ ll     ];
       real & Y2 = Y[ ll + 1 ];
       
       while ( ( kk = DIM * ( * row ) ) >= 0 ) {
         Y[ kk   ] += (*col) * X1;
         Y[ kk+1 ] += (*col) * X2;
         
         if ( kk != ll ) {
           Y1 += (*col) * X[ kk   ];
           Y2 += (*col) * X[ kk+1 ];
         }
         
         ++row;
         ++col;
       }
     }
}
#endif

#if ( DIM == 3 )

//----------------------------------------------------------------------
void MatrixSparseSymmetric::multiplyVectorAddIsotropic( const real * X, real * Y ) const
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

       real & Y1 = Y[ ll     ];
       real & Y2 = Y[ ll + 1 ];
       real & Y3 = Y[ ll + 2 ];
       
       while ( ( kk = DIM * ( * row ) ) >= 0 ) {
         Y[ kk   ] += (*col) * X1;
         Y[ kk+1 ] += (*col) * X2;
         Y[ kk+2 ] += (*col) * X3;
         
         if ( kk != ll ) {
           Y1 += (*col) * X[ kk   ];
           Y2 += (*col) * X[ kk+1 ];
           Y3 += (*col) * X[ kk+2 ];
         }
         
         ++row;
         ++col;
       }
     }
}
#endif

#endif

