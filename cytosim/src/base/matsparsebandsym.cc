//RCS: $Id: matsparsebandsym.cc,v 2.3 2005/03/01 13:56:44 nedelec Exp $

#include "main.h"
#include "matsparsebandsym.h"
#include "cblas.h"
#include "smath.h"

#define FAST_VECTMULT

#define VAL_NOT_USED -1
#define END_OF_COLUMN -2

//----------------------------------------------------------------------
void MatrixSparseBandSymmetric::allocate() { 
  size      = 0;
  allocated = 0;
  diag      = 1;
  allocated_diag = 0;
  val       = 0;
  Vcol      = 0;
  Vrow      = 0;
}


//----------------------------------------------------------------------
void MatrixSparseBandSymmetric::allocate( const int sz, const int di ) {
  assert( sz >= 0 );
  size = sz;
  
  if ( size > allocated ) {
    
    if ( val ) delete[] val;
    val = 0;
    
    real ** Vcol_new = new real * [ size ];	  
    int  ** Vrow_new = new int * [ size ];	 
    
    if ((Vcol_new == 0) || (Vrow_new == 0)) {
      fprintf(stderr, "MatrixSparseBandSymmetric::allocate(int):: memory allocation failed\n");
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
  
  assert(( 0 <= di ) && ( di < 20 ));
  if ( ( val == 0 ) || ( ( di + 1 ) * size > allocated_diag ) ) {
    if ( val )  delete[] val;
    allocated_diag = ( di + 1 ) * allocated;
    val = new real[ allocated_diag ];	
  }
  diag = di;
}


//----------------------------------------------------------------------
void MatrixSparseBandSymmetric::deallocate()
{
  if ( val ) delete[] val;
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
void MatrixSparseBandSymmetric::allocateColumn( const int jj, const int sz )
{
  assert( ( 0 <= jj ) && ( jj < size ));
  assert( sz > 0 );
  //printf("new S-COL %i %i\n", jj, sz );
  
  real * Vcol_new  = new real[ sz ];
  int  * Vrow_new  = new  int[ sz ];
  
  if ((Vcol_new == 0) || (Vrow_new == 0)) {
    fprintf(stderr, "MatrixSparseBandSymmetric::allocateColumn(int):: memory allocation failed\n");
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
//allocate the position if necessary:
real & MatrixSparseBandSymmetric::operator()( int x, int y )
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

  if ( ny <= nx+diag ) 
    return val[ nx + diag * (ny+1) ];
  
  if ( Vrow[ ny ] ) {
    
    int ii = 0;
    for( ii = 0; Vrow[ ny ][ ii ] >= 0; ++ii )
      if ( Vrow[ ny ][ ii ] == nx ) 
        return Vcol[ ny ][ ii ];
    
    if( Vrow[ ny ][ ii ] == END_OF_COLUMN ) 
      allocateColumn( ny, ii + 1 + SPCHUNCK );
    
    assert( Vcol[ ny ] != 0 );
    assert( Vrow[ ny ] != 0 );
    assert( Vrow[ ny ][ ii ] == VAL_NOT_USED );
    Vrow[ ny ][ ii ] = nx;
    Vcol[ ny ][ ii ] = 0;
    return Vcol[ ny ][ ii ];
  }
  
  allocateColumn( ny, SPCHUNCK );
  assert( Vcol[ ny ] != 0 );
  assert( Vrow[ ny ] != 0 );
  assert( Vrow[ ny ][ 0 ] == VAL_NOT_USED );
  Vrow[ ny ][ 0 ] = nx;
  Vcol[ ny ][ 0 ] = 0;
  return Vcol[ ny ][ 0 ];
}


//----------------------------------------------------------------------
real * MatrixSparseBandSymmetric::addr( int x, int y) const 
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
  
  if ( ny <= nx+diag ) 
    return &( val[ nx + diag * (ny+1) ] );
  
  if ( Vrow[ ny ] )
    for( int ii = 0; Vrow[ ny ][ ii ] >= 0; ++ii )
      if ( Vrow[ ny ][ ii ] == nx )
        return &( Vcol[ ny ][ ii ] );
  
  return 0;
}


//----------------------------------------------------------------------
void MatrixSparseBandSymmetric::setToZero() {
  int ii;
  for( ii = 0; ii < (diag+1)*size; ++ii ) 
    val[ ii ] = 0;
  
  for( ii = 0; ii < size; ++ii ) if ( Vrow[ ii ] ) 
    for(int jj = 0; Vrow[ ii ][ jj ] != END_OF_COLUMN; ++jj )
      Vrow[ ii ][ jj ] = VAL_NOT_USED;
}

//----------------------------------------------------------------------
void MatrixSparseBandSymmetric::scale( const real a ) {
  int ii;
  for( ii = 0; ii < (diag+1)*size; ++ii ) 
    val[ ii ] *= a;
  
  for( ii = 0; ii < size; ++ii ) if ( Vrow[ ii ] ) 
    for(int jj = 0; Vrow[ ii ][ jj ] >= 0; ++jj )
      Vcol[ ii ][ jj ] *= a;
}

//----------------------------------------------------------------------
// M <- the upper block contained in [x, x+sx, x, x+sx]
void MatrixSparseBandSymmetric::addDiagonalBlock(real* M, const int x, const int sx ) const
{
  assert( x + sx <= size );
  int ii, jj;
  
  //set the diagonal, colum by column :
  //set the off-diagonals, line by line :
  for( ii = 0; ii < sx ; ++ii ) {
    
    M[ ii + sx * ii ] += val[ ii + x + diag * ( ii + x + 1 ) ];
    
    for( jj = ii+1; (jj <= ii + diag) && (jj < sx) ; ++jj )
      {
      M[ ii + sx * jj ] += val[ ii + x + diag * ( jj + x + 1 ) ];
      M[ jj + sx * ii ] += val[ ii + x + diag * ( jj + x + 1 ) ];
      }
  }
  
  //set the sparse values :
  int * row;
  for( jj = 0; jj < sx; ++jj )
    {
    row = Vrow[ jj + x ];
    if ( row != 0 ) {
      real * col = Vcol[ jj + x ];
      for( ; *row >= 0 ; ++row, ++col)
        {
	    ii = *row - x;
	    if ( ( ii >= 0 ) && ( ii <  sx ) )
	      {
          assert( ii < jj );
          M[ ii + sx * jj ] += * col;
          M[ jj + sx * ii ] += * col;
          //printf("Sp %4i %4i %.4f\n", ii, jj, a );
	      }
        }
    }
    }
}

//----------------------------------------------------------------------
// M <- the upper block contained in [x, x+sx, x, x+sx]
void MatrixSparseBandSymmetric::addTriangularBlock(real* M, const int x, const int sx ) const
{
  assert( x + sx <= size );
  int ii, jj;
  
  //set the diagonal, colum by column :
  //set the off-diagonals, line by line :
  for( ii = 0; ii < sx ; ++ii )
    for( jj = ii; (jj <= ii + diag) && (jj < sx); ++jj )
      M[ ii + sx * jj ] += val[ ii + x + diag * ( jj + x + 1 ) ];
  
  //set the sparse values :
  int * row;
  
  for( jj = 0; jj < sx; ++jj ) {
    row = Vrow[ jj + x ];
    if ( row != 0 ) {
      real * col = Vcol[ jj + x ];
      for( ; * row >= 0 ; ++row, ++col) {
        int ii = *row - x;
        if ( ( ii >= 0 ) && ( ii <  sx ) ) {
          assert( ii < jj );
          M[ ii + sx * jj ] += * col;
          //printf("Sp %4i %4i %.4f\n", ii, jj, a );
        }
	    }
    }
  }
}



//----------------------------------------------------------------------
void MatrixSparseBandSymmetric::checkDataConsistency() const
{
  assert( size >= 0 );
  for(int jj = 0; jj < size; ++jj ) if ( Vrow[ jj ] )
    for (int ii = 0; Vrow[ jj ][ ii ] >= 0; ++ii )
      assert( Vrow[ jj ][ ii ] < size );
}


//----------------------------------------------------------------------
void MatrixSparseBandSymmetric::printSparse(FILE * f) const
{
  int ii, jj;
  //print the diagonal, colum by column :
  for( ii = 0; ii < size ; ++ii )
    for( jj = ii; (jj <= ii + diag) && (jj < size); ++jj )
      fprintf(f, "%i %i %.8e\n", ii, jj, val[ ii + diag * ( jj + 1 ) ] );
  
  for( jj = 0; jj < size; ++jj ) if ( Vrow[ jj ] )
    for ( ii = 0; Vrow[ jj ][ ii ] >= 0; ++ii ) {
      assert( Vrow[ jj ][ ii ] >= 0 );
      assert( Vrow[ jj ][ ii ] < size );
      fprintf(f, "%i %i %.8e\n", Vrow[ jj ][ ii ], jj, Vcol[ jj ][ ii ] );
    }
}

//----------------------------------------------------------------------
bool MatrixSparseBandSymmetric::nonZero() const
{
  //check the diagonal terms:
  for( int ii = 0; ii < (diag+1)*size; ++ii ) 
    if ( 0 != val[ ii ] )
      return true;
  
  //check the sparse terms:
  for( int jj = 0; jj < size; ++jj ) if ( Vrow[ jj ] )
    for ( int ii = 0; Vrow[ jj ][ ii ] >= 0; ++ii )
      if ( Vcol[ jj ][ ii ] != 0 )
        return true;
  
  //empty matrix:
  return false;
}

//----------------------------------------------------------------------
int MatrixSparseBandSymmetric::nbNonZeroElements() const
{
  //FIX: the diagonal elements are not counted here
  int cnt = 0;
  for( int jj = 0; jj < size; ++jj ) if ( Vrow[ jj ] )
    for ( int ii = 0; Vrow[ jj ][ ii ] >= 0; ++ii )
      if ( Vcol[ jj ][ ii ] != 0 )
        ++cnt;
  return cnt;
}

//----------------------------------------------------------------------
char * MatrixSparseBandSymmetric::shortDescription()  const
{
  static char msg[16];
  snprintf(msg, sizeof(msg), "SPBS%i nnz %i", diag, nbNonZeroElements() );
  return msg;
}


//========================================================================
//========================================================================
//========================================================================


void MatrixSparseBandSymmetric::multiplyVectorAdd( const real * X, real * Y )  const
{
  
  blas_xsbmv('U', size, diag, 1.0, val, diag+1, X, 1, 1.0, Y, 1);
  
  int kk;
  for(int jj = 0; jj < size; ++jj ) if ( Vrow[ jj ] )
    for (int ii = 0; ( kk = Vrow[ jj ][ ii ] ) >= 0; ++ii ) {
      assert( jj < size );
      assert( kk < size );
      assert( kk != jj );
      Y[ kk ] += Vcol[ jj ][ ii ] * X[ jj ];
      Y[ jj ] += Vcol[ jj ][ ii ] * X[ kk ];
    }
}

#ifndef FAST_VECTMULT

//----------------------------------------------------------------------
void MatrixSparseBandSymmetric::multiplyVectorAddIsotropic( const real * X, real * Y ) const 
{
  
  for(int d = 0; d < DIM ; ++d )
    blas_xsbmv('U', size, diag, 1.0, val, diag+1, X+d, DIM, 1.0, Y+d, DIM);
  
  int kk;
  for(int jj = 0; jj < size; ++jj ) if ( Vrow[ jj ] )
    for (int ii = 0; ( kk = Vrow[ jj ][ ii ] ) >= 0; ++ii ) {
      assert( jj < size );
      assert( kk < size );
      assert( kk != jj );
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
void MatrixSparseBandSymmetric::multiplyVectorAddIsotropic( const real * X, real * Y ) const
{      
  multiplyVectorAdd( X, Y );
}

#endif


#if ( DIM == 2 )

//----------------------------------------------------------------------
void MatrixSparseBandSymmetric::multiplyVectorAddIsotropic( const real * X, real * Y ) const
{
  
  int kk, * row;
  real * col;
  real * Y1, * Y2, X1, X2;
  
  assert( size >= 0 );
  
  for(int d = 0; d < DIM ; ++d )
    blas_xsbmv('U', size, diag, 1.0, val, diag+1, X+d, DIM, 1.0, Y+d, DIM);
  
  for(int jj = 0, ll = 0; jj < size; ++jj, ll += DIM )
    if ( ( row = Vrow[ jj ] ) != 0 ) {
      col = Vcol[ jj ];
      
      X1 = X[ ll   ];
      X2 = X[ ll+1 ];
      
      Y1 = &Y[ ll   ];
      Y2 = &Y[ ll+1 ];
      
      while ( ( kk = DIM * ( * row ) ) >= 0 ) {
        assert( kk < DIM * size );
        assert( kk != ll );
        
        Y[ kk   ]   += (*col) * X1;
        Y[ kk+1 ]   += (*col) * X2;
        
        *Y1    += (*col) * X[ kk   ];
        *Y2    += (*col) * X[ kk+1 ];
        
        ++row;
        ++col;
      }
    }
}
#endif

#if ( DIM == 3 )

//----------------------------------------------------------------------
void MatrixSparseBandSymmetric::multiplyVectorAddIsotropic( const real * X, real * Y ) const
{      
 
  int kk, * row;
  real * col;
  real * Y1, * Y2, *Y3, X1, X2, X3;
  
  assert( size >= 0 );
  
  for(int d = 0; d < DIM ; ++d )
    blas_xsbmv('U', size, diag, 1.0, val, diag+1, X+d, DIM, 1.0, Y+d, DIM);
  
  for(int jj = 0, ll = 0; jj < size; ++jj, ll += DIM )
    if ( ( row = Vrow[ jj ] ) != 0 ) {
      col = Vcol[ jj ];

      X1 = X[ ll   ];
      X2 = X[ ll+1 ];
      X3 = X[ ll+2 ];
      
      Y1 = &Y[ ll   ];
      Y2 = &Y[ ll+1 ];
      Y3 = &Y[ ll+2 ];
      
      while ( ( kk = DIM * ( * row ) ) >= 0 ) {
        assert( kk < DIM * size );
        assert( kk != ll );
        
        Y[ kk   ]   += (*col) * X1;
        Y[ kk+1 ]   += (*col) * X2;
        Y[ kk+2 ]   += (*col) * X3;
        
        *Y1    += (*col) * X[ kk   ];
        *Y2    += (*col) * X[ kk+1 ];
        *Y3    += (*col) * X[ kk+2 ];
        
        ++row;
        ++col;
      }
    }
}
#endif

#endif

