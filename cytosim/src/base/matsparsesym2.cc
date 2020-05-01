//RCS: $Id: matsparsesym2.cc,v 2.10 2005/04/24 22:47:46 nedelec Exp $
#include <cstdio>

#include "matsparsesym2.h"
#include "cblas.h"
#include "main.h"
#include "smath.h"

#define FAST_VECTMULT

//----------------------------------------------------------------------
void MatrixSparseSymmetric2::allocate()
{ 
  size       = 0;
  allocated  = 0;
  Vcol       = 0;
  Vsize      = 0;
  Vallo      = 0;
  nextCol    = 0;
  firstCol   = 0;
}


//----------------------------------------------------------------------
void MatrixSparseSymmetric2::allocate( const int sz )
{
  size = sz;
  if ( size > allocated ) {

    Element **  Vcol_new   = new Element * [ size ];	  
    int     *  Vsize_new   = new int[ size ];
    int     *  Vallo_new   = new int[ size ];
    int     *  nextCol_new = new int[ size ];
    
    if ((Vcol_new == 0) || (Vsize_new == 0) || (Vallo_new == 0)) {
      fprintf(stderr, "MatrixSparseSymmetric2::allocate(int):: memory allocation failed\n");
      exit(1);
    }
    
    int ii = 0;
    if ( Vcol ) {
      for( ; ii < allocated; ++ii ) {
        Vcol_new[ ii ]    = Vcol[ ii ];
        Vsize_new[ ii ]   = Vsize[ ii ];
        Vallo_new[ ii ]   = Vallo[ ii ];
        nextCol_new[ ii ] = nextCol[ ii ];
      }
      delete[] Vcol;
      delete[] Vsize;
      delete[] Vallo;
      delete[] nextCol;
    }
    
    for( ; ii < size; ++ii ) {
      Vcol_new[ ii ]    = 0;
      Vsize_new[ ii ]   = 0;
      Vallo_new[ ii ]   = 0;
      nextCol_new[ ii ] = ii+1;
    }

    Vcol      = Vcol_new;
    Vsize     = Vsize_new;
    Vallo     = Vallo_new;
    nextCol   = nextCol_new;
    allocated = size;
  }
}

//----------------------------------------------------------------------
void MatrixSparseSymmetric2::deallocate()
{
  if ( Vcol ) {
    for (int ii = 0; ii < allocated; ++ii )
      if ( Vcol[ ii ] ) {
        delete[] Vcol[ ii ];
      };
    delete[] Vcol;
    delete[] Vsize;
    delete[] Vallo;
  }
}

//----------------------------------------------------------------------
void MatrixSparseSymmetric2::allocateColumn( const int jj, int sz )
{
  assert((jj >= 0) && (jj < size));
  assert( sz > 0 );
  //printf("new S-COL %i %i\n", jj, sz );

  if ( sz > Vallo[ jj ] ) {
    
    sz += SPCHUNCK;
    Element * Vcol_new  = new Element[ sz ];
  
    if (Vcol_new == 0) {
      fprintf(stderr, "MatrixSparseSymmetric2::allocateColumn(int):: memory allocation failed\n");
      exit(1);
    }
    
    if ( Vcol[ jj ] ) {
      
      //copy what is there
      for(int ii = 0; ii < Vallo[jj]; ++ii ) {
        Vcol_new[ ii ] = Vcol[ jj ][ ii ];
      }
      delete[] Vcol[ jj ];
      
    }  
    Vcol[ jj ]  = Vcol_new;
    Vallo[ jj ] = sz;
  }
}

//----------------------------------------------------------------------
real & MatrixSparseSymmetric2::operator()( int ii, int jj )
{
  //this allocate the position if necessary
  assert(( ii >= 0 ) && ( ii < size ));
  assert(( jj >= 0 ) && ( jj < size ));

  Element * col, * last;

  //we swap to get the order right
  if ( jj < ii ) {
    int tmp = ii;
    ii  = jj;
    jj  = tmp;
  }
  
  //check if the column is empty:
  if ( Vsize[ jj ] == 0 ) {
    allocateColumn( jj, 2 );
    col = Vcol[ jj ];
    //add the diagonal term first:
    col->indx = jj;
    col->val  = 0.;
    Vsize[ jj ] = 1;
    if ( ii == jj )
      return col->val;
    
    //add the requested term:
    ++col;
    col->indx = ii;
    col->val  = 0.;
    Vsize[ jj ] = 2;  
    return col->val;    
  }

  col = Vcol[ jj ];
  assert( col->indx == jj ); //the first term should be the diagonal

  //check if diagonal term is requested
  if ( ii == jj )
    return col->val;
  
  //TODO: optimize this search in MatrixSparse::operator()
  last = col + Vsize[ jj ];
  ++col; //we have checked the first term already, which is the diagonal
  while ( col < last ) {
    if ( col->indx == ii )
      return col->val;
    ++col;
  }
  
  //we will have to create/allocate a new Element
  if ( Vallo[jj] <= Vsize[jj] ) {
    allocateColumn( jj, Vsize[jj]+1 );
    col = Vcol[jj] + Vsize[jj];
  }
   
  assert( Vallo[ jj ] > Vsize[ jj ] );
  
  //add the requested term:
  col->indx = ii;
  col->val  = 0.;
  ++Vsize[ jj ];  
  return col->val;    
}


//----------------------------------------------------------------------
real * MatrixSparseSymmetric2::addr( int x, int y) const
{
  //this does not allocate the position:
  //may swap to get the order right
  int nx, ny;
  if ( y < x ) {
    nx = y;
    ny = x;
  } else {
    nx = x;
    ny = y;
  }
  
  for( int kk = 0; kk < Vsize[ ny ]; ++kk )
    if ( Vcol[ ny ][ kk ].indx == nx )
      return &( Vcol[ ny ][ kk ].val );
  return 0;
}
  

//----------------------------------------------------------------------
void MatrixSparseSymmetric2::setToZero()
{
  for(int ii = 0; ii < size; ++ii ) 
    Vsize[ ii ] = 0;
}


//----------------------------------------------------------------------
void MatrixSparseSymmetric2::scale( real a )
{
  for(int ii = 0; ii < size; ++ii )
    for(int jj = 0; jj < Vsize[ ii ]; ++jj )
        Vcol[ ii ][ jj ].val *= a;
}

//----------------------------------------------------------------------
// M <- M + the upper block contained in [x, x+sx, x, x+sx]
void MatrixSparseSymmetric2::addTriangularBlock(real* M, const int x, const int sx ) const
{
  assert( x + sx <= size );

  for(int jj = 0; jj < sx; ++jj )
    for(int kk = 0; kk < Vsize[ jj+x ]; ++kk ) {
      int ii = Vcol[ jj+x ][ kk ].indx - x;
      if (( ii >= 0 ) && ( ii < sx )) {
        assert( ii <= jj );	
        M[ ii + sx * jj ] += Vcol[ jj+x ][ kk ].val;
        //printf("Sp %4i %4i %.4f\n", ii, jj, a );
      }
    }      
}

//----------------------------------------------------------------------
// M <- M + the upper block contained in [x, x+sx, x, x+sx]
void MatrixSparseSymmetric2::addDiagonalBlock(real* M, const int x, const int sx ) const
{
  assert( x + sx <= size );
  
  for(int jj = 0; jj < sx; ++jj )
    for(int kk = 0 ; kk < Vsize[ jj+x ] ; ++kk ) {
        int ii = Vcol[ jj+x ][ kk ].indx - x;
        if (( ii >= 0 ) && ( ii < sx )) {
          assert( ii <= jj );	
          M[ ii + sx * jj ] += Vcol[ jj+x ][ kk ].val;
          if ( ii != jj )
            M[ jj + sx * ii ] += Vcol[ jj+x ][ kk ].val;
          //printf("Sp %4i %4i %.4f\n", ii, jj, a );
        }
      }
 }


//----------------------------------------------------------------------
void MatrixSparseSymmetric2::checkDataConsistency() const
{
  assert( size > 0 );
  for(int jj = 0; jj < size; ++jj )
    for(int kk = 0 ; kk < Vsize[ jj ] ; ++kk ) {
      assert( Vcol[ jj ][ kk ].indx < size );
      assert( Vcol[ jj ][ kk ].indx > jj );
    }
}


//----------------------------------------------------------------------
void MatrixSparseSymmetric2::printSparse(FILE * f) const
{
  for(int jj = 0; jj < size; ++jj )
    for(int kk = 0 ; kk < Vsize[ jj ] ; ++kk )
      fprintf(f, "%i %i %16.8e\n", Vcol[ jj ][ kk ].indx, jj, Vcol[ jj ][ kk ].val );
}


//----------------------------------------------------------------------
bool MatrixSparseSymmetric2::nonZero() const
{
  //check for any non-zero sparse term:
  for(int jj = 0; jj < size; ++jj )
    for(int kk = 0 ; kk < Vsize[ jj ] ; ++kk )
      if ( Vcol[ jj ][ kk ].val )
        return true;
      
  //if here, the matrix is empty
  return false;
}

//----------------------------------------------------------------------
int MatrixSparseSymmetric2::nbNonZeroElements() const
{
  //all allocated elements are counted, even if zero
  //the diagonal is not counted
  int cnt = 0;
  for( int jj = 0; jj < size; ++jj ) 
    cnt += Vsize[ jj ];
  return cnt;
}

//----------------------------------------------------------------------
char * MatrixSparseSymmetric2::shortDescription() const
{
  static char msg[16];
  snprintf(msg, sizeof(msg), "SPS2 nnz %i", nbNonZeroElements() );
  return msg;
}

//========================================================================
//========================================================================

int MatrixSparseSymmetric2::mat_comp( const void * a, const void * b )
{
  return ((Element*) a) -> indx - ((Element*)b) -> indx;  
}

//this optional optimization may accelerate the multiplications below
//by making memory access more ordered (check with profiling: no garanty)
void MatrixSparseSymmetric2::optimizeForMultiply()
{
  //update nextCol[], a pointer to the next non-empty column:
  firstCol = size;
  for(int jj = size-1; jj >= 0; --jj ) {
    nextCol[jj] = firstCol;
    if ( Vsize[jj] > 0 ) 
      firstCol = jj;
  }
  
  //the second optimization is not so clear, in terms of CPU gain:
  return;
  
  //we reorder the columns to facilitate memory access
  for(int jj = firstCol; jj < size; jj = nextCol[jj] )
    if ( Vsize[ jj ] > 1 ) {
      
      //printf("\nColumn %i:\n", jj);
      //for(int kk=0; kk<Vsize[jj]; ++kk) printf("%i:%.2f ", Vcol[jj][kk].indx, Vcol[jj][kk].val);
      
      //optimize by putting the diagonal term first
      //the diagonal term is usually not very far, so we do a simple search
      for(int kk = 0; kk < Vsize[ jj ]; ++kk)
        if ( Vcol[jj][kk].indx == jj ) {
          if ( kk > 0 ) {
            Element tmp  = Vcol[jj][0];
            Vcol[jj][0]  = Vcol[jj][kk];
            Vcol[jj][kk] = tmp;
          }
          break;
        }
      
      qsort( Vcol[jj]+1, Vsize[jj]-1, sizeof(Element), &mat_comp );
      
      //printf("\n");printf("\nColumn %i:\n", jj);
      //for(int kk=0; kk<Vsize[jj]; ++kk) printf("%i:%.2f ", Vcol[jj][kk].indx, Vcol[jj][kk].val);
    }
}


#ifndef FAST_VECTMULT

//----------------------------------------------------------------------
void MatrixSparseSymmetric2::multiplyVectorAdd( const real * X, real * Y ) const
{
  int ii;
  
  for(int jj = 0; jj < size; ++jj )
    for(int kk = 0 ; kk < Vsize[ jj ] ; ++kk ) {
      ii = Vcol[ jj ][ kk ].indx;
      Y[ ii ] += (Vcol[ jj ][ kk ].val) * X[ jj ];
      if ( ii != jj ) 
        Y[ jj ] += (Vcol[ jj ][ kk ].val) * X[ ii ];
    }
}

//----------------------------------------------------------------------
void MatrixSparseSymmetric2::multiplyVectorAddIsotropic( const real * X, real * Y ) const
{
  int ii;
    
  for(int jj = 0; jj < size; ++jj ) 
    for(int kk = 0 ; kk < Vsize[ jj ] ; ++kk ) {
      ii = Vcol[ jj ][ kk ].indx;
      if ( ii == jj )
        blas_xaxpy( DIM, Vcol[ jj ][ kk ].val, X + DIM*ii, 1, Y + DIM*ii, 1 );
      else {
        blas_xaxpy( DIM, Vcol[ jj ][ kk ].val, X + DIM*jj, 1, Y + DIM*ii, 1 );
        blas_xaxpy( DIM, Vcol[ jj ][ kk ].val, X + DIM*ii, 1, Y + DIM*jj, 1 );
      }
    }
}

#else

//----------------------------------------------------------------------
void MatrixSparseSymmetric2::multiplyVectorAdd( const real * X, real * Y ) const
{
  register int kk, ii, siz;
  register real val, X1, Y1, P1;
  Element * col;

  //we only visit the non-empty columns
  for( int jj = firstCol; jj < size; jj = nextCol[jj] ) {
    siz = Vsize[ jj ];
    assert( siz > 0 );
    
    col = Vcol[ jj ];
    X1 = X[ jj ];
    
    //the diagonal term should be first
    assert( col->indx == jj );
    val = col->val;
    ++col;
    kk = 1;
    Y1 = val * X1;
    
    while( kk < siz ) {
      val = col->val;
      ii  = col->indx;
      assert( ii != jj );
      P1 = val * X1;
      Y1 += val * X[ ii ];
      ++kk;
      ++col;
      Y[ ii ] += P1;
    }
    Y[ jj ] += Y1;
  }
}
//========================================================================
//                        FAST_VECTMULT
//========================================================================
#if ( DIM == 1 )

//----------------------------------------------------------------------
void MatrixSparseSymmetric2::multiplyVectorAddIsotropic( const real * X, real * Y ) const
{      
  multiplyVectorAdd( X, Y );
}
#endif

#if ( DIM == 2 )
/*
//----------------------------------------------------------------------
void MatrixSparseSymmetric2::multiplyVectorAddIsotropic( const real * X, real * Y ) const
{
  register int ll, ii;
  register real val;
  
  for(int jj = 0; jj < size; ++jj ) 
    if ( Vsize[ jj ] ) {

       ll = DIM * jj;
      
       for(int kk = 0 ; kk < Vsize[ jj ] ; ++kk ) {
         
         val = Vcol[ jj ][ kk ].val;
         ii  = DIM * Vcol[ jj ][ kk ].indx;
         
         Y[ ii   ] += val * X[ ll   ];
         Y[ ii+1 ] += val * X[ ll+1 ];
         
         if ( ii != ll ) {
           Y[ ll   ] += val * X[ ii   ];
           Y[ ll+1 ] += val * X[ ii+1 ];
         }
       }
     }
}
*/

//----------------------------------------------------------------------
void MatrixSparseSymmetric2::multiplyVectorAddIsotropic( const real * X, real * Y ) const
{
  register int kk, ii, siz;
  register real val, X1, X2, Y1, Y2, P1, P2;
  Element * col;
  
  //we only visit the non-empty columns
  for( int jj = firstCol; jj < size; jj = nextCol[jj] ) {
    int ll = DIM * jj;
    siz = Vsize[ jj ];
    assert( siz > 0 );

    col = Vcol[ jj ];

    X1 = X[ ll   ];
    X2 = X[ ll+1 ];
    
    //the diagonal term should be first
    assert( col->indx == jj );
    val = col->val;
    ++col;
    kk = 1;
    Y1 = val * X1;
    Y2 = val * X2;
    
    while( kk < siz ) {
      
      val = col->val;
      ii  = DIM * col->indx;
      
      assert( ii != ll );
      
      P1 = val * X1;
      P2 = val * X2;
      
      Y1 += val * X[ ii   ];
      Y2 += val * X[ ii+1 ];
      
      ++kk;
      ++col;
      
      Y[ ii   ] += P1;
      Y[ ii+1 ] += P2;
    }
    
    Y[ ll   ] += Y1;
    Y[ ll+1 ] += Y2;
  }
}

#endif

#if ( DIM == 3 )

//----------------------------------------------------------------------
void MatrixSparseSymmetric2::multiplyVectorAddIsotropic( const real * X, real * Y ) const
{
  register int kk, ii, siz;
  register real val, X1, X2, X3, Y1, Y2, Y3;
  register real P1, P2, P3;
  Element * col;
  
  for( int jj = firstCol; jj < size; jj = nextCol[jj] ) {
    int ll = DIM * jj;
    siz = Vsize[ jj ];
    assert( siz > 0 );
    
    col = Vcol[ jj ];
      
    X1 = X[ ll   ];
    X2 = X[ ll+1 ];
    X3 = X[ ll+2 ];
      
    //the diagonal term should be first
    assert( col->indx == jj );
    val = col->val;
    ++col;
    kk = 1;
    Y1 = val * X1;
    Y2 = val * X2;
    Y3 = val * X3;
    
    while( kk < siz ) {
      
      val = col->val;
      ii  = DIM * col->indx;
      
      assert( ii != ll );
      
      P1 = val * X1;
      P2 = val * X2;
      P3 = val * X3;
      
      Y1 += val * X[ ii   ];
      Y2 += val * X[ ii+1 ];
      Y3 += val * X[ ii+2 ];
      
      ++kk;
      ++col;
      
      Y[ ii   ] += P1;
      Y[ ii+1 ] += P2;
      Y[ ii+2 ] += P3;
    }
    
    Y[ ll   ] += Y1;
    Y[ ll+1 ] += Y2;
    Y[ ll+2 ] += Y3;
  }
}
#endif

#endif

